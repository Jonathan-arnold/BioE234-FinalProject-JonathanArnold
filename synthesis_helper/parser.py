"""Parsers for MetaCyc input files."""

from __future__ import annotations

import re
import sys
from pathlib import Path

from synthesis_helper.models import Chemical, Reaction


def _strip_proton(inchi: str) -> str:
    """Strip the /p ionization layer only."""
    return re.sub(r"/p[+-]\d+", "", inchi)


def _strip_stereo(inchi: str) -> str:
    """Strip /p plus /t /m /s stereo layers for loose matching."""
    inchi = re.sub(r"/p[+-]\d+", "", inchi)
    inchi = re.sub(r"/t[^/]+", "", inchi)
    inchi = re.sub(r"/m\d+", "", inchi)
    inchi = re.sub(r"/s\d+", "", inchi)
    return inchi


def _read_lines(filepath: str | Path) -> list[str]:
    """Read a file, handling both \\n and \\r line endings."""
    text = Path(filepath).read_text(errors="replace")
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return [line.strip() for line in text.split("\n")]


def parse_chemicals(filepath: str | Path) -> dict[int, Chemical]:
    """Parse good_chems.txt → {id: Chemical}.

    Expected format (tab-separated, with header):
        id  name  inchi  smiles
    """
    chemicals: dict[int, Chemical] = {}
    for line in _read_lines(filepath):
        if not line or line.startswith("id"):
            continue
        parts = line.split("\t")
        chem_id = int(parts[0])
        name = parts[1] if len(parts) > 1 else ""
        inchi = parts[2].strip('"') if len(parts) > 2 else ""
        smiles = parts[3] if len(parts) > 3 else ""
        chemicals[chem_id] = Chemical(id=chem_id, name=name, inchi=inchi, smiles=smiles)
    return chemicals


def parse_reactions(
    filepath: str | Path, chemicals: dict[int, Chemical]
) -> list[Reaction]:
    """Parse good_reactions.txt → list[Reaction].

    Expected format (tab-separated, with header):
        rxnid  ecnum  substrates(space-separated ids)  products(space-separated ids)
    """
    reactions: list[Reaction] = []
    for line in _read_lines(filepath):
        if not line or line.startswith("rxnid"):
            continue
        parts = line.split("\t")
        rxn_id = int(parts[0])
        ecnum = parts[1] if len(parts) > 1 else ""
        substrate_ids = (
            [int(x) for x in parts[2].split()] if len(parts) > 2 and parts[2].strip() else []
        )
        product_ids = (
            [int(x) for x in parts[3].split()] if len(parts) > 3 and parts[3].strip() else []
        )
        substrates = frozenset(
            chemicals[sid] for sid in substrate_ids if sid in chemicals
        )
        products = frozenset(
            chemicals[pid] for pid in product_ids if pid in chemicals
        )
        reactions.append(
            Reaction(id=rxn_id, substrates=substrates, products=products, ecnum=ecnum)
        )
    return reactions


def parse_metabolite_list(
    filepath: str | Path,
    chemicals: dict[int, Chemical],
    descriptor_filter: set[str] | None = None,
) -> set[Chemical]:
    """Parse a metabolite list file (minimal_metabolites.txt or ubiquitous_metabolites.txt).

    Expected format (tab-separated, with header, possibly \\r line endings):
        name  inchi  descriptor

    When ``descriptor_filter`` is provided, only rows whose descriptor is in
    the set are matched (empty/missing descriptors are dropped). ``None``
    (default) keeps every row — preserving the baseline behaviour.

    Matching strategy (in order):
      1. Strict InChI match after stripping the /p proton layer.
      2. Loose InChI match also stripping /t /m /s stereo layers — adds ALL chemicals
         with that connectivity so both stereoisomers reach shell 0 (handles SAM-like
         cases where MetaCyc encodes the same cofactor with inconsistent stereo).
      3. Name match (case-insensitive).
    """
    # Build two lookups: strict (/p stripped) and loose (/p+stereo stripped)
    strict_lookup: dict[str, Chemical] = {}
    loose_lookup: dict[str, list[Chemical]] = {}
    name_to_chem: dict[str, Chemical] = {}

    for chem in chemicals.values():
        raw = chem.inchi.strip('"')
        strict_key = _strip_proton(raw)
        loose_key = _strip_stereo(raw)
        if strict_key:
            strict_lookup[strict_key] = chem
        if loose_key:
            loose_lookup.setdefault(loose_key, []).append(chem)
        if chem.name:
            name_to_chem[chem.name.lower()] = chem

    metabolites: set[Chemical] = set()
    unmatched: list[str] = []

    for line in _read_lines(filepath):
        if not line or line.startswith("name"):
            continue
        parts = line.split("\t")
        name = parts[0] if len(parts) > 0 else ""
        raw_inchi = parts[1].strip('"') if len(parts) > 1 else ""
        descriptor = parts[2].strip() if len(parts) > 2 else ""

        if descriptor_filter is not None and descriptor not in descriptor_filter:
            continue

        matched = False

        # Pass 1: strict match (ionization-normalized)
        strict_key = _strip_proton(raw_inchi)
        if strict_key and strict_key in strict_lookup:
            metabolites.add(strict_lookup[strict_key])
            matched = True

        # Pass 2: loose match (stereo-normalized) — always runs so all stereoisomers
        # of a cofactor reach shell 0, not just the one that happened to match strictly.
        loose_key = _strip_stereo(raw_inchi)
        if loose_key and loose_key in loose_lookup:
            for chem in loose_lookup[loose_key]:
                metabolites.add(chem)
            matched = True

        # Pass 3: name match
        if not matched and name.lower() in name_to_chem:
            metabolites.add(name_to_chem[name.lower()])
            matched = True

        if not matched:
            unmatched.append(name)

    if unmatched:
        print(
            f"  Warning: {len(unmatched)} metabolites not matched to chemicals: "
            f"{unmatched[:5]}{'...' if len(unmatched) > 5 else ''}",
            file=sys.stderr,
        )

    return metabolites


def parse_reaction_organisms(filepath: str | Path) -> dict[int, list[str]]:
    """Parse the EnzymeMap organism sidecar → {rxnid: [organism, ...]}.

    Expected format (tab-separated, with header):
        rxnid  organisms       (comma-separated organism strings)
    """
    out: dict[int, list[str]] = {}
    for line in _read_lines(filepath):
        if not line or line.startswith("rxnid"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        rxn_id = int(parts[0])
        out[rxn_id] = [o for o in parts[1].split(",") if o]
    return out
