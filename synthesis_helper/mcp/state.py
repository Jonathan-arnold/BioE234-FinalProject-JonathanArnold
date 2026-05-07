"""Lazy-loaded singleton holding parsed data + the baseline HyperGraph."""

from __future__ import annotations

import os
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from synthesis_helper.models import Chemical, HyperGraph, Reaction
from synthesis_helper.parser import (
    parse_chemicals,
    parse_metabolite_list,
    parse_reactions,
)
from synthesis_helper.synthesize import synthesize


DEFAULT_DATA_DIR = Path(__file__).resolve().parents[2] / "data"


CURRENCY_DESCRIPTORS: set[str] = {"inorganics", "cofactor", "nucleotide"}

# Canonical cofactor labels → alternative names to recognise. Order matters:
# NADPH / NADH before NADP+ / NAD+ so reduced forms win when a stereoisomer
# overlaps. Each chemical is assigned to the FIRST matching label only.
#
# Matching rule: `chem.name.lower()` EQUALS one of the alternatives, OR
# starts with "<alt> (" for variant suffixes like "NAD+ (0)" or "SAM (27-)".
# Substring matching is intentionally avoided: "atp" would false-positive on
# "dATP" and "coa" would false-positive on "acetyl-coa".
_COFACTOR_MATCHERS: list[tuple[str, frozenset[str]]] = [
    ("NADPH", frozenset({"nadph", "beta-nadph"})),
    ("NADH",  frozenset({"nadh", "beta-nadh"})),
    ("NADP+", frozenset({"nadp+", "nadp", "nadp(+)"})),
    ("NAD+",  frozenset({"nad+", "nad", "nad(+)"})),
    ("ATP",   frozenset({"atp"})),
    ("CoA",   frozenset({"coa", "coenzyme a", "coenzyme-a"})),
    ("SAM",   frozenset({"sam", "s-adenosyl-l-methionine", "s-adenosylmethionine"})),
]
COFACTOR_KEYS: tuple[str, ...] = tuple(k for k, _ in _COFACTOR_MATCHERS)


def _name_matches_cofactor(lname: str, alts: frozenset[str]) -> bool:
    if lname in alts:
        return True
    for alt in alts:
        if lname.startswith(alt + " ("):
            return True
    return False


# Canonical toxic / reactive intermediates known to stall heterologous
# pathways in E. coli — reactive aldehydes (genotoxic / crosslinking),
# phenolics prone to auto-oxidation / quinone formation, Michael acceptors,
# and ROS-generating quinones. A pathway touching any of these flags a
# downstream production risk regardless of whether the individual step is
# catalytically clean. Intentionally conservative: false positives here
# would blame innocent pathways. Match semantics mirror `_COFACTOR_MATCHERS`
# (equality or "<alt> (" prefix for MetaCyc charge / stereo variants).
_TOXIC_INTERMEDIATE_MATCHERS: list[tuple[str, frozenset[str]]] = [
    ("formaldehyde",           frozenset({"formaldehyde", "methanal"})),
    ("acetaldehyde",           frozenset({"acetaldehyde", "ethanal"})),
    ("glycolaldehyde",         frozenset({"glycolaldehyde"})),
    ("methylglyoxal",          frozenset({"methylglyoxal"})),
    ("catechol",               frozenset({"catechol", "pyrocatechol",
                                          "1,2-benzenediol", "1,2-dihydroxybenzene"})),
    ("hydroquinone",           frozenset({"hydroquinone", "p-hydroquinone",
                                          "1,4-benzenediol", "1,4-dihydroxybenzene"})),
    ("acrylate",               frozenset({"acrylate", "acrylic acid", "prop-2-enoate"})),
    ("4-hydroxybenzaldehyde",  frozenset({"4-hydroxybenzaldehyde"})),
    ("salicylaldehyde",        frozenset({"salicylaldehyde", "2-hydroxybenzaldehyde"})),
    ("menadione",              frozenset({"menadione", "vitamin k3",
                                          "2-methyl-1,4-naphthoquinone"})),
]


@dataclass
class _State:
    data_dir: Path = DEFAULT_DATA_DIR
    chemicals: dict[int, Chemical] | None = None
    reactions: list[Reaction] | None = None
    natives: set[Chemical] | None = None
    universals: set[Chemical] | None = None
    currency: set[Chemical] | None = None
    hypergraph: HyperGraph | None = None
    fed_hypergraphs: dict[tuple[int, ...], HyperGraph] = field(default_factory=dict)
    ec_names: dict[str, str] | None = None
    named_cofactor_ids: dict[str, frozenset[int]] | None = None
    toxic_intermediate_map: dict[int, str] | None = None
    fingerprint_index: list[tuple[int, object]] | None = None
    fingerprint_build_ms: float = 0.0
    fingerprint_skipped: int = 0
    load_time_ms: float = 0.0
    rhea_carrier_map: dict[int, frozenset[str]] | None = None


_state = _State()
_lock = threading.Lock()


def configure(data_dir: Path | str | None = None) -> None:
    """Override the data directory. Must be called before first access."""
    env = os.environ.get("SYNTHESIS_DATA_DIR")
    if data_dir is not None:
        _state.data_dir = Path(data_dir)
    elif env:
        _state.data_dir = Path(env)


def reset() -> None:
    """Clear all cached state. Intended for tests."""
    with _lock:
        _state.chemicals = None
        _state.reactions = None
        _state.natives = None
        _state.universals = None
        _state.currency = None
        _state.hypergraph = None
        _state.fed_hypergraphs.clear()
        _state.ec_names = None
        _state.named_cofactor_ids = None
        _state.toxic_intermediate_map = None
        _state.fingerprint_index = None
        _state.fingerprint_build_ms = 0.0
        _state.fingerprint_skipped = 0
        _state.load_time_ms = 0.0
        _state.rhea_carrier_map = None


def _parse_ec_names(path: Path) -> dict[str, str]:
    """Parse optional ec_names.tsv (columns: ecnum<TAB>name). Missing file → {}."""
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.lower().startswith("ecnum"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            ec, name = parts[0].strip(), parts[1].strip()
            if ec and name:
                out[ec] = name
    return out


def _bootstrap() -> None:
    """Parse files, seed shell 0, run BFS. Idempotent under lock."""
    if _state.hypergraph is not None:
        return

    configure()
    data = _state.data_dir
    t0 = time.perf_counter()

    chemicals = parse_chemicals(data / "good_chems.txt")
    reactions = parse_reactions(data / "good_reactions.txt", chemicals)
    natives = parse_metabolite_list(data / "minimal_metabolites.txt", chemicals)
    universals = parse_metabolite_list(
        data / "ubiquitous_metabolites.txt", chemicals
    )
    currency = parse_metabolite_list(
        data / "ubiquitous_metabolites.txt",
        chemicals,
        descriptor_filter=CURRENCY_DESCRIPTORS,
    )

    hg = synthesize(reactions, natives, universals, verbose=False)

    _state.chemicals = chemicals
    _state.reactions = reactions
    _state.natives = natives
    _state.universals = universals
    _state.currency = currency
    _state.hypergraph = hg
    _state.ec_names = _parse_ec_names(data / "ec_names.tsv")
    _state.load_time_ms = (time.perf_counter() - t0) * 1000.0


def get_hypergraph() -> HyperGraph:
    with _lock:
        if _state.hypergraph is None:
            _bootstrap()
        assert _state.hypergraph is not None
        return _state.hypergraph


def get_chemicals() -> dict[int, Chemical]:
    with _lock:
        if _state.chemicals is None:
            _bootstrap()
        assert _state.chemicals is not None
        return _state.chemicals


def get_reactions() -> list[Reaction]:
    with _lock:
        if _state.reactions is None:
            _bootstrap()
        assert _state.reactions is not None
        return _state.reactions


def get_natives() -> set[Chemical]:
    with _lock:
        if _state.natives is None:
            _bootstrap()
        assert _state.natives is not None
        return _state.natives


def get_universals() -> set[Chemical]:
    with _lock:
        if _state.universals is None:
            _bootstrap()
        assert _state.universals is not None
        return _state.universals


def get_currency_chemicals() -> set[Chemical]:
    """Return the subset of ubiquitous metabolites classified as currency
    (inorganics, cofactors, nucleotides) — the "hub" chemicals that tend to
    connect to most reactions and clutter visualizations.
    """
    with _lock:
        if _state.currency is None:
            _bootstrap()
        assert _state.currency is not None
        return _state.currency


def get_fingerprint_index() -> list[tuple[int, object]]:
    """Lazy-built Morgan-fingerprint index over the chemical corpus.

    First call imports RDKit and fingerprints every chemical (~3–5 s for
    9361 mols). Subsequent calls return the cached list in-process. This
    is intentionally outside of ``_bootstrap()`` so baseline MCP startup
    never pays the RDKit cost — only users of ``find_by_structure`` do.

    Each chemical is attempted as SMILES first, then InChI. Failures are
    silently skipped and counted in ``_state.fingerprint_skipped``.
    """
    with _lock:
        if _state.fingerprint_index is None:
            if _state.chemicals is None:
                _bootstrap()
            assert _state.chemicals is not None
            # Local import keeps RDKit out of the hot-path startup.
            from synthesis_helper.mcp import similarity

            index, elapsed_ms, skipped = similarity.build_index(
                _state.chemicals.values()
            )
            _state.fingerprint_index = index
            _state.fingerprint_build_ms = elapsed_ms
            _state.fingerprint_skipped = skipped
            print(
                f"[synthesis_helper] fingerprint index built: "
                f"{len(index)}/{len(_state.chemicals)} chemicals "
                f"({skipped} skipped) in {elapsed_ms:.0f} ms",
                file=sys.stderr,
            )
        return _state.fingerprint_index


def get_named_cofactor_ids() -> dict[str, frozenset[int]]:
    """Canonical cofactor label → ids of chems whose name matches.

    Keys: ``NADPH``, ``NADH``, ``NADP+``, ``NAD+``, ``ATP``, ``CoA``, ``SAM``.
    Each chemical is assigned to the FIRST matching label only (no double
    counting — NADPH is matched before NADP+ because "nadp" is a substring
    of "nadph"). Every key is always present in the returned dict, even if
    empty, so downstream callers can iterate keys safely without KeyError.

    Source pool: ``get_currency_chemicals()``. MetaCyc stores multiple
    ionisation/stereo variants of the same cofactor (NAD+ (0), NAD+ (p+1),
    SAM (27-), SAM (27+) …) so we use substring matching to pull them all
    under one canonical label.
    """
    with _lock:
        if _state.named_cofactor_ids is None:
            if _state.currency is None:
                _bootstrap()
            assert _state.currency is not None
            assigned: set[int] = set()
            buckets: dict[str, set[int]] = {key: set() for key, _ in _COFACTOR_MATCHERS}
            for chem in _state.currency:
                lname = (chem.name or "").lower()
                if not lname:
                    continue
                for label, alts in _COFACTOR_MATCHERS:
                    if _name_matches_cofactor(lname, alts):
                        if chem.id not in assigned:
                            buckets[label].add(chem.id)
                            assigned.add(chem.id)
                        break
            _state.named_cofactor_ids = {k: frozenset(v) for k, v in buckets.items()}
        return _state.named_cofactor_ids


def get_toxic_intermediate_map() -> dict[int, str]:
    """Chem id → canonical toxic-intermediate label for every chemical in the
    corpus whose name matches ``_TOXIC_INTERMEDIATE_MATCHERS``.

    Scans the full chemical corpus (not just currency) since toxic
    intermediates are, by definition, NOT cofactors — they're pathway-specific
    metabolites. Each chem is assigned to at most one label (first match
    wins, same equality-or-paren-variant semantics as cofactor matching).
    """
    with _lock:
        if _state.toxic_intermediate_map is None:
            if _state.chemicals is None:
                _bootstrap()
            assert _state.chemicals is not None
            assigned: dict[int, str] = {}
            for chem in _state.chemicals.values():
                lname = (chem.name or "").lower()
                if not lname:
                    continue
                for label, alts in _TOXIC_INTERMEDIATE_MATCHERS:
                    if _name_matches_cofactor(lname, alts):
                        assigned.setdefault(chem.id, label)
                        break
            _state.toxic_intermediate_map = assigned
        return _state.toxic_intermediate_map


def get_ec_names() -> dict[str, str]:
    """Return EC → enzyme-name lookup. Empty dict if data/ec_names.tsv absent."""
    with _lock:
        if _state.ec_names is None:
            _bootstrap()
        assert _state.ec_names is not None
        return _state.ec_names


def _parse_rhea_carrier_map(path: Path) -> dict[int, frozenset[str]]:
    """Parse rhea_carrier_map.tsv → {rxn_id: frozenset of normalized InChI strings}."""
    result: dict[int, frozenset[str]] = {}
    if not path.exists():
        return result
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("rxnid"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                rxn_id = int(parts[0])
            except ValueError:
                continue
            inchis = frozenset(i for i in parts[1].split(",") if i)
            result[rxn_id] = inchis
    return result


def get_rhea_carrier_map() -> dict[int, frozenset[str]]:
    """Lazy-loaded Rhea carrier map: reaction id → frozenset of carrier InChI strings.

    Returns an empty dict if ``data/rhea_carrier_map.tsv`` has not been built yet
    (run ``scripts/build_rhea_carrier_map.py`` to generate it). When the map is
    absent, ``rhea_atom_filter`` falls back to ``shell_zero_filter`` per-reaction,
    so the MCP server starts normally either way.
    """
    with _lock:
        if _state.rhea_carrier_map is None:
            path = _state.data_dir / "rhea_carrier_map.tsv"
            _state.rhea_carrier_map = _parse_rhea_carrier_map(path)
            if _state.rhea_carrier_map:
                print(
                    f"[synthesis_helper] Rhea carrier map loaded: "
                    f"{len(_state.rhea_carrier_map)} reactions",
                    file=sys.stderr,
                )
            else:
                print(
                    "[synthesis_helper] rhea_carrier_map.tsv not found — "
                    "rhea_atom_filter will fall back to shell_zero_filter. "
                    "Run scripts/build_rhea_carrier_map.py to generate it.",
                    file=sys.stderr,
                )
        return _state.rhea_carrier_map


def get_stats() -> dict[str, Any]:
    hg = get_hypergraph()
    shells = max(hg.reaction_to_shell.values(), default=0)
    out: dict[str, Any] = {
        "shells": shells,
        "total_chemicals": len(_state.chemicals or {}),
        "reachable_chemicals": len(hg.chemical_to_shell),
        "total_reactions": len(_state.reactions or []),
        "reachable_reactions": len(hg.reaction_to_shell),
        "load_time_ms": round(_state.load_time_ms, 1),
        "data_dir": str(_state.data_dir),
    }
    # Only surface fingerprint stats if the index has actually been built
    # (i.e. find_by_structure has been called at least once in this process).
    if _state.fingerprint_index is not None:
        indexed = len(_state.fingerprint_index)
        total = len(_state.chemicals or {})
        out["fingerprint_coverage"] = {
            "indexed": indexed,
            "total": total,
            "skipped": _state.fingerprint_skipped,
            "build_ms": round(_state.fingerprint_build_ms, 1),
            "fraction": round(indexed / total, 4) if total else 0.0,
        }
    return out


def get_or_build_fed_hypergraph(fed_ids: tuple[int, ...]) -> HyperGraph:
    """Return a cached (or freshly computed) HyperGraph with fed chemicals."""
    with _lock:
        if fed_ids in _state.fed_hypergraphs:
            return _state.fed_hypergraphs[fed_ids]

        if _state.hypergraph is None:
            _bootstrap()
        assert _state.chemicals is not None
        assert _state.reactions is not None
        assert _state.natives is not None
        assert _state.universals is not None

        extra = {_state.chemicals[cid] for cid in fed_ids if cid in _state.chemicals}
        hg = synthesize(
            _state.reactions,
            native_metabolites=_state.natives | extra,
            universal_metabolites=_state.universals,
            verbose=False,
        )
        _state.fed_hypergraphs[fed_ids] = hg
        return hg
