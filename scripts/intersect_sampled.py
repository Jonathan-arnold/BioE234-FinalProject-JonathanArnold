"""Intersect shell-cutoff samples with enzymemap-shell0 reachables.

Reads data/sampled_reachables.tsv (produced by scripts/sample_reachables.py
on the shell-cutoff branch), builds the enzymemap-shell0 hypergraph (two-pass
BFS as in scripts/evaluate.py), and keeps only samples whose stripped InChI
matches a reachable (shell >= 0) in the enzymemap-shell0 hypergraph. Writes
the intersection to data/comparison_targets.tsv.

Usage:
    python scripts/intersect_sampled.py
"""

from __future__ import annotations

from pathlib import Path

from synthesis_helper.parser import (
    _strip_stereo,
    parse_chemicals,
    parse_metabolite_list,
    parse_reaction_organisms,
    parse_reactions,
)
from synthesis_helper.synthesize import synthesize

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
SAMPLES_PATH = DATA_DIR / "sampled_reachables.tsv"
OUT_PATH = DATA_DIR / "comparison_targets.tsv"
ECOLI_TOKEN = "escherichia coli"


def build_hypergraph():
    chemicals = parse_chemicals(DATA_DIR / "enzymemap_chems.tsv")
    reactions = parse_reactions(DATA_DIR / "enzymemap_reactions.tsv", chemicals)
    rxn_orgs = parse_reaction_organisms(
        DATA_DIR / "enzymemap_reaction_organisms.tsv"
    )
    native = parse_metabolite_list(DATA_DIR / "minimal_metabolites.txt", chemicals)
    universal = parse_metabolite_list(
        DATA_DIR / "ubiquitous_metabolites.txt", chemicals
    )
    print(f"Loaded {len(chemicals)} chemicals, {len(reactions)} reactions")
    ecoli_reactions = [
        r
        for r in reactions
        if any(ECOLI_TOKEN in o.lower() for o in rxn_orgs.get(r.id, ()))
    ]
    print(f"Pass 1 — E. coli reactions: {len(ecoli_reactions)}")
    hg_seed = synthesize(ecoli_reactions, native, universal, verbose=False)
    seed = set(hg_seed.chemical_to_shell.keys())
    print(f"Pass 2 — full EnzymeMap, shell 0 = {len(seed)} chemicals")
    hg = synthesize(reactions, seed, set(), verbose=False)
    return hg, chemicals


def main() -> None:
    hg, chemicals = build_hypergraph()

    # Build stripped-inchi -> chemical index over reachables only.
    inchi_to_chem: dict[str, object] = {}
    collisions = 0
    for c in hg.chemical_to_shell:
        if not c.inchi:
            continue
        key = _strip_stereo(c.inchi)
        if key in inchi_to_chem:
            collisions += 1
            continue
        inchi_to_chem[key] = c
    print(
        f"Indexed {len(inchi_to_chem)} reachable InChIs "
        f"({collisions} stripped-InChI collisions skipped)"
    )

    matched_rows = []
    n_total = 0
    n_no_inchi = 0
    n_no_match = 0
    n_shell_zero = 0

    with SAMPLES_PATH.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        i_id = header.index("id")
        i_shell = header.index("shell")
        i_name = header.index("name")
        i_stripped = header.index("inchi_stripped")
        for line in f:
            n_total += 1
            parts = line.rstrip("\n").split("\t")
            shell_cutoff_id = parts[i_id]
            shell_cutoff_shell = parts[i_shell]
            name = parts[i_name]
            stripped = parts[i_stripped]
            if not stripped:
                n_no_inchi += 1
                continue
            chem = inchi_to_chem.get(stripped)
            if chem is None:
                n_no_match += 1
                continue
            if hg.chemical_to_shell[chem] == 0:
                n_shell_zero += 1
                continue
            matched_rows.append(
                (
                    shell_cutoff_id,
                    shell_cutoff_shell,
                    chem.id,
                    hg.chemical_to_shell[chem],
                    name,
                    chem.name,
                    stripped,
                )
            )

    with OUT_PATH.open("w") as f:
        f.write(
            "shell_cutoff_id\tshell_cutoff_shell\tenzymemap_id\tenzymemap_shell\t"
            "shell_cutoff_name\tenzymemap_name\tinchi_stripped\n"
        )
        for row in matched_rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(
        f"\nSamples: {n_total} | matched: {len(matched_rows)} | "
        f"no InChI: {n_no_inchi} | no match in enzymemap: {n_no_match} | "
        f"dropped enzymemap shell 0: {n_shell_zero}"
    )
    print(f"Wrote {len(matched_rows)} comparison targets to {OUT_PATH}")


if __name__ == "__main__":
    main()
