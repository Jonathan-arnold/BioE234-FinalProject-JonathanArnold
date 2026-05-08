"""Entry point: load EnzymeMap data, run two-pass BFS, query results.

Per CLAUDE.md (approach #2 — enzymemap-shell0):

  1. Pass 1: starting from the curated minimal + ubiquitous metabolites,
     run BFS using only EnzymeMap reactions annotated as natively occurring
     in *E. coli*. Every reachable from this pass becomes the new shell 0.

  2. Pass 2: re-run BFS over the full EnzymeMap reaction set, treating the
     pass-1 reachables as shell 0. Cascade work downstream operates on the
     resulting hypergraph.
"""

from __future__ import annotations

from pathlib import Path

from synthesis_helper.parser import (
    parse_chemicals,
    parse_metabolite_list,
    parse_reaction_organisms,
    parse_reactions,
)
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback
from synthesis_helper.pathways import enumerate_pathways
from synthesis_helper.io import dump_reachables

DATA_DIR = Path(__file__).parent / "data"

ECOLI_TOKEN = "Escherichia coli"  # case-insensitive substring match


def main() -> None:
    chemicals = parse_chemicals(DATA_DIR / "enzymemap_chems.tsv")
    reactions = parse_reactions(DATA_DIR / "enzymemap_reactions.tsv", chemicals)
    rxn_orgs = parse_reaction_organisms(
        DATA_DIR / "enzymemap_reaction_organisms.tsv"
    )
    native = parse_metabolite_list(DATA_DIR / "minimal_metabolites.txt", chemicals)
    universal = parse_metabolite_list(
        DATA_DIR / "ubiquitous_metabolites.txt", chemicals
    )

    print(
        f"Loaded {len(chemicals)} chemicals, {len(reactions)} reactions "
        f"({len(rxn_orgs)} with organism annotations)"
    )
    print(f"Curated seed: {len(native)} minimal + {len(universal)} ubiquitous")

    token = ECOLI_TOKEN.lower()
    ecoli_reactions = [
        r for r in reactions
        if any(token in o.lower() for o in rxn_orgs.get(r.id, ()))
    ]
    print(f"\nPass 1 — E. coli-only reactions: {len(ecoli_reactions)}")
    hg_ecoli = synthesize(ecoli_reactions, native, universal, verbose=True)

    seed = set(hg_ecoli.chemical_to_shell.keys())
    out_path = DATA_DIR / "enzymemap_ecoli_reachables.txt"
    dump_reachables(hg_ecoli, out_path)
    print(f"  E. coli reachables → {out_path}")

    print(f"\nPass 2 — full EnzymeMap reactions, shell 0 = {len(seed)} chemicals")
    hg = synthesize(reactions, seed, set(), verbose=True)

    target_id = 317157
    if target_id in chemicals:
        target = chemicals[target_id]
        print(
            f"\nTraceback for {target.name} (shell {hg.chemical_to_shell.get(target)}):"
        )
        cascade = traceback(hg, target)
        print(f"  Cascade contains {len(cascade.reactions)} reactions")
        pathways = enumerate_pathways(cascade, hg)
        print(f"  Found {len(pathways)} pathways")


if __name__ == "__main__":
    main()
