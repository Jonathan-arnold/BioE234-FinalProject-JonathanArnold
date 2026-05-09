"""Entry point: load data, run synthesis, query results."""

from __future__ import annotations

from pathlib import Path

from synthesis_helper.parser import (
    parse_chemicals,
    parse_metabolite_list,
    parse_reactions,
)
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback
from synthesis_helper.pathways import enumerate_pathways
from synthesis_helper.composition import pathway_to_composition
from synthesis_helper.io import dump_reachables

DATA_DIR = Path(__file__).parent / "data"


def main() -> None:
    # 1. Parse input files
    chemicals = parse_chemicals(DATA_DIR / "good_chems.txt")
    reactions = parse_reactions(DATA_DIR / "good_reactions.txt", chemicals)
    native = parse_metabolite_list(DATA_DIR / "minimal_metabolites.txt", chemicals)
    universal = parse_metabolite_list(
        DATA_DIR / "ubiquitous_metabolites.txt", chemicals
    )

    print(f"Loaded {len(chemicals)} chemicals, {len(reactions)} reactions")
    print(f"Native metabolites: {len(native)}, Universal: {len(universal)}")

    # 2. Run BFS waveform expansion
    hg = synthesize(reactions, native, universal, verbose=True)

    # 3. Dump reachables
    output_path = DATA_DIR / "metacyc_L2_reachables.txt"
    dump_reachables(hg, output_path)
    print(f"\nReachables written to {output_path}")

    # 4. Example: query a target chemical
    #    Uncomment and set target_id to a chemical of interest:
    #
    target_id = 317157
    if target_id in chemicals:
        target = chemicals[target_id]
        print(
            f"\nTraceback for {target.name} (shell {hg.chemical_to_shell.get(target)}):"
        )
        cascade = traceback(hg, target, max_producers_per_chemical=25)
        print(f"  Cascade contains {len(cascade.reactions)} reactions")
        pathways = enumerate_pathways(cascade, hg)
        print(f"  Found {len(pathways)} pathways")


if __name__ == "__main__":
    main()
