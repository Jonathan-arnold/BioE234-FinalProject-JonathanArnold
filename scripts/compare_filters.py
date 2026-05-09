"""Compare cascade sizes across filter approaches for a sample of target chemicals.

Runs traceback on N targets sampled from each shell with five filter variants
and prints a TSV to stdout.

Usage:
    uv run python scripts/compare_filters.py
    uv run python scripts/compare_filters.py --data-dir /path/to/data --per-shell 10
    uv run python scripts/compare_filters.py --targets 317157 265539  # specific ids

Columns:
    id, name, shell,
    cascade_baseline, cascade_thresh1, cascade_thresh2,
    cascade_inchi_norm, cascade_composed,
    dropped_thresh1, dropped_thresh2, dropped_inchi_norm, dropped_composed
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running from repo root or scripts/ dir
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from synthesis_helper.filters import (
    compose,
    inchi_normalized_filter,
    shell_threshold_filter,
    shell_zero_filter,
)
from synthesis_helper.models import Chemical, HyperGraph
from synthesis_helper.parser import parse_chemicals, parse_metabolite_list, parse_reactions
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback


DEFAULT_DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def load_data(data_dir: Path) -> tuple[dict[int, Chemical], HyperGraph, set[Chemical]]:
    print("Loading chemicals ...", file=sys.stderr)
    chemicals = parse_chemicals(data_dir / "good_chems.txt")
    reactions = parse_reactions(data_dir / "good_reactions.txt", chemicals)
    native = parse_metabolite_list(data_dir / "minimal_metabolites.txt", chemicals)

    univ_path = data_dir / "ecoli_reachables_shell0.txt"
    if not univ_path.exists():
        univ_path = data_dir / "ubiquitous_metabolites.txt"
    universal = parse_metabolite_list(univ_path, chemicals)

    print(
        f"  {len(chemicals):,} chemicals  {len(reactions):,} reactions  "
        f"{len(native)} native  {len(universal)} universal",
        file=sys.stderr,
    )
    print("Running BFS expansion ...", file=sys.stderr)
    hg = synthesize(reactions, native, universal, verbose=False)
    print(
        f"  {len(hg.chemical_to_shell):,} reachable  {len(hg.reaction_to_shell):,} reactions  "
        f"{max(hg.chemical_to_shell.values(), default=0)} shells",
        file=sys.stderr,
    )
    return chemicals, hg, universal


def sample_targets(
    hg: HyperGraph,
    per_shell: int,
    min_shell: int,
    max_shell: int,
) -> list[Chemical]:
    """Return up to per_shell chemicals evenly spread across each shell level."""
    by_shell: dict[int, list[Chemical]] = {}
    for chem, shell in hg.chemical_to_shell.items():
        if min_shell <= shell <= max_shell:
            by_shell.setdefault(shell, []).append(chem)

    targets: list[Chemical] = []
    for shell in sorted(by_shell):
        bucket = sorted(by_shell[shell], key=lambda c: c.id)
        step = max(1, len(bucket) // per_shell)
        targets.extend(bucket[::step][:per_shell])
    return targets


def cascade_size(
    hg: HyperGraph,
    target: Chemical,
    substrate_filter,
    max_producers: int,
) -> int:
    try:
        c = traceback(hg, target, max_producers_per_chemical=max_producers,
                      substrate_filter=substrate_filter)
        return len(c.reactions)
    except ValueError:
        return -1


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=str(DEFAULT_DATA_DIR),
                    help="Path to directory containing good_chems.txt etc.")
    ap.add_argument("--per-shell", type=int, default=5,
                    help="Targets sampled per shell level (default: 5)")
    ap.add_argument("--min-shell", type=int, default=2,
                    help="Minimum shell to include (default: 2)")
    ap.add_argument("--max-shell", type=int, default=8,
                    help="Maximum shell to include (default: 8)")
    ap.add_argument("--max-producers", type=int, default=10,
                    help="max_producers_per_chemical passed to traceback (default: 10)")
    ap.add_argument("--targets", nargs="+", type=int, metavar="ID",
                    help="Specific chemical ids to compare (overrides sampling)")
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    if not (data_dir / "good_chems.txt").exists():
        sys.exit(f"Data not found at {data_dir}. Pass --data-dir or run download_data.py.")

    chemicals, hg, universal = load_data(data_dir)

    if args.targets:
        targets = [chemicals[i] for i in args.targets if i in chemicals]
        missing = [i for i in args.targets if i not in chemicals]
        if missing:
            print(f"Warning: unknown ids {missing}", file=sys.stderr)
    else:
        targets = sample_targets(hg, args.per_shell, args.min_shell, args.max_shell)

    print(f"Comparing filters on {len(targets)} targets ...", file=sys.stderr)

    filters = {
        "baseline":   shell_zero_filter,
        "thresh1":    shell_threshold_filter(1),
        "thresh2":    shell_threshold_filter(2),
        "inchi_norm": inchi_normalized_filter(universal),
        "composed":   compose(inchi_normalized_filter(universal), shell_threshold_filter(1)),
    }
    filter_keys = list(filters)

    header = (
        "id\tname\tshell\t"
        + "\t".join(f"cascade_{k}" for k in filter_keys)
        + "\t"
        + "\t".join(f"dropped_{k}" for k in filter_keys[1:])
    )
    print(header)

    for target in targets:
        shell = hg.chemical_to_shell.get(target, -1)
        sizes = {k: cascade_size(hg, target, f, args.max_producers)
                 for k, f in filters.items()}
        baseline = sizes["baseline"]
        dropped = {k: (baseline - sizes[k] if baseline >= 0 and sizes[k] >= 0 else "")
                   for k in filter_keys[1:]}

        row = (
            f"{target.id}\t{target.name}\t{shell}\t"
            + "\t".join(str(sizes[k]) for k in filter_keys)
            + "\t"
            + "\t".join(str(dropped[k]) for k in filter_keys[1:])
        )
        print(row)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
