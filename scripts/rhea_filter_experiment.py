"""Measure cascade-size reduction from the Rhea carrier-map filter.

Loads data/rhea_carrier_map.tsv and runs traceback on a stratified sample of
target chemicals with four filter configurations, reporting how much each one
shrinks cascades relative to the shell-zero baseline.

Also reports carrier-map *coverage*: what fraction of the reactions visited
during a baseline traceback have an entry in the Rhea map (vs. falling back to
shell_zero_filter). Low coverage means the Rhea filter's benefit is limited by
match rate, not carrier-detection quality.

Filter configurations
---------------------
  baseline   -- shell_zero_filter only
  rhea       -- rhea_atom_filter (falls back to shell_zero for unmatched rxns)
  inchi_norm -- inchi_normalized_filter (strips /p /q layers against ubiquitous set)
  rhea+inchi -- compose(rhea_atom_filter, inchi_normalized_filter)

Usage
-----
    uv run python scripts/rhea_filter_experiment.py
    uv run python scripts/rhea_filter_experiment.py --per-shell 10 --max-shell 10
    uv run python scripts/rhea_filter_experiment.py --targets 317157 265539
    uv run python scripts/rhea_filter_experiment.py --verbose   # per-target coverage
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from synthesis_helper.filters import (
    compose,
    inchi_normalized_filter,
    rhea_atom_filter,
    shell_zero_filter,
)
from synthesis_helper.models import Chemical, HyperGraph, Reaction
from synthesis_helper.parser import parse_chemicals, parse_metabolite_list, parse_reactions
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback


DEFAULT_DATA_DIR = Path(__file__).resolve().parents[1] / "data"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _parse_carrier_map(path: Path) -> dict[int, frozenset[str]]:
    if not path.exists():
        print(f"  Warning: {path} not found — Rhea filter will always fall back to shell_zero",
              file=sys.stderr)
        return {}
    result: dict[int, frozenset[str]] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("rxnid"):
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            try:
                rxn_id = int(parts[0])
            except ValueError:
                continue
            inchis = frozenset(i for i in parts[1].split(",") if i)
            result[rxn_id] = inchis
    return result


def load_data(data_dir: Path) -> tuple[dict[int, Chemical], HyperGraph, set[Chemical], dict[int, frozenset[str]]]:
    print("Loading chemicals and reactions ...", file=sys.stderr)
    chemicals = parse_chemicals(data_dir / "good_chems.txt")
    reactions  = parse_reactions(data_dir / "good_reactions.txt", chemicals)
    native     = parse_metabolite_list(data_dir / "minimal_metabolites.txt", chemicals)

    univ_path = data_dir / "ubiquitous_metabolites.txt"
    universal = parse_metabolite_list(univ_path, chemicals) if univ_path.exists() else set()

    print(
        f"  {len(chemicals):,} chemicals  {len(reactions):,} reactions  "
        f"{len(native)} native  {len(universal)} universal",
        file=sys.stderr,
    )

    print("Running BFS expansion ...", file=sys.stderr)
    hg = synthesize(reactions, native, universal, verbose=False)
    print(
        f"  {len(hg.chemical_to_shell):,} reachable  "
        f"{len(hg.reaction_to_shell):,} reactions  "
        f"{max(hg.chemical_to_shell.values(), default=0)} shells",
        file=sys.stderr,
    )

    carrier_map = _parse_carrier_map(data_dir / "rhea_carrier_map.tsv")
    n_rxns_in_map = sum(1 for rxn in hg.reaction_to_shell if rxn.id in carrier_map)
    pct = 100 * n_rxns_in_map / max(len(hg.reaction_to_shell), 1)
    print(f"  Rhea map: {len(carrier_map):,} entries — "
          f"{n_rxns_in_map:,} / {len(hg.reaction_to_shell):,} "
          f"hypergraph reactions matched ({pct:.1f}%)",
          file=sys.stderr)

    return chemicals, hg, universal, carrier_map


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------

def sample_targets(
    hg: HyperGraph,
    per_shell: int,
    min_shell: int,
    max_shell: int,
) -> list[Chemical]:
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


# ---------------------------------------------------------------------------
# Measurement helpers
# ---------------------------------------------------------------------------

def run_traceback(hg: HyperGraph, target: Chemical, substrate_filter, max_producers: int) -> int:
    try:
        c = traceback(hg, target, max_producers_per_chemical=max_producers,
                      substrate_filter=substrate_filter)
        return len(c.reactions)
    except ValueError:
        return -1


def reaction_coverage(hg: HyperGraph, target: Chemical, carrier_map: dict, max_producers: int) -> tuple[int, int]:
    """(n_matched, n_total) reactions visited in baseline cascade."""
    try:
        cascade = traceback(hg, target, max_producers_per_chemical=max_producers,
                            substrate_filter=shell_zero_filter)
    except ValueError:
        return 0, 0
    total   = len(cascade.reactions)
    matched = sum(1 for r in cascade.reactions if r.id in carrier_map)
    return matched, total


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=str(DEFAULT_DATA_DIR))
    ap.add_argument("--per-shell", type=int, default=5,
                    help="Targets sampled per shell (default: 5)")
    ap.add_argument("--min-shell", type=int, default=2)
    ap.add_argument("--max-shell", type=int, default=8)
    ap.add_argument("--max-producers", type=int, default=10,
                    help="max_producers_per_chemical passed to traceback (default: 10)")
    ap.add_argument("--targets", nargs="+", type=int, metavar="ID",
                    help="Specific chemical ids (overrides sampling)")
    ap.add_argument("--verbose", action="store_true",
                    help="Print per-target Rhea coverage stats to stderr")
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    if not (data_dir / "good_chems.txt").exists():
        sys.exit(f"Data not found at {data_dir}. Pass --data-dir.")

    chemicals, hg, universal, carrier_map = load_data(data_dir)

    if args.targets:
        targets = [chemicals[i] for i in args.targets if i in chemicals]
        missing = [i for i in args.targets if i not in chemicals]
        if missing:
            print(f"Warning: unknown ids {missing}", file=sys.stderr)
    else:
        targets = sample_targets(hg, args.per_shell, args.min_shell, args.max_shell)

    print(f"\nComparing filters on {len(targets)} targets ...\n", file=sys.stderr)

    f_baseline   = shell_zero_filter
    f_rhea       = rhea_atom_filter(carrier_map)
    f_inchi      = inchi_normalized_filter(universal)
    f_rhea_inchi = compose(rhea_atom_filter(carrier_map), inchi_normalized_filter(universal))

    # TSV header
    print("\t".join([
        "id", "name", "shell",
        "cascade_baseline", "cascade_rhea", "cascade_inchi", "cascade_rhea+inchi",
        "dropped_rhea", "dropped_inchi", "dropped_rhea+inchi",
        "pct_dropped_rhea", "pct_dropped_inchi", "pct_dropped_rhea+inchi",
        "rhea_matched_rxns", "rhea_total_rxns", "rhea_coverage_pct",
    ]))

    agg_baseline = agg_rhea = agg_inchi = agg_combo = 0
    agg_matched = agg_total = 0
    n = 0

    for target in targets:
        shell = hg.chemical_to_shell.get(target, -1)
        b  = run_traceback(hg, target, f_baseline,   args.max_producers)
        r  = run_traceback(hg, target, f_rhea,        args.max_producers)
        i  = run_traceback(hg, target, f_inchi,       args.max_producers)
        ri = run_traceback(hg, target, f_rhea_inchi,  args.max_producers)

        matched, total = reaction_coverage(hg, target, carrier_map, args.max_producers)

        def _drop(x):  return (b - x) if b >= 0 and x >= 0 else ""
        def _pct(x):   return f"{100*(b-x)/b:.1f}" if b > 0 and x >= 0 else ""
        cov_pct = f"{100*matched/total:.1f}" if total > 0 else "0.0"

        if args.verbose:
            print(f"  {target.name}: Rhea coverage {matched}/{total} rxns ({cov_pct}%)",
                  file=sys.stderr)

        print("\t".join(str(v) for v in [
            target.id, target.name, shell,
            b, r, i, ri,
            _drop(r), _drop(i), _drop(ri),
            _pct(r), _pct(i), _pct(ri),
            matched, total, cov_pct,
        ]))

        if b >= 0:
            agg_baseline += b
            agg_rhea     += r if r >= 0 else b
            agg_inchi    += i if i >= 0 else b
            agg_combo    += ri if ri >= 0 else b
            agg_matched  += matched
            agg_total    += total
            n += 1

    if n:
        print(f"\n--- Aggregate over {n} targets ---", file=sys.stderr)
        print(f"  Baseline total cascade reactions : {agg_baseline}", file=sys.stderr)
        def _fmt(label, agg):
            drop = agg_baseline - agg
            pct  = 100 * drop / agg_baseline if agg_baseline else 0
            print(f"  {label:30s}: {agg:6d}  ({drop:+d}, {pct:.1f}% reduction)", file=sys.stderr)
        _fmt("rhea_atom_filter", agg_rhea)
        _fmt("inchi_normalized_filter", agg_inchi)
        _fmt("rhea + inchi (composed)", agg_combo)
        cov = 100 * agg_matched / agg_total if agg_total else 0
        print(f"  Rhea map coverage               : {agg_matched}/{agg_total} ({cov:.1f}% of baseline cascade rxns)", file=sys.stderr)


if __name__ == "__main__":
    main()
