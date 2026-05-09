"""Audit InChI matching quality between ubiquitous_metabolites.txt and good_chems.txt.

For each entry in the universal metabolite list, tries four progressively looser
InChI normalization strategies and reports which chemicals are matched or missed
at each level. This quantifies how much of the "cofactor in hypergraph" problem
is caused by protonation/charge/stereo mismatches vs. genuine gaps.

Usage:
    uv run python scripts/inchi_mismatch_audit.py
    uv run python scripts/inchi_mismatch_audit.py --data-dir /path/to/data --verbose
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from synthesis_helper.parser import _read_lines, parse_chemicals


DEFAULT_DATA_DIR = Path(__file__).resolve().parents[1] / "data"


# ---------------------------------------------------------------------------
# InChI normalization strategies (progressively looser)
# ---------------------------------------------------------------------------

def exact(inchi: str) -> str:
    return inchi.strip('"')

def strip_proton(inchi: str) -> str:
    return re.sub(r"/p[+-]\d+", "", exact(inchi))

def strip_proton_charge(inchi: str) -> str:
    s = re.sub(r"/p[+-]\d+", "", exact(inchi))
    return re.sub(r"/q[+-]\d+", "", s)

def strip_all_variable(inchi: str) -> str:
    s = strip_proton_charge(inchi)
    s = re.sub(r"/t[^/]+", "", s)
    s = re.sub(r"/m\d+", "", s)
    s = re.sub(r"/s\d+", "", s)
    return s

STRATEGIES: list[tuple[str, callable]] = [
    ("exact",              exact),
    ("strip_p",            strip_proton),
    ("strip_p_q",          strip_proton_charge),
    ("strip_p_q_stereo",   strip_all_variable),
]


# ---------------------------------------------------------------------------
# Parsing the metabolite list (minimal version — no chemical lookup needed)
# ---------------------------------------------------------------------------

@dataclass
class MetaboliteEntry:
    name: str
    inchi: str
    descriptor: str


def parse_metabolite_entries(filepath: Path) -> list[MetaboliteEntry]:
    entries = []
    for line in _read_lines(filepath):
        if not line or line.startswith("name"):
            continue
        parts = line.split("\t")
        name = parts[0] if len(parts) > 0 else ""
        inchi = parts[1].strip('"') if len(parts) > 1 else ""
        descriptor = parts[2].strip() if len(parts) > 2 else ""
        entries.append(MetaboliteEntry(name=name, inchi=inchi, descriptor=descriptor))
    return entries


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=str(DEFAULT_DATA_DIR))
    ap.add_argument("--univ-file", default="ubiquitous_metabolites.txt",
                    help="Universal metabolite file relative to data-dir")
    ap.add_argument("--verbose", action="store_true",
                    help="Print per-entry match details for each strategy")
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    if not (data_dir / "good_chems.txt").exists():
        sys.exit(f"Data not found at {data_dir}. Pass --data-dir or run download_data.py.")

    print("Loading chemicals ...", file=sys.stderr)
    chemicals = parse_chemicals(data_dir / "good_chems.txt")
    print(f"  {len(chemicals):,} chemicals loaded", file=sys.stderr)

    univ_path = data_dir / args.univ_file
    entries = parse_metabolite_entries(univ_path)
    print(f"  {len(entries)} universal metabolite entries from {univ_path.name}", file=sys.stderr)

    # Build one lookup per strategy: normalized_inchi -> [Chemical, ...]
    lookups: dict[str, dict[str, list]] = {name: {} for name, _ in STRATEGIES}
    for chem in chemicals.values():
        raw = chem.inchi.strip('"')
        for strat_name, norm_fn in STRATEGIES:
            key = norm_fn(raw)
            if key:
                lookups[strat_name].setdefault(key, []).append(chem)

    # Also build name lookup as a last-resort baseline
    name_lookup = {c.name.lower(): c for c in chemicals.values() if c.name}

    # ---------------------------------------------------------------------------
    # Per-entry matching: find first strategy that hits, track miss cascade
    # ---------------------------------------------------------------------------

    @dataclass
    class Result:
        entry: MetaboliteEntry
        first_hit_strategy: str | None     # None = unmatched by any strategy
        first_hit_count: int               # number of chemicals matched at that level
        name_match: bool                   # whether name lookup also finds it
        hit_by: dict[str, int]             # strategy -> count (0 = miss)

    results: list[Result] = []
    for entry in entries:
        hit_by: dict[str, int] = {}
        first_hit: str | None = None
        first_count = 0

        for strat_name, norm_fn in STRATEGIES:
            key = norm_fn(entry.inchi)
            matches = lookups[strat_name].get(key, [])
            hit_by[strat_name] = len(matches)
            if matches and first_hit is None:
                first_hit = strat_name
                first_count = len(matches)

        name_match = entry.name.lower() in name_lookup
        results.append(Result(entry, first_hit, first_count, name_match, hit_by))

    # ---------------------------------------------------------------------------
    # Summary table
    # ---------------------------------------------------------------------------

    total = len(results)
    print(f"\n{'='*60}")
    print(f"Universal metabolite match summary  (n={total})")
    print(f"{'='*60}")
    print(f"{'Strategy':<22}  {'Matched':>8}  {'Miss':>8}  {'Cumulative %':>13}")
    print(f"{'-'*22}  {'-'*8}  {'-'*8}  {'-'*13}")

    cumulative_matched: set[str] = set()
    for strat_name, _ in STRATEGIES:
        newly_matched = {
            r.entry.name for r in results
            if r.first_hit_strategy == strat_name
        }
        cumulative_matched |= newly_matched
        missed = total - len(cumulative_matched)
        pct = 100.0 * len(cumulative_matched) / total if total else 0
        print(f"  {strat_name:<20}  {len(newly_matched):>8}  {missed:>8}  {pct:>12.1f}%")

    name_only = [r for r in results if r.first_hit_strategy is None and r.name_match]
    truly_missing = [r for r in results if r.first_hit_strategy is None and not r.name_match]

    print(f"\n  Name-only fallback matches:  {len(name_only)}")
    print(f"  Completely unmatched:        {len(truly_missing)}")

    # ---------------------------------------------------------------------------
    # Incremental-gain detail: what does each looser strategy add?
    # ---------------------------------------------------------------------------

    print(f"\n{'='*60}")
    print("Incremental gain per relaxation step")
    print(f"{'='*60}")
    for strat_name, _ in STRATEGIES[1:]:
        gained = [
            r for r in results
            if r.first_hit_strategy == strat_name
        ]
        if not gained:
            print(f"  {strat_name}: no additional matches")
            continue
        print(f"\n  {strat_name}  (+{len(gained)} new matches):")
        for r in sorted(gained, key=lambda x: x.entry.name):
            n_chems = r.hit_by[strat_name]
            flag = " [multi]" if n_chems > 1 else ""
            print(f"    {r.entry.name:<40}  → {n_chems} chemical(s){flag}")

    # ---------------------------------------------------------------------------
    # Unmatched entries
    # ---------------------------------------------------------------------------

    if name_only:
        print(f"\n{'='*60}")
        print(f"Name-only matches (InChI missed by all strategies, n={len(name_only)})")
        print(f"{'='*60}")
        for r in sorted(name_only, key=lambda x: x.entry.name):
            print(f"  {r.entry.name:<40}  inchi: {r.entry.inchi[:60]}")

    if truly_missing:
        print(f"\n{'='*60}")
        print(f"Completely unmatched (not found by InChI or name, n={len(truly_missing)})")
        print(f"{'='*60}")
        for r in sorted(truly_missing, key=lambda x: x.entry.name):
            print(f"  {r.entry.name:<40}  inchi: {r.entry.inchi[:60]}")

    # ---------------------------------------------------------------------------
    # Verbose per-entry detail
    # ---------------------------------------------------------------------------

    if args.verbose:
        print(f"\n{'='*60}")
        print("Per-entry detail")
        print(f"{'='*60}")
        header = f"  {'name':<40}  " + "  ".join(f"{s:>14}" for s, _ in STRATEGIES)
        print(header)
        for r in sorted(results, key=lambda x: x.entry.name):
            counts = "  ".join(f"{r.hit_by[s]:>14}" for s, _ in STRATEGIES)
            print(f"  {r.entry.name:<40}  {counts}")


if __name__ == "__main__":
    main()
