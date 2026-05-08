"""Evaluate the shell-cutoff + max-producers combo across a 2-D sweep.

Builds the hypergraph once, samples reachables per shell (seeded) — or loads
matched targets via --use-matches — then for each (shell_cutoff,
max_producers_per_chemical) combination runs traceback on every target and
records cascade-shape metrics per shell.

Usage:
    python scripts/evaluate.py --max-cutoff 2 [--use-matches]
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
import time
from collections import Counter, defaultdict
from pathlib import Path

from synthesis_helper.models import Cascade, Chemical, HyperGraph, Reaction
from synthesis_helper.parser import (
    parse_chemicals,
    parse_metabolite_list,
    parse_reactions,
)
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "eval_results"
MATCHES_PATH = REPO_ROOT / "data" / "comparison_targets.tsv"
MATCHES_ID_COLUMN = "shell_cutoff_id"
SEED = 20260509
TOTAL_SAMPLES = 5
CAP_SWEEP: tuple[int, ...] = (5, 10, 15)


def build_hypergraph() -> tuple[HyperGraph, dict[int, Chemical]]:
    chemicals = parse_chemicals(DATA_DIR / "good_chems.txt")
    reactions = parse_reactions(DATA_DIR / "good_reactions.txt", chemicals)
    native = parse_metabolite_list(DATA_DIR / "minimal_metabolites.txt", chemicals)
    universal = parse_metabolite_list(
        DATA_DIR / "ubiquitous_metabolites.txt", chemicals
    )
    print(f"Loaded {len(chemicals)} chemicals, {len(reactions)} reactions")
    print(f"Native: {len(native)}, Universal: {len(universal)}")
    hg = synthesize(reactions, native, universal, verbose=True)
    return hg, chemicals


def sample_targets_per_shell(hg: HyperGraph) -> dict[int, list[Chemical]]:
    """Sample TOTAL_SAMPLES reachables uniformly from all non-shell-0 chemicals.

    Returned grouped by shell so downstream per-shell reporting still works.
    """
    pool = sorted(
        (c for c, s in hg.chemical_to_shell.items() if s > 0),
        key=lambda c: c.id,
    )
    rng = random.Random(SEED)
    k = min(TOTAL_SAMPLES, len(pool))
    chosen = rng.sample(pool, k)

    sampled: dict[int, list[Chemical]] = defaultdict(list)
    for c in chosen:
        sampled[hg.chemical_to_shell[c]].append(c)
    return dict(sorted(sampled.items()))


def load_match_targets(
    hg: HyperGraph, chemicals: dict[int, Chemical]
) -> dict[int, list[Chemical]]:
    """Load targets from data/comparison_targets.tsv, grouped by current shell.

    Uses MATCHES_ID_COLUMN to pick the correct branch-local chemical id.
    Skips rows whose id isn't a known chemical or isn't reachable in hg.
    """
    sampled: dict[int, list[Chemical]] = defaultdict(list)
    skipped_unknown = 0
    skipped_unreachable = 0
    with MATCHES_PATH.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = header.index(MATCHES_ID_COLUMN)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            try:
                chem_id = int(parts[idx])
            except (ValueError, IndexError):
                continue
            chem = chemicals.get(chem_id)
            if chem is None:
                skipped_unknown += 1
                continue
            shell = hg.chemical_to_shell.get(chem)
            if shell is None:
                skipped_unreachable += 1
                continue
            sampled[shell].append(chem)
    for chems in sampled.values():
        chems.sort(key=lambda c: c.id)
    print(
        f"Loaded {sum(len(v) for v in sampled.values())} match targets from "
        f"{MATCHES_PATH.name} (skipped {skipped_unknown} unknown id, "
        f"{skipped_unreachable} unreachable)"
    )
    return dict(sorted(sampled.items()))


def build_global_producers(hg: HyperGraph) -> dict[int, list[Reaction]]:
    """For every chemical, list every enabled reaction that has it as a product."""
    producers: dict[int, list[Reaction]] = defaultdict(list)
    for r in hg.reaction_to_shell:
        for p in r.products:
            producers[p.id].append(r)
    return producers


def target_diagnostics(
    hg: HyperGraph,
    cascade: Cascade,
    global_producers: dict[int, list[Reaction]],
    top_n: int = 10,
) -> dict:
    """Rich per-target diagnostics for one cascade.

    "Necessary" chemicals = the target plus every substrate of any cascade
    reaction (pure byproducts excluded). For each non-shell-0 necessary
    chemical we report:
      * producer_rxns_in_cascade — # cascade reactions that list it as a
        product. (Inflated by cofactor byproducts; not bounded by the cap.)
      * producer_rxns_in_hypergraph — # candidate producers globally.
    """
    necessary: set[Chemical] = {cascade.target}
    for rxn in cascade.reactions:
        necessary.update(rxn.substrates)

    shell_histogram = Counter(hg.chemical_to_shell[c] for c in necessary)

    in_cascade_producers: dict[int, int] = defaultdict(int)
    for rxn in cascade.reactions:
        for p in rxn.products:
            if p in necessary:
                in_cascade_producers[p.id] += 1

    deep_chemicals = [c for c in necessary if hg.chemical_to_shell[c] > 0]
    rows = [
        {
            "chemical_id": c.id,
            "chemical_name": c.name,
            "shell": hg.chemical_to_shell[c],
            "producer_rxns_in_cascade": in_cascade_producers.get(c.id, 0),
            "producer_rxns_in_hypergraph": len(global_producers.get(c.id, [])),
        }
        for c in deep_chemicals
    ]
    rows.sort(
        key=lambda r: (
            -r["producer_rxns_in_cascade"],
            -r["producer_rxns_in_hypergraph"],
            r["chemical_id"],
        )
    )

    branching = Counter(
        sum(1 for s in rxn.substrates if hg.chemical_to_shell[s] > 0)
        for rxn in cascade.reactions
    )

    return {
        "n_reactions": len(cascade.reactions),
        "n_necessary_chemicals": len(necessary),
        "n_necessary_non_shell_0": len(deep_chemicals),
        "necessary_chemicals_by_shell": dict(sorted(shell_histogram.items())),
        "reactions_by_n_non_shell_0_substrates": dict(sorted(branching.items())),
        "reactions_per_non_shell_0_chemical": (
            len(cascade.reactions) / len(deep_chemicals) if deep_chemicals else 0.0
        ),
        "top_chemicals": rows[:top_n],
    }


def evaluate_target(
    hg: HyperGraph,
    target: Chemical,
    global_producers: dict[int, list[Reaction]],
    shell_cutoff: int,
    max_producers: int | None,
) -> dict:
    """Run traceback on a single target and record cascade diagnostics."""
    t0 = time.perf_counter()
    cascade = traceback(
        hg,
        target,
        shell_cutoff=shell_cutoff,
        max_producers_per_chemical=max_producers,
    )
    t_cascade = time.perf_counter() - t0

    diag = target_diagnostics(hg, cascade, global_producers)
    return {
        "chemical_id": target.id,
        "chemical_name": target.name,
        "target_shell": hg.chemical_to_shell[target],
        "cascade_seconds": t_cascade,
        **diag,
    }


def print_target_report(target: Chemical, target_shell: int, r: dict) -> None:
    """Mirror scripts/debug_blowup.py's per-target printout."""
    print(
        f"\n  TARGET: {target.name} (id={target.id}, shell={target_shell})"
    )
    print(f"    cascade: {r['n_reactions']} reactions")
    print(
        f"    necessary chemicals: {r['n_necessary_chemicals']}, "
        f"by shell: {r['necessary_chemicals_by_shell']}"
    )
    print(
        f"    non-shell-0 necessary chemicals: "
        f"{r['n_necessary_non_shell_0']} / {r['n_necessary_chemicals']}"
    )
    print(
        f"    reactions per non-shell-0 chemical: "
        f"{r['reactions_per_non_shell_0_chemical']:.2f}"
    )
    print(
        f"    reactions by # non-shell-0 substrates: "
        f"{r['reactions_by_n_non_shell_0_substrates']}"
    )
    if r["top_chemicals"]:
        print(
            f"    top {len(r['top_chemicals'])} non-shell-0 necessary "
            f"chemicals by producer reactions in cascade:"
        )
        print(
            f"      {'in_cascade':>10} {'available':>10}  shell  name (id)"
        )
        for row in r["top_chemicals"]:
            print(
                f"      {row['producer_rxns_in_cascade']:>10} "
                f"{row['producer_rxns_in_hypergraph']:>10}  "
                f"{row['shell']:>5}  {row['chemical_name'][:50]} "
                f"(id={row['chemical_id']})"
            )


def aggregate(target_results: list[dict]) -> dict:
    """Summary stats across a list of per-target results."""
    if not target_results:
        return {"count": 0}
    n_rxns = [r["n_reactions"] for r in target_results]
    n_chems = [r["n_necessary_chemicals"] for r in target_results]
    n_deep = [r["n_necessary_non_shell_0"] for r in target_results]
    return {
        "count": len(target_results),
        "mean_n_reactions": statistics.mean(n_rxns),
        "median_n_reactions": statistics.median(n_rxns),
        "max_n_reactions": max(n_rxns),
        "mean_n_necessary_chemicals": statistics.mean(n_chems),
        "median_n_necessary_chemicals": statistics.median(n_chems),
        "mean_n_necessary_non_shell_0": statistics.mean(n_deep),
        "max_n_necessary_non_shell_0": max(n_deep),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--max-cutoff",
        type=int,
        default=None,
        help=(
            "If given, sweep shell_cutoff from -1 up to this value (inclusive). "
            "If omitted, run a single point with no shell cutoff."
        ),
    )
    ap.add_argument(
        "--use-matches",
        action="store_true",
        help=(
            "Use targets from data/comparison_targets.tsv (cross-branch matched "
            "reachables) instead of fresh random sampling."
        ),
    )
    args = ap.parse_args()

    if args.max_cutoff is not None and args.max_cutoff < -1:
        ap.error("--max-cutoff must be >= -1")

    RESULTS_DIR.mkdir(exist_ok=True)
    caps = list(CAP_SWEEP)

    hg, chemicals = build_hypergraph()
    global_producers = build_global_producers(hg)

    if args.max_cutoff is None:
        max_shell = max(hg.chemical_to_shell.values(), default=0)
        cutoffs: list[int] = [max_shell]
        cutoff_labels = {max_shell: "none"}
        print("No --max-cutoff given; running with no shell cutoff only.")
    else:
        cutoffs = list(range(-1, args.max_cutoff + 1))
        cutoff_labels = {n: str(n) for n in cutoffs}
        print(f"Sweeping shell_cutoff over: {cutoffs}")
    print(f"Sweeping max_producers_per_chemical over: {caps}")

    if args.use_matches:
        sampled = load_match_targets(hg, chemicals)
    else:
        sampled = sample_targets_per_shell(hg)
    print("\nSampled targets per shell:")
    for shell, chems in sorted(sampled.items()):
        print(f"  shell {shell}: {len(chems)} targets")

    results: dict = {
        "approach": "shell-cutoff",
        "target_source": "matches" if args.use_matches else "random_sample",
        "seed": SEED,
        "total_samples": (
            sum(len(v) for v in sampled.values())
            if args.use_matches
            else TOTAL_SAMPLES
        ),
        "max_shell": max(sampled) if sampled else 0,
        "cap_sweep": caps,
        "sampled_targets": {
            str(shell): [{"id": c.id, "name": c.name} for c in chems]
            for shell, chems in sampled.items()
        },
        "sweep": {},
    }

    for n in cutoffs:
        label = cutoff_labels[n]
        results["sweep"][label] = {}
        for cap in caps:
            print(f"\n=== shell_cutoff = {label}, max_producers = {cap} ===")
            per_shell: dict[int, list[dict]] = {}
            for shell, chems in sorted(sampled.items()):
                target_results = []
                for chem in chems:
                    r = evaluate_target(
                        hg,
                        chem,
                        global_producers,
                        shell_cutoff=n,
                        max_producers=cap,
                    )
                    target_results.append(r)
                    print_target_report(chem, shell, r)
                per_shell[shell] = target_results

            results["sweep"][label][str(cap)] = {
                "per_shell": {
                    str(shell): {
                        "targets": tr,
                        "summary": aggregate(tr),
                    }
                    for shell, tr in per_shell.items()
                },
                "overall": aggregate([r for tr in per_shell.values() for r in tr]),
            }

    suffix = "_matches" if args.use_matches else ""
    json_path = RESULTS_DIR / f"shell-cutoff{suffix}.json"
    json_path.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {json_path}")


if __name__ == "__main__":
    main()
