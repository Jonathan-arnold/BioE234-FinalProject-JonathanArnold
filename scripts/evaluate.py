"""Evaluate the enzymemap-shell0 approach (single run, unlimited cutoff).

Builds the hypergraph via the two-pass BFS from main.py (E. coli reactions
seed shell 0; full EnzymeMap expansion runs from there), samples reachables
uniformly across all non-shell-0 chemicals (or loads matched targets via
--use-matches), then runs traceback on every target with shell_cutoff left
unlimited and records cascade-shape metrics.

Usage:
    python scripts/evaluate.py [--use-matches]
"""

from __future__ import annotations

import argparse
import json
import statistics
import time
import random
from collections import defaultdict
from pathlib import Path

from synthesis_helper.models import Cascade, Chemical, HyperGraph
from synthesis_helper.parser import (
    parse_chemicals,
    parse_metabolite_list,
    parse_reaction_organisms,
    parse_reactions,
)
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "eval_results"
MATCHES_PATH = REPO_ROOT / "data" / "comparison_targets.tsv"
MATCHES_ID_COLUMN = "enzymemap_id"
SEED = 20260509
TOTAL_SAMPLES = 5
ECOLI_TOKEN = "escherichia coli"


def build_hypergraph() -> tuple[HyperGraph, dict[int, Chemical]]:
    chemicals = parse_chemicals(DATA_DIR / "enzymemap_chems.tsv")
    reactions = parse_reactions(DATA_DIR / "enzymemap_reactions.tsv", chemicals)
    rxn_orgs = parse_reaction_organisms(DATA_DIR / "enzymemap_reaction_organisms.tsv")
    native = parse_metabolite_list(DATA_DIR / "minimal_metabolites.txt", chemicals)
    universal = parse_metabolite_list(
        DATA_DIR / "ubiquitous_metabolites.txt", chemicals
    )
    print(f"Loaded {len(chemicals)} chemicals, {len(reactions)} reactions")
    print(f"Native: {len(native)}, Universal: {len(universal)}")

    ecoli_reactions = [
        r
        for r in reactions
        if any(ECOLI_TOKEN in o.lower() for o in rxn_orgs.get(r.id, ()))
    ]
    print(f"Pass 1 — E. coli reactions: {len(ecoli_reactions)}")
    hg_seed = synthesize(ecoli_reactions, native, universal, verbose=True)
    seed = set(hg_seed.chemical_to_shell.keys())
    print(f"Pass 2 — full EnzymeMap, shell 0 = {len(seed)} chemicals")
    hg = synthesize(reactions, seed, set(), verbose=True)
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


def cascade_shape(cascade: Cascade) -> dict:
    """Count chemicals and find the max in/out reaction counts in the cascade.

    "Necessary" chemicals are the target plus every substrate of any
    cascade reaction — pure byproducts (products not consumed downstream
    and not the target) are excluded. In/out counts are computed only over
    necessary chemicals.
    """
    necessary: set[Chemical] = {cascade.target}
    for rxn in cascade.reactions:
        necessary.update(rxn.substrates)

    in_count: dict[Chemical, int] = {c: 0 for c in necessary}
    out_count: dict[Chemical, int] = {c: 0 for c in necessary}
    for rxn in cascade.reactions:
        for p in rxn.products:
            if p in necessary:
                in_count[p] += 1
        for s in rxn.substrates:
            out_count[s] += 1

    max_in_chem = max(in_count, key=lambda c: (in_count[c], -c.id))
    max_out_chem = max(out_count, key=lambda c: (out_count[c], -c.id))
    return {
        "n_chemicals": len(necessary),
        "max_in_count": in_count[max_in_chem],
        "max_in_chemical_id": max_in_chem.id,
        "max_in_chemical_name": max_in_chem.name,
        "max_out_count": out_count[max_out_chem],
        "max_out_chemical_id": max_out_chem.id,
        "max_out_chemical_name": max_out_chem.name,
    }


def evaluate_target(hg: HyperGraph, target: Chemical) -> dict:
    """Run traceback (unlimited cutoff) and record cascade-shape metrics."""
    t0 = time.perf_counter()
    cascade = traceback(hg, target)
    t_cascade = time.perf_counter() - t0

    shape = cascade_shape(cascade)
    return {
        "chemical_id": target.id,
        "chemical_name": target.name,
        "target_shell": hg.chemical_to_shell[target],
        "n_reactions": len(cascade.reactions),
        "cascade_seconds": t_cascade,
        **shape,
    }


def aggregate(target_results: list[dict]) -> dict:
    """Summary stats across a list of per-target results."""
    if not target_results:
        return {"count": 0}
    n_rxns = [r["n_reactions"] for r in target_results]
    n_chems = [r["n_chemicals"] for r in target_results]
    max_ins = [r["max_in_count"] for r in target_results]
    max_outs = [r["max_out_count"] for r in target_results]
    return {
        "count": len(target_results),
        "mean_n_reactions": statistics.mean(n_rxns),
        "median_n_reactions": statistics.median(n_rxns),
        "mean_n_chemicals": statistics.mean(n_chems),
        "median_n_chemicals": statistics.median(n_chems),
        "max_in_count_overall": max(max_ins),
        "max_out_count_overall": max(max_outs),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--use-matches",
        action="store_true",
        help=(
            "Use targets from data/comparison_targets.tsv (cross-branch matched "
            "reachables) instead of fresh random sampling."
        ),
    )
    args = ap.parse_args()

    RESULTS_DIR.mkdir(exist_ok=True)

    hg, chemicals = build_hypergraph()

    if args.use_matches:
        sampled = load_match_targets(hg, chemicals)
    else:
        sampled = sample_targets_per_shell(hg)
    print("\nSampled targets per shell:")
    for shell, chems in sorted(sampled.items()):
        print(f"  shell {shell}: {len(chems)} targets")

    per_shell: dict[int, list[dict]] = {}
    for shell, chems in sorted(sampled.items()):
        target_results = []
        for chem in chems:
            r = evaluate_target(hg, chem)
            target_results.append(r)
            print(
                f"  shell {shell} [{chem.name[:40]:40}] "
                f"rxns={r['n_reactions']:5d} chems={r['n_chemicals']:5d} "
                f"max_in={r['max_in_count']:3d} max_out={r['max_out_count']:3d}"
            )
        per_shell[shell] = target_results

    results: dict = {
        "approach": "enzymemap-shell0",
        "shell_cutoff": "unlimited",
        "target_source": "matches" if args.use_matches else "random_sample",
        "seed": SEED,
        "total_samples": (
            sum(len(v) for v in sampled.values()) if args.use_matches else TOTAL_SAMPLES
        ),
        "max_shell": max(sampled) if sampled else 0,
        "sampled_targets": {
            str(shell): [{"id": c.id, "name": c.name} for c in chems]
            for shell, chems in sampled.items()
        },
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
    json_path = RESULTS_DIR / f"enzymemap-shell0{suffix}.json"
    json_path.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {json_path}")


if __name__ == "__main__":
    main()
