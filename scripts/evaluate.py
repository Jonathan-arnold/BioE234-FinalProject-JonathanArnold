"""Evaluate the shell-cutoff approach across a sweep of cutoff values.

Builds the hypergraph once, samples 5 reachables per shell (seeded), then for
each n in [-1, max_cutoff] runs traceback + enumerate_pathways on every
sampled target and records per-shell metrics. Also renders one shell-1 and
one shell-2 target's cascade as a nested dict for manual inspection.

Usage:
    python scripts/evaluate.py --max-cutoff 2
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
import time
from collections import defaultdict
from pathlib import Path

from synthesis_helper.composition import pathway_to_composition  # noqa: F401
from synthesis_helper.models import Cascade, Chemical, HyperGraph, Pathway, Reaction
from synthesis_helper.parser import (
    parse_chemicals,
    parse_metabolite_list,
    parse_reactions,
)
from synthesis_helper.pathways import enumerate_pathways
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "eval_results"
SEED = 20260507
SAMPLES_PER_SHELL = 1
MAX_PATHWAYS = 100_000


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
    """Group reachables by shell (excluding shell 0), sample 5 per shell."""
    by_shell: dict[int, list[Chemical]] = defaultdict(list)
    for chem, shell in hg.chemical_to_shell.items():
        if shell > 0:
            by_shell[shell].append(chem)

    rng = random.Random(SEED)
    sampled: dict[int, list[Chemical]] = {}
    for shell in sorted(by_shell):
        pool = sorted(by_shell[shell], key=lambda c: c.id)  # determinism
        k = min(SAMPLES_PER_SHELL, len(pool))
        sampled[shell] = rng.sample(pool, k)
    return sampled


def pathway_depth(pathway: Pathway, hg: HyperGraph) -> int:
    """Longest chain from a shell-0 leaf to the target through this pathway.

    `len(pathway.reactions)` counts every node in the pathway DAG (each
    branch of a multi-substrate reaction adds nodes), so it overstates
    "steps." Depth is what most people mean by pathway length.
    """
    producers_in_pathway: dict[int, list[Reaction]] = {}
    for rxn in pathway.reactions:
        for prod in rxn.products:
            producers_in_pathway.setdefault(prod.id, []).append(rxn)

    rxn_depth: dict[int, int] = {}

    def depth_of_rxn(rxn: Reaction) -> int:
        if rxn.id in rxn_depth:
            return rxn_depth[rxn.id]
        rxn_depth[rxn.id] = 0  # cycle guard
        max_sub = 0
        for s in rxn.substrates:
            if hg.chemical_to_shell.get(s) == 0:
                continue
            sub_rxns = producers_in_pathway.get(s.id, [])
            for sr in sub_rxns:
                d = depth_of_rxn(sr)
                if d > max_sub:
                    max_sub = d
        rxn_depth[rxn.id] = 1 + max_sub
        return rxn_depth[rxn.id]

    target_rxns = producers_in_pathway.get(pathway.target.id, [])
    if not target_rxns:
        return 0
    return max(depth_of_rxn(r) for r in target_rxns)


def evaluate_target(hg: HyperGraph, target: Chemical, shell_cutoff: int) -> dict:
    """Run traceback + pathway enumeration on a single target."""
    t0 = time.perf_counter()
    cascade = traceback(hg, target, shell_cutoff=shell_cutoff)
    t_cascade = time.perf_counter() - t0

    t0 = time.perf_counter()
    pathways = enumerate_pathways(cascade, hg, max_pathways=MAX_PATHWAYS)
    t_pathways = time.perf_counter() - t0
    hit_cap = len(pathways) >= MAX_PATHWAYS

    depths = [pathway_depth(p, hg) for p in pathways]
    sizes = [len(p.reactions) for p in pathways]
    return {
        "chemical_id": target.id,
        "chemical_name": target.name,
        "target_shell": hg.chemical_to_shell[target],
        "n_reactions": len(cascade.reactions),
        "n_pathways": len(pathways),
        "hit_pathway_cap": hit_cap,
        "min_pathway_depth": min(depths) if depths else None,
        "mean_pathway_depth": statistics.mean(depths) if depths else None,
        "max_pathway_depth": max(depths) if depths else None,
        "mean_pathway_size": statistics.mean(sizes) if sizes else None,
        "max_pathway_size": max(sizes) if sizes else None,
        "cascade_seconds": t_cascade,
        "pathways_seconds": t_pathways,
    }


def aggregate(target_results: list[dict]) -> dict:
    """Summary stats across a list of per-target results."""
    if not target_results:
        return {"count": 0}
    n_rxns = [r["n_reactions"] for r in target_results]
    n_paths = [r["n_pathways"] for r in target_results]
    mean_depths = [
        r["mean_pathway_depth"]
        for r in target_results
        if r["mean_pathway_depth"] is not None
    ]
    max_depths = [
        r["max_pathway_depth"]
        for r in target_results
        if r["max_pathway_depth"] is not None
    ]
    coverage = sum(1 for r in target_results if r["n_pathways"] > 0)
    capped = sum(1 for r in target_results if r.get("hit_pathway_cap"))
    return {
        "count": len(target_results),
        "coverage": coverage / len(target_results),
        "frac_hit_pathway_cap": capped / len(target_results),
        "mean_n_reactions": statistics.mean(n_rxns),
        "median_n_reactions": statistics.median(n_rxns),
        "mean_n_pathways": statistics.mean(n_paths),
        "median_n_pathways": statistics.median(n_paths),
        "mean_pathway_depth_overall": (
            statistics.mean(mean_depths) if mean_depths else None
        ),
        "max_pathway_depth_overall": (max(max_depths) if max_depths else None),
    }


def render_cascade_nested(
    hg: HyperGraph,
    cascade: Cascade,
) -> dict:
    """Render a cascade as a nested dict for human inspection.

    Cycles are broken: a chemical already expanded higher in the current
    branch is rendered as the string "<seen above>".
    """
    producers: dict[int, list[Reaction]] = {}
    for rxn in cascade.reactions:
        for prod in rxn.products:
            producers.setdefault(prod.id, []).append(rxn)
    for rxns in producers.values():
        rxns.sort(key=lambda r: r.id)

    def chem_label(c: Chemical) -> str:
        return f"[shell {hg.chemical_to_shell.get(c, '?')}] {c.name} (id={c.id})"

    def expand_chem(chem: Chemical, in_path: set[int]) -> object:
        shell = hg.chemical_to_shell.get(chem)
        if shell == 0:
            return "<shell 0>"
        if chem.id in in_path:
            return "<seen above>"
        rxns = producers.get(chem.id, [])
        if not rxns:
            return "<no producer in cascade>"
        in_path = in_path | {chem.id}
        return {
            f"rxn {r.id}"
            + (f" (EC {r.ecnum})" if r.ecnum else ""): {
                chem_label(s): expand_chem(s, in_path)
                for s in sorted(r.substrates, key=lambda x: x.id)
            }
            for r in rxns
        }

    return {chem_label(cascade.target): expand_chem(cascade.target, set())}


def render_cascade_indented(nested: dict, indent: int = 0) -> str:
    """Flat indented form of the nested-dict cascade."""
    lines: list[str] = []
    pad = "  " * indent
    if isinstance(nested, str):
        lines.append(f"{pad}{nested}")
        return "\n".join(lines)
    for key, val in nested.items():
        lines.append(f"{pad}{key}")
        lines.append(render_cascade_indented(val, indent + 1))
    return "\n".join(lines)


def pick_inspection_target(
    sampled: dict[int, list[Chemical]], shell: int
) -> Chemical | None:
    """Pick a deterministic target from the sample for manual inspection."""
    if shell not in sampled or not sampled[shell]:
        return None
    return sorted(sampled[shell], key=lambda c: c.id)[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--max-cutoff",
        type=int,
        default=2,
        help="Sweep shell_cutoff from -1 up to this value (inclusive).",
    )
    args = ap.parse_args()

    if args.max_cutoff < -1:
        ap.error("--max-cutoff must be >= -1")

    RESULTS_DIR.mkdir(exist_ok=True)
    cutoffs = list(range(-1, args.max_cutoff + 1))
    print(f"Sweeping shell_cutoff over: {cutoffs}")

    hg, _ = build_hypergraph()

    sampled = sample_targets_per_shell(hg)
    print("\nSampled targets per shell:")
    for shell, chems in sorted(sampled.items()):
        print(f"  shell {shell}: {len(chems)} targets")

    inspection_shell1 = pick_inspection_target(sampled, 1)
    inspection_shell2 = pick_inspection_target(sampled, 2)

    results: dict = {
        "approach": "shell-cutoff",
        "seed": SEED,
        "samples_per_shell": SAMPLES_PER_SHELL,
        "max_shell": max(sampled) if sampled else 0,
        "sampled_targets": {
            str(shell): [{"id": c.id, "name": c.name} for c in chems]
            for shell, chems in sampled.items()
        },
        "sweep": {},
    }

    inspection_lines: list[str] = []

    for n in cutoffs:
        print(f"\n=== shell_cutoff = {n} ===")
        per_shell: dict[int, list[dict]] = {}
        for shell, chems in sorted(sampled.items()):
            target_results = []
            for chem in chems:
                r = evaluate_target(hg, chem, shell_cutoff=n)
                target_results.append(r)
                print(
                    f"  shell {shell} [{chem.name[:40]:40}] "
                    f"rxns={r['n_reactions']:5d} paths={r['n_pathways']:6d} "
                    f"depth(mean/max)={r['mean_pathway_depth']}/{r['max_pathway_depth']}"
                )
            per_shell[shell] = target_results

        results["sweep"][str(n)] = {
            "per_shell": {
                str(shell): {
                    "targets": tr,
                    "summary": aggregate(tr),
                }
                for shell, tr in per_shell.items()
            },
            "overall": aggregate([r for tr in per_shell.values() for r in tr]),
        }

        # Inspection cascades for shell 1 + shell 2.
        for label, target in (
            ("shell-1", inspection_shell1),
            ("shell-2", inspection_shell2),
        ):
            if target is None:
                continue
            cascade = traceback(hg, target, shell_cutoff=n)
            nested = render_cascade_nested(hg, cascade)
            inspection_lines.append(
                f"\n{'=' * 70}\n"
                f"shell_cutoff = {n}, inspection target = {label}\n"
                f"  {target.name} (id={target.id}, "
                f"shell={hg.chemical_to_shell[target]})\n"
                f"  {len(cascade.reactions)} reactions in cascade\n"
                f"{'=' * 70}\n"
            )
            inspection_lines.append(render_cascade_indented(nested))
            inspection_lines.append("")
            inspection_lines.append("--- nested dict ---")
            inspection_lines.append(json.dumps(nested, indent=2))

    json_path = RESULTS_DIR / "shell-cutoff.json"
    text_path = RESULTS_DIR / "shell-cutoff_inspection.txt"
    json_path.write_text(json.dumps(results, indent=2))
    text_path.write_text("\n".join(inspection_lines))
    print(f"\nWrote {json_path}")
    print(f"Wrote {text_path}")


if __name__ == "__main__":
    main()
