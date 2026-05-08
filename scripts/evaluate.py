"""Evaluate the shell-cutoff + max-producers combo across a 2-D sweep.

Builds the hypergraph once, samples reachables per shell (seeded), then for
each (shell_cutoff, max_producers_per_chemical) combination runs traceback
on every sampled target and records cascade-shape metrics per shell.

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
SEED = 20260508
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


def evaluate_target(
    hg: HyperGraph,
    target: Chemical,
    shell_cutoff: int,
    max_producers: int | None,
) -> dict:
    """Run traceback on a single target and record cascade-shape metrics."""
    t0 = time.perf_counter()
    cascade = traceback(
        hg,
        target,
        shell_cutoff=shell_cutoff,
        max_producers_per_chemical=max_producers,
    )
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


def pick_inspection_targets(
    sampled: dict[int, list[Chemical]],
) -> list[tuple[str, Chemical]]:
    """Pick the lowest-shell and highest-shell sampled targets for inspection."""
    shells = sorted(s for s, chems in sampled.items() if chems)
    if not shells:
        return []
    lo = sorted(sampled[shells[0]], key=lambda c: c.id)[0]
    if len(shells) == 1:
        return [(f"shell-{shells[0]}", lo)]
    hi = sorted(sampled[shells[-1]], key=lambda c: c.id)[0]
    return [(f"shell-{shells[0]}", lo), (f"shell-{shells[-1]}", hi)]


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
    caps = list(CAP_SWEEP)
    print(f"Sweeping shell_cutoff over: {cutoffs}")
    print(f"Sweeping max_producers_per_chemical over: {caps}")

    hg, _ = build_hypergraph()

    sampled = sample_targets_per_shell(hg)
    print("\nSampled targets per shell:")
    for shell, chems in sorted(sampled.items()):
        print(f"  shell {shell}: {len(chems)} targets")

    inspection_targets = pick_inspection_targets(sampled)

    results: dict = {
        "approach": "shell-cutoff",
        "seed": SEED,
        "total_samples": TOTAL_SAMPLES,
        "max_shell": max(sampled) if sampled else 0,
        "cap_sweep": caps,
        "sampled_targets": {
            str(shell): [{"id": c.id, "name": c.name} for c in chems]
            for shell, chems in sampled.items()
        },
        "sweep": {},
    }

    inspection_lines: list[str] = []

    for n in cutoffs:
        results["sweep"][str(n)] = {}
        for cap in caps:
            print(f"\n=== shell_cutoff = {n}, max_producers = {cap} ===")
            per_shell: dict[int, list[dict]] = {}
            for shell, chems in sorted(sampled.items()):
                target_results = []
                for chem in chems:
                    r = evaluate_target(hg, chem, shell_cutoff=n, max_producers=cap)
                    target_results.append(r)
                    print(
                        f"  shell {shell} [{chem.name[:40]:40}] "
                        f"rxns={r['n_reactions']:5d} chems={r['n_chemicals']:5d} "
                        f"max_in={r['max_in_count']} (id={r['max_in_chemical_id']}) "
                        f"max_out={r['max_out_count']} (id={r['max_out_chemical_id']})"
                    )
                per_shell[shell] = target_results

            results["sweep"][str(n)][str(cap)] = {
                "per_shell": {
                    str(shell): {
                        "targets": tr,
                        "summary": aggregate(tr),
                    }
                    for shell, tr in per_shell.items()
                },
                "overall": aggregate([r for tr in per_shell.values() for r in tr]),
            }

            # Inspection cascades for picked targets at this (cutoff, cap).
            # for label, target in inspection_targets:
            #     cascade = traceback(
            #         hg,
            #         target,
            #         shell_cutoff=n,
            #         max_producers_per_chemical=cap,
            #     )
            #     nested = render_cascade_nested(hg, cascade)
            #     inspection_lines.append(
            #         f"\n{'=' * 70}\n"
            #         f"shell_cutoff = {n}, max_producers = {cap}, "
            #         f"inspection target = {label}\n"
            #         f"  {target.name} (id={target.id}, "
            #         f"shell={hg.chemical_to_shell[target]})\n"
            #         f"  {len(cascade.reactions)} reactions in cascade\n"
            #         f"{'=' * 70}\n"
            #     )
            #     inspection_lines.append(render_cascade_indented(nested))
            #     inspection_lines.append("")
            #     inspection_lines.append("--- nested dict ---")
            #     inspection_lines.append(json.dumps(nested, indent=2))

    json_path = RESULTS_DIR / "shell-cutoff.json"
    text_path = RESULTS_DIR / "shell-cutoff_inspection.txt"
    json_path.write_text(json.dumps(results, indent=2))
    text_path.write_text("\n".join(inspection_lines))
    print(f"\nWrote {json_path}")
    print(f"Wrote {text_path}")


if __name__ == "__main__":
    main()
