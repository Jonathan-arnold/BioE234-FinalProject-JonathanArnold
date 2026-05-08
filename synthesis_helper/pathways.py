"""Enumerate individual Pathways from a Cascade.

A pathway is a set of reactions that produces the cascade's target from
shell-0 (native) chemicals. Enumeration uses a shared-frontier choice-function
model: a partial pathway is a map (chemical -> chosen producer reaction) plus
a frontier of chemicals still needing a producer. Distinct choice-functions
yield distinct reaction sets (modulo multi-product overlap, which is deduped
at yield time), so there is no cartesian-product explosion to collapse.
"""

from __future__ import annotations

from typing import Iterator

from synthesis_helper.models import Cascade, Chemical, HyperGraph, Pathway, Reaction


def _topological_order(
    rxn_set: set[Reaction], hypergraph: HyperGraph
) -> list[Reaction] | None:
    """Order reactions so each fires only after its substrates are produced.

    Returns None if some reaction in the set is not forward-reachable from
    shell-0 chemicals through the rest of the set — the canonical
    cycle-disconnected-from-natives case.
    """
    produced: set[Chemical] = set()
    remaining = list(rxn_set)
    order: list[Reaction] = []
    while remaining:
        ready = [
            r
            for r in remaining
            if all(
                hypergraph.chemical_to_shell.get(s) == 0 or s in produced
                for s in r.substrates
            )
        ]
        if not ready:
            return None
        ready.sort(key=lambda r: r.id)
        for r in ready:
            order.append(r)
            produced.update(r.products)
        ready_ids = {r.id for r in ready}
        remaining = [r for r in remaining if r.id not in ready_ids]
    return order


def enumerate_pathways(
    cascade: Cascade,
    hypergraph: HyperGraph,
    max_pathways: int = 1000,
) -> list[Pathway]:
    """Enumerate distinct pathways implied by a cascade.

    Walks a choice-function search tree: at each step, pop the lowest-id
    chemical from the frontier and branch over its producer reactions. Each
    chosen producer adds its non-native substrates that aren't already decided
    to the frontier. A complete pathway is yielded when the frontier is
    empty.

    Pathways disconnected from shell-0 chemicals (cycle-only reaction sets)
    are dropped via the topological-order check. Multi-product reactions can
    occasionally produce two distinct choice-functions that map to the same
    reaction set; those are deduped at yield time via a frozenset key.
    """
    target = cascade.target

    if hypergraph.chemical_to_shell.get(target) == 0:
        return [Pathway(target=target, reactions=[], metabolites={target})]

    producers: dict[Chemical, list[Reaction]] = {}
    for rxn in cascade.reactions:
        for product in rxn.products:
            producers.setdefault(product, []).append(rxn)
    for rxns in producers.values():
        rxns.sort(key=lambda r: r.id)

    seen_keys: set[frozenset[int]] = set()

    def walk(
        choices: dict[Chemical, Reaction], frontier: frozenset[Chemical]
    ) -> Iterator[set[Reaction]]:
        if not frontier:
            rxn_set = set(choices.values())
            key = frozenset(r.id for r in rxn_set)
            if key in seen_keys:
                return
            seen_keys.add(key)
            yield rxn_set
            return

        chem = min(frontier, key=lambda c: c.id)
        rest = frontier - {chem}
        for rxn in producers.get(chem, ()):
            new_choices = {**choices, chem: rxn}
            needs = frozenset(
                s
                for s in rxn.substrates
                if hypergraph.chemical_to_shell.get(s) != 0
                and s not in new_choices
            )
            yield from walk(new_choices, rest | needs)

    results: list[Pathway] = []
    for rxn_set in walk({}, frozenset({target})):
        ordered = _topological_order(rxn_set, hypergraph)
        if ordered is None:
            continue
        metabolites: set[Chemical] = set()
        for r in ordered:
            metabolites.update(r.substrates)
            metabolites.update(r.products)
        results.append(
            Pathway(target=target, reactions=ordered, metabolites=metabolites)
        )
        if len(results) >= max_pathways:
            break

    return results
