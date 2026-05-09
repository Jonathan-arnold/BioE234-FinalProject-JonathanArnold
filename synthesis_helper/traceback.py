"""Traceback: build a Cascade for a target chemical from the HyperGraph."""

from __future__ import annotations

from synthesis_helper.filters import SubstrateFilter, shell_zero_filter
from synthesis_helper.models import Cascade, Chemical, HyperGraph, Reaction


def traceback(
    hypergraph: HyperGraph,
    target: Chemical,
    max_producers_per_chemical: int | None = None,
    substrate_filter: SubstrateFilter | None = None,
) -> Cascade:
    """Walk backward from target to native metabolites, building a Cascade.

    For each reaction producing the target, recurse on that reaction's
    substrates until shell-0 metabolites are reached.

    max_producers_per_chemical caps how many producing reactions are followed
    per non-shell-0 chemical. When a chemical has more producers than the cap,
    the reactions with the lowest shell are kept (shortest routes first;
    ties broken by reaction id). Shell-0 chemicals are unaffected. None
    disables the cap.

    substrate_filter controls which substrates are recursed into for each
    reaction. Defaults to shell_zero_filter (stop at shell-0 chemicals).
    Use synthesis_helper.filters to build and compose alternative filters.
    """
    if target not in hypergraph.chemical_to_shell:
        raise ValueError(f"Chemical {target.name!r} (id={target.id}) is not reachable.")

    f = substrate_filter if substrate_filter is not None else shell_zero_filter
    producers_index = _build_producers_index(hypergraph)

    cascade = Cascade(target=target)
    _collect_reactions(
        hypergraph,
        target,
        cascade,
        producers_index,
        max_producers_per_chemical,
        f,
        visited_rxns=set(),
        visited_chems=set(),
    )
    return cascade


def _build_producers_index(hg: HyperGraph) -> dict[int, list[Reaction]]:
    """Map chemical id -> list of enabled reactions that produce it.

    Sorted by (reaction shell, reaction id) so callers can cheaply take the
    first N to favor the shortest routes when applying a cap.
    """
    index: dict[int, list[Reaction]] = {}
    for rxn, shell in hg.reaction_to_shell.items():
        for product in rxn.products:
            index.setdefault(product.id, []).append(rxn)
    for producers in index.values():
        producers.sort(key=lambda r: (hg.reaction_to_shell[r], r.id))
    return index


def _collect_reactions(
    hg: HyperGraph,
    chemical: Chemical,
    cascade: Cascade,
    producers_index: dict[int, list[Reaction]],
    max_producers: int | None,
    substrate_filter: SubstrateFilter,
    visited_rxns: set[int],
    visited_chems: set[int],
) -> None:
    """Recursively collect all reactions that contribute to producing *chemical*."""
    if hg.chemical_to_shell.get(chemical) == 0:
        return  # base case: always stop at shell-0 chemicals

    if chemical.id in visited_chems:
        return  # break cycles in the reaction graph
    visited_chems.add(chemical.id)

    producers = producers_index.get(chemical.id, ())
    if max_producers is not None and len(producers) > max_producers:
        producers = producers[:max_producers]

    for rxn in producers:
        if rxn.id in visited_rxns:
            continue
        visited_rxns.add(rxn.id)
        cascade.reactions.add(rxn)
        for substrate in substrate_filter(rxn, hg):
            _collect_reactions(
                hg,
                substrate,
                cascade,
                producers_index,
                max_producers,
                substrate_filter,
                visited_rxns,
                visited_chems,
            )
