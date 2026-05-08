"""Traceback: build a Cascade for a target chemical from the HyperGraph."""

from __future__ import annotations

from synthesis_helper.models import Cascade, Chemical, HyperGraph, Reaction


def traceback(
    hypergraph: HyperGraph,
    target: Chemical,
    shell_cutoff: int | None = None,
) -> Cascade:
    """Walk backward from target to native metabolites, building a Cascade.

    For each reaction producing the target, recurse on that reaction's
    substrates until shell-0 metabolites are reached.

    shell_cutoff (n >= -1) admits a producing reaction only when every one
    of its substrates lies in shell <= (chemical's shell + n). Tighter
    values prune more aggressively: n=-1 keeps only producers whose
    substrates are all strictly closer to shell 0 than the chemical itself.
    Pass ``None`` (the default) to disable the cutoff entirely — every
    producing reaction is admitted regardless of substrate shell.
    """
    if shell_cutoff is not None and shell_cutoff < -1:
        raise ValueError(f"shell_cutoff must be >= -1 or None, got {shell_cutoff}")
    if target not in hypergraph.chemical_to_shell:
        raise ValueError(f"Chemical {target.name!r} (id={target.id}) is not reachable.")

    producers_index = _build_producers_index(hypergraph)

    cascade = Cascade(target=target)
    _collect_reactions(
        hypergraph,
        target,
        cascade,
        producers_index,
        shell_cutoff,
        visited_rxns=set(),
        visited_chems=set(),
    )
    return cascade


def _build_producers_index(hg: HyperGraph) -> dict[int, list[Reaction]]:
    """Map chemical id -> list of enabled reactions that produce it, by id."""
    index: dict[int, list[Reaction]] = {}
    for rxn in hg.reaction_to_shell:
        for product in rxn.products:
            index.setdefault(product.id, []).append(rxn)
    for producers in index.values():
        producers.sort(key=lambda r: r.id)
    return index


def _collect_reactions(
    hg: HyperGraph,
    chemical: Chemical,
    cascade: Cascade,
    producers_index: dict[int, list[Reaction]],
    shell_cutoff: int | None,
    visited_rxns: set[int],
    visited_chems: set[int],
) -> None:
    """Recursively collect all reactions that contribute to producing *chemical*."""
    if hg.chemical_to_shell.get(chemical) == 0:
        return  # base case: native/universal metabolite

    if chemical.id in visited_chems:
        return  # break cycles in the reaction graph
    visited_chems.add(chemical.id)

    threshold: int | None
    if shell_cutoff is None:
        threshold = None
    else:
        threshold = hg.chemical_to_shell[chemical] + shell_cutoff

    for rxn in producers_index.get(chemical.id, ()):
        if rxn.id in visited_rxns:
            continue
        if threshold is not None and any(
            hg.chemical_to_shell[s] > threshold for s in rxn.substrates
        ):
            continue
        visited_rxns.add(rxn.id)
        cascade.reactions.add(rxn)
        for substrate in rxn.substrates:
            _collect_reactions(
                hg,
                substrate,
                cascade,
                producers_index,
                shell_cutoff,
                visited_rxns,
                visited_chems,
            )
