"""Traceback: build a Cascade for a target chemical from the HyperGraph."""

from __future__ import annotations

from synthesis_helper.models import Cascade, Chemical, HyperGraph, Reaction


def traceback(
    hypergraph: HyperGraph,
    target: Chemical,
    shell_cutoff: int = -1,
    max_producers_per_chemical: int | None = None,
) -> Cascade:
    """Walk backward from target to native metabolites, building a Cascade.

    For each reaction producing the target, recurse on that reaction's
    substrates until shell-0 metabolites are reached.

    Two pruning knobs, applied in order at every chemical:

    1. shell_cutoff (n >= -1) admits a producing reaction only when every
       one of its substrates lies in shell <= (chemical's shell + n).
       n=-1 keeps only producers whose substrates are all strictly closer
       to shell 0 than the chemical itself; raise n to admit more branches.
    2. max_producers_per_chemical caps how many of the surviving producers
       are followed. The lowest-shell producers are kept first (shortest
       routes; ties broken by reaction id). None disables the cap.
       Shell-0 chemicals are unaffected by either knob.
    """
    if shell_cutoff < -1:
        raise ValueError(f"shell_cutoff must be >= -1, got {shell_cutoff}")
    if max_producers_per_chemical is not None and max_producers_per_chemical < 1:
        raise ValueError(
            "max_producers_per_chemical must be None or >= 1, "
            f"got {max_producers_per_chemical}"
        )
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
        max_producers_per_chemical,
        visited_rxns=set(),
        visited_chems=set(),
    )
    return cascade


def _build_producers_index(hg: HyperGraph) -> dict[int, list[Reaction]]:
    """Map chemical id -> list of enabled reactions that produce it.

    Sorted by (reaction shell, reaction id) so the cap can cheaply take the
    first N to favor the shortest routes.
    """
    index: dict[int, list[Reaction]] = {}
    for rxn in hg.reaction_to_shell:
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
    shell_cutoff: int,
    max_producers: int | None,
    visited_rxns: set[int],
    visited_chems: set[int],
) -> None:
    """Recursively collect all reactions that contribute to producing *chemical*."""
    if hg.chemical_to_shell.get(chemical) == 0:
        return  # base case: native/universal metabolite

    if chemical.id in visited_chems:
        return  # break cycles in the reaction graph
    visited_chems.add(chemical.id)

    chem_shell = hg.chemical_to_shell[chemical]
    threshold = chem_shell + shell_cutoff

    surviving: list[Reaction] = [
        rxn
        for rxn in producers_index.get(chemical.id, ())
        if all(hg.chemical_to_shell[s] <= threshold for s in rxn.substrates)
    ]
    if max_producers is not None and len(surviving) > max_producers:
        surviving = surviving[:max_producers]

    for rxn in surviving:
        if rxn.id in visited_rxns:
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
                max_producers,
                visited_rxns,
                visited_chems,
            )
