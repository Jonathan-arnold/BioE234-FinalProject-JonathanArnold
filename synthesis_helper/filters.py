"""Composable substrate filters for cascade traceback.

Each filter answers: given a reaction and the hypergraph, which of its
substrates should the traceback recurse into? Returning a smaller frozenset
means more chemicals are treated as background (cofactor-like).

Compose multiple filters with ``compose()`` — the result is their
intersection, so each additional filter can only further restrict recursion.
"""

from __future__ import annotations

from typing import Callable

from synthesis_helper.models import Chemical, HyperGraph, Reaction

SubstrateFilter = Callable[[Reaction, HyperGraph], frozenset[Chemical]]


def shell_zero_filter(rxn: Reaction, hg: HyperGraph) -> frozenset[Chemical]:
    """Baseline: recurse into any substrate not in shell 0."""
    return frozenset(s for s in rxn.substrates if hg.chemical_to_shell.get(s) != 0)


def shell_threshold_filter(threshold: int) -> SubstrateFilter:
    """Recurse only into substrates whose shell number exceeds *threshold*.

    threshold=0 is equivalent to shell_zero_filter.
    threshold=1 additionally suppresses chemicals reachable in one step from
    native metabolites — useful for trimming near-universal cofactors that
    weren't in the universal_metabolites list.
    """
    def _filter(rxn: Reaction, hg: HyperGraph) -> frozenset[Chemical]:
        return frozenset(
            s for s in rxn.substrates
            if hg.chemical_to_shell.get(s, 0) > threshold
        )
    return _filter


def _normalize_inchi(inchi: str) -> str:
    """Strip proton (/p) and charge (/q) layers for fuzzy matching.

    Handles the common case where the same cofactor appears in different data
    sources with different protonation or formal-charge states.
    """
    if not inchi.startswith("InChI="):
        return inchi
    parts = inchi.split("/")
    return "/".join(p for p in parts if not p.startswith(("p", "q")))


def inchi_normalized_filter(universal: set[Chemical]) -> SubstrateFilter:
    """Treat a substrate as background if its stripped InChI matches any
    chemical in *universal*, regardless of protonation or charge layer.

    Shell-0 chemicals are always background regardless of InChI match.
    """
    normalized_universal = {_normalize_inchi(c.inchi) for c in universal if c.inchi}

    def _filter(rxn: Reaction, hg: HyperGraph) -> frozenset[Chemical]:
        def _is_background(s: Chemical) -> bool:
            if hg.chemical_to_shell.get(s) == 0:
                return True
            return bool(s.inchi) and _normalize_inchi(s.inchi) in normalized_universal

        return frozenset(s for s in rxn.substrates if not _is_background(s))

    return _filter


def rhea_atom_filter(carrier_map: dict[int, frozenset[str]]) -> SubstrateFilter:
    """Drop substrates identified as carriers by the Rhea atom-tracking pipeline.

    *carrier_map* maps reaction id → frozenset of normalized InChI strings
    (proton/charge layers stripped) for substrates that are carriers in that
    reaction. Build this map with ``scripts/build_rhea_carrier_map.py``.

    For reactions not in the map, falls back to ``shell_zero_filter`` so the
    filter degrades gracefully on unmatched reactions rather than suppressing
    nothing or everything.
    """
    import re as _re

    def _norm(inchi: str) -> str:
        s = _re.sub(r"/p[+-]\d+", "", inchi.strip('"'))
        return _re.sub(r"/q[+-]\d+", "", s)

    def _filter(rxn: Reaction, hg: HyperGraph) -> frozenset[Chemical]:
        carriers = carrier_map.get(rxn.id)
        if carriers is None:
            return shell_zero_filter(rxn, hg)
        return frozenset(
            s for s in rxn.substrates
            if hg.chemical_to_shell.get(s) != 0
            and _norm(s.inchi) not in carriers
        )

    return _filter


def compose(*filters: SubstrateFilter) -> SubstrateFilter:
    """Intersect filters — a substrate must survive all of them to be recursed into."""
    if not filters:
        return lambda rxn, hg: frozenset(rxn.substrates)

    def _composed(rxn: Reaction, hg: HyperGraph) -> frozenset[Chemical]:
        result = filters[0](rxn, hg)
        for f in filters[1:]:
            result &= f(rxn, hg)
        return result

    return _composed
