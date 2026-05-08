"""Universal invariants for synthesize() and traceback() (cascade).

Unlike the small-cell tests (which hard-code expected shell/cascade values for
a specific topology), these tests assert **properties that must hold for any
correctly-implemented synthesize or traceback, on any input**.

The idea: teammate ships a new impl — these invariants catch bugs without
needing to know anything about the impl's internals. Each test is a formal
statement of what the algorithm means.

All tests run against the small_cell fixture (rich enough to exercise most
structural cases) plus a few minimal custom inputs for boundary cases.
"""

from __future__ import annotations

import pytest

from synthesis_helper.models import Chemical, HyperGraph, Reaction
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback


# ---------------------------------------------------------------------------
# synthesize() universal invariants
# ---------------------------------------------------------------------------


def test_invariant_natives_and_universals_all_at_shell_zero(small_cell, small_cell_hg):
    """Every native and every universal metabolite must be present at shell 0.
    This is the seed contract of the BFS — violating it means shell 0 is wrong,
    and every downstream shell is off by at least one."""
    for chem in small_cell.natives:
        assert small_cell_hg.chemical_to_shell.get(chem) == 0, (
            f"native {chem.name} not at shell 0"
        )
    for chem in small_cell.universals:
        assert small_cell_hg.chemical_to_shell.get(chem) == 0, (
            f"universal {chem.name} not at shell 0"
        )


def test_invariant_shells_are_non_negative_integers(small_cell_hg):
    """No negative shells, no fractional shells. Sanity on data types."""
    for chem, shell in small_cell_hg.chemical_to_shell.items():
        assert isinstance(shell, int) and shell >= 0, (
            f"{chem.name} has invalid shell {shell!r}"
        )
    for rxn, shell in small_cell_hg.reaction_to_shell.items():
        assert isinstance(shell, int) and shell >= 1, (
            f"reaction {rxn.id} has invalid shell {shell!r}"
        )


def test_invariant_reaction_shell_equals_max_substrate_shell_plus_one(small_cell_hg):
    """For every enabled reaction: shell(R) == max(shell(s) for s in R.substrates) + 1.

    This is THE defining equation of BFS waveform expansion. Any implementation
    that violates it is wrong — regardless of how the inner loop is structured.
    """
    for rxn, rxn_shell in small_cell_hg.reaction_to_shell.items():
        if not rxn.substrates:
            # Zero-substrate reactions: max of empty set is undefined; skip.
            continue
        sub_shells = [
            small_cell_hg.chemical_to_shell.get(s) for s in rxn.substrates
        ]
        assert all(s is not None for s in sub_shells), (
            f"rxn {rxn.id} has substrate with no shell assignment"
        )
        expected = max(sub_shells) + 1
        assert rxn_shell == expected, (
            f"rxn {rxn.id}: shell {rxn_shell} != max(substrates)+1 = {expected}"
        )


def test_invariant_all_products_of_enabled_reactions_are_reachable(small_cell_hg):
    """If a reaction is enabled (in reaction_to_shell), ALL its products must
    be in chemical_to_shell. A producer exists → its products exist."""
    for rxn in small_cell_hg.reaction_to_shell:
        for prod in rxn.products:
            assert prod in small_cell_hg.chemical_to_shell, (
                f"rxn {rxn.id} is enabled but its product {prod.name} is not reachable"
            )


def test_invariant_all_substrates_of_enabled_reactions_are_reachable(small_cell_hg):
    """Dual of the above: if a reaction is enabled, ALL its substrates are
    reachable too (otherwise the reaction should not have fired)."""
    for rxn in small_cell_hg.reaction_to_shell:
        for sub in rxn.substrates:
            assert sub in small_cell_hg.chemical_to_shell, (
                f"rxn {rxn.id} is enabled but substrate {sub.name} is not reachable"
            )


def test_invariant_every_non_native_chemical_has_a_producer_at_matching_shell(
    small_cell, small_cell_hg
):
    """For every chemical C with shell > 0, at least one enabled reaction must
    produce C and have shell == shell(C). Otherwise we'd have a chemical that
    was magically conjured with no production path — a contradiction.
    """
    for chem, shell in small_cell_hg.chemical_to_shell.items():
        if shell == 0:
            continue
        producers_at_shell = [
            r for r, rshell in small_cell_hg.reaction_to_shell.items()
            if chem in r.products and rshell == shell
        ]
        assert producers_at_shell, (
            f"chemical {chem.name} at shell {shell} has no producer at that shell"
        )


def test_invariant_no_disabled_reaction_has_all_substrates_reachable(small_cell):
    """Completeness / fixed-point invariant: after synthesize returns, there
    must be NO unchosen reaction whose substrates are all reachable. If there
    were, BFS should have fired it. This catches early-exit bugs.
    """
    hg = synthesize(
        small_cell.all_reactions,
        small_cell.natives,
        small_cell.universals,
        verbose=False,
    )
    enabled_ids = {r.id for r in hg.reaction_to_shell}
    for rxn in small_cell.all_reactions:
        if rxn.id in enabled_ids:
            continue
        all_subs_reachable = all(s in hg.chemical_to_shell for s in rxn.substrates)
        assert not all_subs_reachable, (
            f"rxn {rxn.id} has all substrates reachable but was not enabled — "
            "BFS terminated prematurely"
        )


def test_invariant_every_reachable_chemical_can_be_tracebacked(small_cell_hg):
    """For every chemical in chemical_to_shell, traceback() must succeed
    (no exception). This asserts the cascade algorithm is defined for every
    reachable chemical — a cross-module invariant."""
    for chem in small_cell_hg.chemical_to_shell:
        cascade = traceback(small_cell_hg, chem)
        assert cascade.target == chem


# ---------------------------------------------------------------------------
# synthesize() boundary inputs
# ---------------------------------------------------------------------------


def test_synthesize_with_no_reactions_returns_only_natives():
    """Empty reaction corpus: the HyperGraph contains only natives + universals
    at shell 0, and no reactions."""
    a = Chemical(id=1, name="A")
    b = Chemical(id=2, name="B")

    hg = synthesize(
        [], native_metabolites={a, b}, universal_metabolites=set(), verbose=False
    )

    assert hg.chemical_to_shell == {a: 0, b: 0}
    assert hg.reaction_to_shell == {}


def test_synthesize_with_no_natives_or_universals_is_empty():
    """No seeds: no reactions can fire (unless zero-substrate), no chemicals exist."""
    a = Chemical(id=1, name="A")
    b = Chemical(id=2, name="B")
    r = Reaction(id=1, substrates=frozenset([a]), products=frozenset([b]))

    hg = synthesize(
        [r], native_metabolites=set(), universal_metabolites=set(), verbose=False
    )

    assert hg.chemical_to_shell == {}
    assert hg.reaction_to_shell == {}


def test_synthesize_with_only_universals_seeds_shell_zero():
    """Universal-only start: universal metabolites reach shell 0, reactions
    downstream of them fire normally."""
    u = Chemical(id=1, name="U1")
    x = Chemical(id=2, name="X")
    r = Reaction(id=1, substrates=frozenset([u]), products=frozenset([x]))

    hg = synthesize(
        [r], native_metabolites=set(), universal_metabolites={u}, verbose=False
    )

    assert hg.chemical_to_shell[u] == 0
    assert hg.chemical_to_shell[x] == 1


def test_synthesize_native_and_universal_overlap_is_still_shell_zero():
    """If a chemical is in both native_metabolites and universal_metabolites,
    it still ends up at shell 0 exactly once."""
    a = Chemical(id=1, name="A")

    hg = synthesize(
        [], native_metabolites={a}, universal_metabolites={a}, verbose=False
    )

    assert hg.chemical_to_shell == {a: 0}


# ---------------------------------------------------------------------------
# cascade (traceback) universal invariants
# ---------------------------------------------------------------------------


def _all_cascades(hg: HyperGraph) -> dict[Chemical, set]:
    """Helper: traceback every reachable chemical."""
    return {chem: traceback(hg, chem).reactions for chem in hg.chemical_to_shell}


def test_invariant_cascade_reactions_are_subset_of_hypergraph(small_cell_hg):
    """Every reaction in any Cascade must already be an enabled reaction.
    Traceback can never introduce a reaction that BFS didn't enable."""
    for chem, cascade_rxns in _all_cascades(small_cell_hg).items():
        assert cascade_rxns <= set(small_cell_hg.reaction_to_shell), (
            f"cascade for {chem.name} includes reactions not in hypergraph"
        )


def test_invariant_cascade_every_reaction_produces_something_in_transitive_closure(
    small_cell_hg,
):
    """For every reaction R in Cascade(T), R's products must contain either T
    itself or some chemical that is also in the 'trace frontier' (i.e. some
    substrate of another reaction in the cascade).

    In other words: no orphan reactions that don't actually contribute to producing T.
    """
    for target, cascade_rxns in _all_cascades(small_cell_hg).items():
        if not cascade_rxns:
            continue
        frontier: set[Chemical] = {target}
        for r in cascade_rxns:
            frontier.update(r.substrates)
        for r in cascade_rxns:
            produces_useful = any(p in frontier for p in r.products)
            assert produces_useful, (
                f"cascade({target.name}): rxn {r.id} produces nothing "
                "in the substrate/target frontier — orphan reaction"
            )


def test_invariant_cascade_closed_under_substrate_backtracking(small_cell_hg):
    """For every reaction R in Cascade(T), every non-shell-0 substrate of R
    must have at least one producer that is also in Cascade(T).

    This is the formal closure property: cascade is closed under "pull in the
    producers of my substrates". Breaking this means the cascade has a dangling
    non-native substrate with no way to be made — a broken cascade.
    """
    for target, cascade_rxns in _all_cascades(small_cell_hg).items():
        for r in cascade_rxns:
            for sub in r.substrates:
                if small_cell_hg.chemical_to_shell.get(sub) == 0:
                    continue
                has_producer_in_cascade = any(
                    sub in other.products for other in cascade_rxns
                )
                assert has_producer_in_cascade, (
                    f"cascade({target.name}): rxn {r.id} needs substrate "
                    f"{sub.name} but no reaction in the cascade produces it"
                )


def test_invariant_cascade_target_is_produced_or_is_native(small_cell, small_cell_hg):
    """Cascade(T) is either empty (when T is shell-0 native) or contains at
    least one reaction whose products include T."""
    for target in small_cell_hg.chemical_to_shell:
        cascade_rxns = traceback(small_cell_hg, target).reactions
        if small_cell_hg.chemical_to_shell[target] == 0:
            # Native — cascade may or may not be empty (natives can also be
            # products of other reactions), but it's not required to produce T.
            continue
        assert any(target in r.products for r in cascade_rxns), (
            f"cascade({target.name}) contains no reaction producing {target.name}"
        )


def test_invariant_cascade_is_monotonic_under_intermediate_containment(small_cell_hg):
    """Cascade property: if an intermediate M appears in Cascade(T)'s substrate
    frontier (i.e. some reaction in Cascade(T) requires M), then Cascade(M)
    must be a subset of Cascade(T).

    Intuitively: whatever you need to make M when going to T, you already
    have in Cascade(T). Violating this means traceback is inconsistent.
    """
    for target in small_cell_hg.chemical_to_shell:
        if small_cell_hg.chemical_to_shell[target] == 0:
            continue
        cascade_t = traceback(small_cell_hg, target).reactions
        intermediates = set()
        for r in cascade_t:
            for s in r.substrates:
                if small_cell_hg.chemical_to_shell.get(s, 0) > 0:
                    intermediates.add(s)

        for m in intermediates:
            cascade_m = traceback(small_cell_hg, m).reactions
            assert cascade_m <= cascade_t, (
                f"cascade({m.name}) is not a subset of cascade({target.name}) "
                "— traceback is inconsistent"
            )


def test_invariant_cascade_of_shell_zero_target_is_empty(small_cell, small_cell_hg):
    """Every shell-0 chemical has an empty cascade — natives have no upstream
    producers to trace."""
    for chem in small_cell.natives | small_cell.universals:
        cascade = traceback(small_cell_hg, chem)
        assert cascade.reactions == set(), (
            f"shell-0 chemical {chem.name} has non-empty cascade"
        )


def test_invariant_cascade_unreachable_target_raises(small_cell, small_cell_hg):
    """Asking for a cascade on an unreachable chemical must raise, not silently
    return an empty cascade (which would be ambiguous with a native)."""
    # M7 is unreachable in the small cell (R11's substrate X_missing is absent)
    with pytest.raises(ValueError):
        traceback(small_cell_hg, small_cell.chems["M7"])


# ---------------------------------------------------------------------------
# End-to-end integration (parse -> synthesize -> traceback)
# ---------------------------------------------------------------------------


def test_end_to_end_parse_synthesize_traceback(small_cell, tmp_path):
    """Full pipeline: write fixture as TSV, parse, synthesize, traceback —
    confirm the result matches what we get from the in-memory Python fixture.
    Catches integration breakage where individual modules work but their
    interfaces don't compose correctly."""
    from synthesis_helper.parser import (
        parse_chemicals,
        parse_metabolite_list,
        parse_reactions,
    )
    from tests.fixtures.small_cell import write_small_cell_tsvs

    write_small_cell_tsvs(small_cell, tmp_path)

    chems_by_id = parse_chemicals(tmp_path / "enzymemap_chems.tsv")
    reactions = parse_reactions(tmp_path / "enzymemap_reactions.tsv", chems_by_id)
    natives = parse_metabolite_list(tmp_path / "minimal_metabolites.txt", chems_by_id)
    universals = parse_metabolite_list(
        tmp_path / "ubiquitous_metabolites.txt", chems_by_id
    )
    hg_tsv = synthesize(reactions, natives, universals, verbose=False)

    hg_py = synthesize(
        small_cell.all_reactions, small_cell.natives, small_cell.universals, verbose=False
    )

    assert {c.id: s for c, s in hg_tsv.chemical_to_shell.items()} == {
        c.id: s for c, s in hg_py.chemical_to_shell.items()
    }
    assert {r.id: s for r, s in hg_tsv.reaction_to_shell.items()} == {
        r.id: s for r, s in hg_py.reaction_to_shell.items()
    }

    t_tsv = chems_by_id[small_cell.chems["T"].id]
    t_py = small_cell.chems["T"]
    cascade_tsv = traceback(hg_tsv, t_tsv).reactions
    cascade_py = traceback(hg_py, t_py).reactions
    assert {r.id for r in cascade_tsv} == {r.id for r in cascade_py}
