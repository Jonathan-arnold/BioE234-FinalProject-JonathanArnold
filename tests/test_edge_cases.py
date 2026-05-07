"""Edge case tests beyond the primary small-cell fixture.

Each test uses a minimal inline fixture to isolate one specific edge behavior.
Grouped by category; see validation/small_cell.md for the index linking these
tests back to the documented gap list.
"""

from __future__ import annotations

from synthesis_helper.models import Chemical, Reaction
from synthesis_helper.parser import (
    parse_chemicals,
    parse_metabolite_list,
    parse_reactions,
)
from synthesis_helper.pathways import enumerate_pathways
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback


# ---------------------------------------------------------------------------
# Synthesize — structural edge cases
# ---------------------------------------------------------------------------


def test_self_referential_reaction_does_not_loop():
    """R: {A} -> {A}. Matches real MetaCyc pattern (e.g. rxn 88: `4 263 -> 4 263`).
    BFS must terminate; A stays at shell 0; no phantom chemicals appear."""
    a = Chemical(id=1, name="A")
    r = Reaction(id=1, substrates=frozenset([a]), products=frozenset([a]))

    hg = synthesize([r], native_metabolites={a}, universal_metabolites=set(), verbose=False)

    assert hg.chemical_to_shell[a] == 0
    assert hg.reaction_to_shell[r] == 1
    assert len(hg.chemical_to_shell) == 1


def test_product_already_in_shell_zero_is_not_overwritten():
    """R: {A} -> {N} where N is native. N's shell must stay 0, not get bumped."""
    n = Chemical(id=1, name="N")
    a = Chemical(id=2, name="A")
    r_make_a = Reaction(id=1, substrates=frozenset([n]), products=frozenset([a]))
    r_back = Reaction(id=2, substrates=frozenset([a]), products=frozenset([n]))

    hg = synthesize(
        [r_make_a, r_back], native_metabolites={n}, universal_metabolites=set(), verbose=False
    )

    assert hg.chemical_to_shell[n] == 0
    assert hg.chemical_to_shell[a] == 1
    assert hg.reaction_to_shell[r_make_a] == 1
    assert hg.reaction_to_shell[r_back] == 2


def test_multi_product_reaction_both_products_reachable():
    """R: {A, B} -> {C, D}. Both C and D get shell 1; traceback from either returns R."""
    a = Chemical(id=1, name="A")
    b = Chemical(id=2, name="B")
    c = Chemical(id=3, name="C")
    d = Chemical(id=4, name="D")
    r = Reaction(id=1, substrates=frozenset([a, b]), products=frozenset([c, d]))

    hg = synthesize([r], native_metabolites={a, b}, universal_metabolites=set(), verbose=False)

    assert hg.chemical_to_shell[c] == 1
    assert hg.chemical_to_shell[d] == 1
    assert traceback(hg, c).reactions == {r}
    assert traceback(hg, d).reactions == {r}


def test_synthesize_deep_chain_20_shells():
    """20-step linear chain. Verifies BFS terminates cleanly on deeper topologies
    than the small cell (max shell 3) exercises."""
    chems = [Chemical(id=i, name=f"C{i}") for i in range(21)]
    reactions = [
        Reaction(id=i, substrates=frozenset([chems[i - 1]]), products=frozenset([chems[i]]))
        for i in range(1, 21)
    ]

    hg = synthesize(
        reactions, native_metabolites={chems[0]}, universal_metabolites=set(), verbose=False
    )

    for i in range(21):
        assert hg.chemical_to_shell[chems[i]] == i
    assert len(hg.reaction_to_shell) == 20


def test_empty_substrate_reaction_enables_at_shell_one():
    """R: {} -> {A}. With no substrates, all([]) is True, so reaction enables
    immediately at shell 1. This is the current behavior; flagged here so any
    change to reject zero-substrate reactions breaks this test explicitly."""
    a = Chemical(id=1, name="A")
    r = Reaction(id=1, substrates=frozenset(), products=frozenset([a]))

    hg = synthesize(
        [r], native_metabolites=set(), universal_metabolites=set(), verbose=False
    )

    assert hg.chemical_to_shell[a] == 1
    assert hg.reaction_to_shell[r] == 1


def test_empty_product_reaction_is_noop():
    """R: {A} -> {}. Reaction enables but adds no chemicals."""
    a = Chemical(id=1, name="A")
    r = Reaction(id=1, substrates=frozenset([a]), products=frozenset())

    hg = synthesize(
        [r], native_metabolites={a}, universal_metabolites=set(), verbose=False
    )

    assert hg.reaction_to_shell[r] == 1
    assert hg.chemical_to_shell == {a: 0}


def test_synthesize_is_deterministic():
    """Running synthesize twice on identical inputs produces identical shell maps."""
    from tests.fixtures.small_cell import build_small_cell

    sc = build_small_cell()
    hg1 = synthesize(sc.all_reactions, sc.natives, sc.universals, verbose=False)
    hg2 = synthesize(sc.all_reactions, sc.natives, sc.universals, verbose=False)

    shells_1 = {c.id: s for c, s in hg1.chemical_to_shell.items()}
    shells_2 = {c.id: s for c, s in hg2.chemical_to_shell.items()}
    rxns_1 = {r.id: s for r, s in hg1.reaction_to_shell.items()}
    rxns_2 = {r.id: s for r, s in hg2.reaction_to_shell.items()}
    assert shells_1 == shells_2
    assert rxns_1 == rxns_2


def test_synthesize_handles_100_chemical_corpus():
    """Smoke/scale test: 100-chemical linear chain completes and gives correct shell."""
    chems = [Chemical(id=i, name=f"C{i}") for i in range(100)]
    reactions = [
        Reaction(id=i, substrates=frozenset([chems[i - 1]]), products=frozenset([chems[i]]))
        for i in range(1, 100)
    ]

    hg = synthesize(
        reactions, native_metabolites={chems[0]}, universal_metabolites=set(), verbose=False
    )

    assert hg.chemical_to_shell[chems[99]] == 99
    assert len(hg.reaction_to_shell) == 99


# ---------------------------------------------------------------------------
# Traceback / enumerate_pathways — edge cases
# ---------------------------------------------------------------------------


def test_enumerate_pathways_handles_chemical_level_cycle():
    """Cascade contains a chemical-level cycle: B <-> C via different reactions.
    enumerate_pathways must terminate (the cascade cycle is broken by current_rxns
    membership check — revisiting the same Reaction is rejected)."""
    a = Chemical(id=1, name="A")
    b = Chemical(id=2, name="B")
    c = Chemical(id=3, name="C")
    r_ab = Reaction(id=1, substrates=frozenset([a]), products=frozenset([b]))
    r_bc = Reaction(id=2, substrates=frozenset([b]), products=frozenset([c]))
    r_cb = Reaction(id=3, substrates=frozenset([c]), products=frozenset([b]))

    hg = synthesize(
        [r_ab, r_bc, r_cb], native_metabolites={a}, universal_metabolites=set(), verbose=False
    )
    # shell_cutoff=1 admits the C->B back-edge so the cycle is fully captured;
    # the default (-1) would prune r_cb since C sits at the same shell as the
    # B it produces. The test target is enumerate_pathways termination.
    cascade = traceback(hg, c, shell_cutoff=1)

    # Cascade should contain all three reactions; r_cb is pulled in because b is a
    # non-shell-0 substrate of r_bc and r_cb is a producer of b.
    assert cascade.reactions == {r_ab, r_bc, r_cb}

    # enumerate_pathways must NOT loop. At minimum, the direct path A -> B -> C exists.
    pathways = enumerate_pathways(cascade, hg)
    reaction_sets = {frozenset(p.reactions) for p in pathways}
    assert frozenset({r_ab, r_bc}) in reaction_sets


def test_enumerate_pathways_respects_max_pathways_cap_in_deep_cascade():
    """max_pathways bounds pathway enumeration for recursive (deep) cascades.

    Note: the cap is checked at each _backtrack entry. Very wide, shallow cascades
    may slightly exceed the cap; deep ones (as here) respect it.
    """
    natives = [Chemical(id=i, name=f"N{i}") for i in range(5)]
    intermediates = [Chemical(id=10 + i, name=f"I{i}") for i in range(5)]
    t = Chemical(id=100, name="T")

    rxns: list[Reaction] = []
    for i in range(5):
        rxns.append(
            Reaction(id=i, substrates=frozenset([natives[i]]), products=frozenset([intermediates[i]]))
        )
        rxns.append(
            Reaction(
                id=100 + i,
                substrates=frozenset([intermediates[i]]),
                products=frozenset([t]),
            )
        )

    hg = synthesize(rxns, native_metabolites=set(natives), universal_metabolites=set(), verbose=False)
    cascade = traceback(hg, t)

    pathways = enumerate_pathways(cascade, hg, max_pathways=3)
    assert len(pathways) <= 3
    assert len(pathways) >= 1  # at least one route enumerated before the cap fired


# ---------------------------------------------------------------------------
# Cascade algorithm (traceback) — invariants
# ---------------------------------------------------------------------------


def test_traceback_on_native_chemical_returns_empty_cascade():
    """Target is a shell-0 native. Natives have no producers to trace — the
    cascade is correctly empty."""
    a = Chemical(id=1, name="A")
    b = Chemical(id=2, name="B")
    r = Reaction(id=1, substrates=frozenset([a]), products=frozenset([b]))

    hg = synthesize([r], native_metabolites={a}, universal_metabolites=set(), verbose=False)
    cascade = traceback(hg, a)

    assert cascade.target == a
    assert cascade.reactions == set()


def test_traceback_deep_cascade_includes_full_chain():
    """8-step chain. Cascade(C8) must include all 8 reactions — full transitive closure."""
    chems = [Chemical(id=i, name=f"C{i}") for i in range(9)]
    reactions = [
        Reaction(id=i, substrates=frozenset([chems[i - 1]]), products=frozenset([chems[i]]))
        for i in range(1, 9)
    ]

    hg = synthesize(
        reactions, native_metabolites={chems[0]}, universal_metabolites=set(), verbose=False
    )
    cascade = traceback(hg, chems[8])

    assert cascade.reactions == set(reactions)


def test_traceback_diamond_dependency_deduplicates_shared_ancestor():
    """Diamond topology:
         R1: A -> B
         R2: B -> C
         R3: B -> D
         R4: C + D -> T
       B's producer R1 is reachable via both the C-branch and D-branch. Because
       Cascade.reactions is a set, R1 must appear exactly once (not twice)."""
    a = Chemical(id=1, name="A")
    b = Chemical(id=2, name="B")
    c = Chemical(id=3, name="C")
    d = Chemical(id=4, name="D")
    t = Chemical(id=5, name="T")

    r1 = Reaction(id=1, substrates=frozenset([a]), products=frozenset([b]))
    r2 = Reaction(id=2, substrates=frozenset([b]), products=frozenset([c]))
    r3 = Reaction(id=3, substrates=frozenset([b]), products=frozenset([d]))
    r4 = Reaction(id=4, substrates=frozenset([c, d]), products=frozenset([t]))

    hg = synthesize(
        [r1, r2, r3, r4], native_metabolites={a}, universal_metabolites=set(), verbose=False
    )
    cascade = traceback(hg, t)

    assert cascade.reactions == {r1, r2, r3, r4}


def test_traceback_shared_substrate_two_producers_of_target():
    """Two different reactions produce T from the same upstream chemical X:
         R1: A -> X
         R2: X -> T
         R3: X -> T   (different Reaction, same X -> T chemistry)
       Cascade(T) must include all three; R1 is traced once despite being
       reachable from both R2 and R3."""
    a = Chemical(id=1, name="A")
    x = Chemical(id=2, name="X")
    t = Chemical(id=3, name="T")

    r1 = Reaction(id=1, substrates=frozenset([a]), products=frozenset([x]))
    r2 = Reaction(id=2, substrates=frozenset([x]), products=frozenset([t]))
    r3 = Reaction(id=3, substrates=frozenset([x]), products=frozenset([t]))

    hg = synthesize(
        [r1, r2, r3], native_metabolites={a}, universal_metabolites=set(), verbose=False
    )
    cascade = traceback(hg, t)

    assert cascade.reactions == {r1, r2, r3}


def test_traceback_excludes_disjoint_subgraph():
    """Two independent sub-pathways both start from the same native A:
         R_main_1: A -> M_main,  R_main_2: M_main -> T_main
         R_other_1: A -> M_other, R_other_2: M_other -> T_other
       Cascade(T_main) must include only the main chain; traceback must NOT
       cross over into T_other's reactions just because they share ancestor A."""
    a = Chemical(id=1, name="A")
    m_main = Chemical(id=2, name="M_main")
    t_main = Chemical(id=3, name="T_main")
    m_other = Chemical(id=4, name="M_other")
    t_other = Chemical(id=5, name="T_other")

    r_main_1 = Reaction(id=1, substrates=frozenset([a]), products=frozenset([m_main]))
    r_main_2 = Reaction(id=2, substrates=frozenset([m_main]), products=frozenset([t_main]))
    r_other_1 = Reaction(id=3, substrates=frozenset([a]), products=frozenset([m_other]))
    r_other_2 = Reaction(id=4, substrates=frozenset([m_other]), products=frozenset([t_other]))

    hg = synthesize(
        [r_main_1, r_main_2, r_other_1, r_other_2],
        native_metabolites={a},
        universal_metabolites=set(),
        verbose=False,
    )
    cascade = traceback(hg, t_main)

    assert cascade.reactions == {r_main_1, r_main_2}
    assert r_other_1 not in cascade.reactions
    assert r_other_2 not in cascade.reactions


def test_cascade_is_superset_of_all_pathway_reactions():
    """Invariant across cascade and enumerate_pathways: every reaction in every
    pathway must already be in the cascade. Currently true by construction
    (enumerate_pathways only draws from cascade.reactions) — pinning it guards
    against future refactors that might pull from the full HyperGraph."""
    from tests.fixtures.small_cell import build_small_cell

    sc = build_small_cell()
    hg = synthesize(sc.all_reactions, sc.natives, sc.universals, verbose=False)
    cascade = traceback(hg, sc.chems["T"])
    pathways = enumerate_pathways(cascade, hg)

    union_of_pathway_reactions = set()
    for p in pathways:
        union_of_pathway_reactions.update(p.reactions)

    assert union_of_pathway_reactions <= cascade.reactions, (
        f"pathway reactions leaked outside cascade: "
        f"{union_of_pathway_reactions - cascade.reactions}"
    )


# ---------------------------------------------------------------------------
# Model invariants
# ---------------------------------------------------------------------------


def test_same_id_different_inchi_treated_as_same_chemical():
    """Chemical equality is by id only. Two objects with same id collide in sets
    regardless of InChI/name differences. Documents the model invariant:
    Chemical.id must be canonical within a given corpus."""
    a1 = Chemical(id=1, name="A", inchi="InChI=X")
    a2 = Chemical(id=1, name="B", inchi="InChI=Y")

    assert a1 == a2
    assert hash(a1) == hash(a2)
    assert {a1, a2} == {a1}


# ---------------------------------------------------------------------------
# Parser — edge cases
# ---------------------------------------------------------------------------


def test_parse_reactions_silently_drops_unknown_substrate_ids(tmp_path):
    """parse_reactions silently drops ids not in the chemicals dict. This is a
    known data-quality risk: a reaction referencing a missing chemical ends up
    with a smaller substrate set than declared, potentially letting BFS enable
    it when it shouldn't.

    Test pins current behavior; if we later decide to warn or raise on unknown
    ids, this assertion will break loudly and force a conscious choice.
    """
    (tmp_path / "good_chems.txt").write_text(
        "id\tname\tinchi\tsmiles\n"
        "1\tA\tInChI=1S/A\tnull\n"
    )
    (tmp_path / "good_reactions.txt").write_text(
        "rxnid\tecnum\tsubstrates\tproducts\n"
        "1\t\t1 99999\t1\n"
    )

    chems = parse_chemicals(tmp_path / "good_chems.txt")
    reactions = parse_reactions(tmp_path / "good_reactions.txt", chems)

    assert len(reactions) == 1
    assert len(reactions[0].substrates) == 1
    assert chems[1] in reactions[0].substrates


def test_parse_metabolite_list_loose_match_adds_all_stereoisomers(tmp_path):
    """A metabolite declared without stereo layers loose-matches BOTH stereoisomers
    in the chemicals corpus. This is the behavior documented in parser.py:94-97
    for handling cofactors like SAM."""
    (tmp_path / "good_chems.txt").write_text(
        "id\tname\tinchi\tsmiles\n"
        "1\tchem-R\tInChI=1S/C6H12/t1+/m0/s1\tnull\n"
        "2\tchem-S\tInChI=1S/C6H12/t1-/m1/s1\tnull\n"
    )
    (tmp_path / "metabolites.txt").write_text(
        "name\tinchi\tdescriptor\n"
        "no-stereo\tInChI=1S/C6H12\tfeedstock\n"
    )

    chems = parse_chemicals(tmp_path / "good_chems.txt")
    metabolites = parse_metabolite_list(tmp_path / "metabolites.txt", chems)

    assert {c.id for c in metabolites} == {1, 2}


def test_parse_metabolite_list_name_match_fallback(tmp_path):
    """When strict and loose InChI matches both miss, parser falls back to a
    case-insensitive name match (parser.py:142-144)."""
    (tmp_path / "good_chems.txt").write_text(
        "id\tname\tinchi\tsmiles\n"
        "1\tGlucose\tInChI=1S/realglucose\tnull\n"
    )
    (tmp_path / "metabolites.txt").write_text(
        "name\tinchi\tdescriptor\n"
        "glucose\tInChI=1S/differentstring\tfeedstock\n"
    )

    chems = parse_chemicals(tmp_path / "good_chems.txt")
    metabolites = parse_metabolite_list(tmp_path / "metabolites.txt", chems)

    assert len(metabolites) == 1
    assert next(iter(metabolites)).id == 1
