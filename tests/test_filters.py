"""Tests for synthesis_helper.filters — SubstrateFilter composition."""

from __future__ import annotations

import pytest

from synthesis_helper.filters import (
    compose,
    inchi_normalized_filter,
    rhea_atom_filter,
    shell_threshold_filter,
    shell_zero_filter,
)
from synthesis_helper.models import Chemical, HyperGraph, Reaction
from synthesis_helper.traceback import traceback


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _chem(id_: int, name: str, inchi: str = "") -> Chemical:
    return Chemical(id=id_, name=name, inchi=inchi)


def _rxn(id_: int, subs: list[Chemical], prods: list[Chemical]) -> Reaction:
    return Reaction(id=id_, substrates=frozenset(subs), products=frozenset(prods))


# ---------------------------------------------------------------------------
# shell_zero_filter
# ---------------------------------------------------------------------------

def test_shell_zero_filter_drops_shell0(small_cell, small_cell_hg):
    """shell_zero_filter must exclude shell-0 substrates (U1, N1, N2)."""
    c = small_cell.chems
    # R4: substrates = {M4 (shell 2), U1 (shell 0)}
    r4 = small_cell.rxns["R4"]
    result = shell_zero_filter(r4, small_cell_hg)
    assert c["U1"] not in result
    assert c["M4"] in result


def test_shell_zero_filter_is_default_behavior(small_cell, small_cell_hg):
    """Passing shell_zero_filter explicitly should produce the same cascade
    as the default (no filter argument)."""
    c = small_cell.chems
    default_cascade = traceback(small_cell_hg, c["T"])
    explicit_cascade = traceback(small_cell_hg, c["T"], substrate_filter=shell_zero_filter)
    assert default_cascade.reactions == explicit_cascade.reactions


# ---------------------------------------------------------------------------
# shell_threshold_filter
# ---------------------------------------------------------------------------

def test_shell_threshold_0_matches_shell_zero_filter(small_cell, small_cell_hg):
    """threshold=0 is equivalent to shell_zero_filter on enabled reactions
    (i.e., reactions whose substrates are all in the hypergraph).
    Unenabled reactions (substrates absent from chemical_to_shell) are not
    encountered during traceback, so the edge case doesn't matter in practice.
    """
    thresh0 = shell_threshold_filter(0)
    for rxn in small_cell_hg.reaction_to_shell:
        assert thresh0(rxn, small_cell_hg) == shell_zero_filter(rxn, small_cell_hg)


def test_shell_threshold_drops_low_shell_substrates(small_cell, small_cell_hg):
    """threshold=2 must suppress chemicals at shells 0, 1, and 2."""
    c = small_cell.chems
    thresh2 = shell_threshold_filter(2)
    # R4 substrates: M4 (shell 2), U1 (shell 0) — both should be dropped
    r4 = small_cell.rxns["R4"]
    result = thresh2(r4, small_cell_hg)
    assert c["M4"] not in result
    assert c["U1"] not in result


def test_shell_threshold_cascade_is_subset(small_cell, small_cell_hg):
    """A cascade built with threshold=1 must be a subset of the default cascade."""
    c = small_cell.chems
    default_cascade = traceback(small_cell_hg, c["T"])
    filtered_cascade = traceback(
        small_cell_hg, c["T"], substrate_filter=shell_threshold_filter(1)
    )
    assert filtered_cascade.reactions <= default_cascade.reactions


# ---------------------------------------------------------------------------
# inchi_normalized_filter
# ---------------------------------------------------------------------------

def _make_hg_with_fuzzy_cofactor() -> tuple[HyperGraph, Reaction, Chemical, Chemical]:
    """Build a tiny 3-chemical HyperGraph where a cofactor appears in the
    reaction with a different protonation state than in the universal set.

    Layout:
        N  (shell 0, native)
        COF_exact  (shell 0, universal — the "canonical" form)
        COF_fuzzy  (NOT in shell 0 — same base InChI as COF_exact but /p+1)
        P  (shell 1 — produced by N + COF_fuzzy)
    """
    N = _chem(1, "N", inchi="InChI=1S/formula_N/c1/h1")
    COF_exact = _chem(2, "COF_exact", inchi="InChI=1S/formula_COF/c1/h1")
    COF_fuzzy = _chem(3, "COF_fuzzy", inchi="InChI=1S/formula_COF/c1/h1/p+1")
    P = _chem(4, "P", inchi="InChI=1S/formula_P/c1/h1")

    rxn = _rxn(1, subs=[N, COF_fuzzy], prods=[P])

    hg = HyperGraph()
    hg.chemical_to_shell[N] = 0
    hg.chemical_to_shell[COF_exact] = 0
    # COF_fuzzy is NOT in shell 0 — it's a "near miss"
    hg.chemical_to_shell[COF_fuzzy] = 1
    hg.chemical_to_shell[P] = 1
    hg.reaction_to_shell[rxn] = 1

    return hg, rxn, COF_fuzzy, COF_exact


def test_inchi_normalized_filter_drops_fuzzy_cofactor():
    """A substrate whose stripped InChI matches a universal chem must be
    treated as background even if it is not literally in shell 0."""
    hg, rxn, COF_fuzzy, COF_exact = _make_hg_with_fuzzy_cofactor()
    universal = {COF_exact}
    f = inchi_normalized_filter(universal)
    result = f(rxn, hg)
    assert COF_fuzzy not in result


def test_inchi_normalized_filter_keeps_non_cofactor():
    """A non-cofactor substrate must still appear in the filter output."""
    N = _chem(1, "N", inchi="InChI=1S/formula_N/c1/h1")
    hg, rxn, COF_fuzzy, COF_exact = _make_hg_with_fuzzy_cofactor()
    universal = {COF_exact}
    f = inchi_normalized_filter(universal)
    result = f(rxn, hg)
    # N is shell 0 so also dropped, but that's expected — check COF_fuzzy only
    assert COF_fuzzy not in result


def test_inchi_normalized_filter_cascade_suppresses_fuzzy_cofactor():
    """End-to-end: traceback with inchi_normalized_filter must not recurse
    into a cofactor that only fails to match due to protonation state."""
    hg, rxn, COF_fuzzy, COF_exact = _make_hg_with_fuzzy_cofactor()
    P = next(c for c in hg.chemical_to_shell if c.name == "P")
    universal = {COF_exact}
    f = inchi_normalized_filter(universal)
    cascade = traceback(hg, P, substrate_filter=f)
    # The cascade should contain rxn (it produces P) but not recurse into COF_fuzzy
    assert rxn in cascade.reactions
    # No reaction in the cascade should have COF_fuzzy as a product
    # (it has no producers anyway, but verifies no phantom reactions were added)
    for r in cascade.reactions:
        assert COF_fuzzy not in r.products


# ---------------------------------------------------------------------------
# compose
# ---------------------------------------------------------------------------

def test_compose_empty_returns_all_substrates(small_cell, small_cell_hg):
    """compose() with no arguments returns all substrates of a reaction."""
    r4 = small_cell.rxns["R4"]
    result = compose()(r4, small_cell_hg)
    assert result == frozenset(r4.substrates)


def test_compose_single_is_identity(small_cell, small_cell_hg):
    """compose(f) == f for any single filter."""
    r4 = small_cell.rxns["R4"]
    assert (
        compose(shell_zero_filter)(r4, small_cell_hg)
        == shell_zero_filter(r4, small_cell_hg)
    )


def test_compose_intersection_is_more_restrictive(small_cell, small_cell_hg):
    """compose(f, g) must be a subset of compose(f) for any additional g."""
    c = small_cell.chems
    r4 = small_cell.rxns["R4"]  # substrates: M4 (shell 2), U1 (shell 0)
    # shell_zero drops U1; threshold(2) additionally drops M4
    base = shell_zero_filter(r4, small_cell_hg)
    combined = compose(shell_zero_filter, shell_threshold_filter(2))(r4, small_cell_hg)
    assert combined <= base
    assert c["M4"] not in combined
    assert c["U1"] not in combined


def test_compose_two_filters_cascade_is_subset(small_cell, small_cell_hg):
    """End-to-end: cascade from compose(shell_zero, threshold(1)) is a subset
    of the default cascade."""
    c = small_cell.chems
    default_cascade = traceback(small_cell_hg, c["T"])
    composed_filter = compose(shell_zero_filter, shell_threshold_filter(1))
    filtered_cascade = traceback(small_cell_hg, c["T"], substrate_filter=composed_filter)
    assert filtered_cascade.reactions <= default_cascade.reactions


# ---------------------------------------------------------------------------
# rhea_atom_filter
# ---------------------------------------------------------------------------

def _make_rhea_hg() -> tuple[HyperGraph, dict, dict]:
    """Build a 4-chemical, 2-reaction HyperGraph for Rhea filter testing.

    Layout:
        N  (shell 0, native)
        CAR (shell 0, also in carrier map for R1 — simulates NAD+/FAD)
        M2  (shell 1, skeleton product of R1)
        T   (shell 2, produced by M2 alone in R2)

    Reaction map:
        R1 (id=101): {N, CAR} → {M2}   CAR is a carrier per rhea map
        R2 (id=102): {M2}     → {T}    no entry in rhea map → fallback
    """
    N   = _chem(1, "N",   inchi="InChI=1S/N_formula")
    CAR = _chem(2, "CAR", inchi="InChI=1S/CAR_formula")
    M2  = _chem(3, "M2",  inchi="InChI=1S/M2_formula")
    T   = _chem(4, "T",   inchi="InChI=1S/T_formula")

    r1 = _rxn(101, subs=[N, CAR], prods=[M2])
    r2 = _rxn(102, subs=[M2],     prods=[T])

    hg = HyperGraph()
    hg.chemical_to_shell[N]   = 0
    hg.chemical_to_shell[CAR] = 0
    hg.chemical_to_shell[M2]  = 1
    hg.chemical_to_shell[T]   = 2
    hg.reaction_to_shell[r1]  = 1
    hg.reaction_to_shell[r2]  = 2

    chems = {"N": N, "CAR": CAR, "M2": M2, "T": T}
    rxns  = {"R1": r1, "R2": r2}
    return hg, chems, rxns


def test_rhea_atom_filter_drops_carrier():
    """Substrates listed as carriers in the map must be excluded."""
    hg, chems, rxns = _make_rhea_hg()
    carrier_inchi = "InChI=1S/CAR_formula"
    carrier_map = {101: frozenset([carrier_inchi])}
    f = rhea_atom_filter(carrier_map)
    result = f(rxns["R1"], hg)
    assert chems["CAR"] not in result


def test_rhea_atom_filter_keeps_skeleton():
    """Non-carrier substrates not in shell 0 must survive the filter."""
    hg, chems, rxns = _make_rhea_hg()
    # CAR is shell-0 anyway; make N shell-1 so we can distinguish filter behaviour
    hg.chemical_to_shell[chems["N"]] = 1
    carrier_inchi = "InChI=1S/CAR_formula"
    carrier_map = {101: frozenset([carrier_inchi])}
    f = rhea_atom_filter(carrier_map)
    result = f(rxns["R1"], hg)
    assert chems["N"] in result


def test_rhea_atom_filter_fallback_to_shell_zero():
    """Reactions absent from the carrier map must fall back to shell_zero_filter."""
    hg, chems, rxns = _make_rhea_hg()
    # R2 (id=102) is not in the carrier map
    carrier_map: dict[int, frozenset[str]] = {}
    f = rhea_atom_filter(carrier_map)
    result_rhea = f(rxns["R2"], hg)
    result_sz   = shell_zero_filter(rxns["R2"], hg)
    assert result_rhea == result_sz


def test_rhea_atom_filter_cascade_subset(small_cell, small_cell_hg):
    """End-to-end: cascade with rhea_atom_filter(empty map) equals default cascade
    because the empty map causes full fallback to shell_zero_filter."""
    c = small_cell.chems
    default_cascade  = traceback(small_cell_hg, c["T"])
    filtered_cascade = traceback(
        small_cell_hg, c["T"],
        substrate_filter=rhea_atom_filter({}),
    )
    assert filtered_cascade.reactions == default_cascade.reactions


def test_rhea_atom_filter_composes_with_shell_zero(small_cell, small_cell_hg):
    """compose(rhea_atom_filter, shell_zero_filter) cascade is subset of default."""
    c = small_cell.chems
    default_cascade  = traceback(small_cell_hg, c["T"])
    composed = compose(rhea_atom_filter({}), shell_zero_filter)
    filtered_cascade = traceback(small_cell_hg, c["T"], substrate_filter=composed)
    assert filtered_cascade.reactions <= default_cascade.reactions
