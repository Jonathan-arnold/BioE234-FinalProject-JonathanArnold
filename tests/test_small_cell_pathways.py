"""Scenario (d-pathways): enumerate_pathways must return every complete route from
native metabolites to the target.
"""

from __future__ import annotations

from synthesis_helper.parser import (
    parse_chemicals,
    parse_metabolite_list,
    parse_reactions,
)
from synthesis_helper.pathways import enumerate_pathways
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback
from tests.fixtures.small_cell import write_small_cell_tsvs


def _pathway_reaction_sets(pathways) -> set[frozenset]:
    return {frozenset(p.reactions) for p in pathways}


def test_pathways_finds_route_B(small_cell, small_cell_hg):
    """Route B ({R7}: N2 -> T) is a single-reaction pathway with an all-native
       substrate. It is recovered correctly even under the buggy partial-enumeration
       implementation — this test gives a stable green signal that the easy case works.
    """
    c = small_cell.chems
    r = small_cell.rxns

    cascade = traceback(small_cell_hg, c["T"], shell_cutoff=1)
    pathways = enumerate_pathways(cascade, small_cell_hg)

    reaction_sets = _pathway_reaction_sets(pathways)
    assert frozenset({r["R7"]}) in reaction_sets, (
        f"Route B ({{R7}}) missing from pathways; got {reaction_sets}"
    )


def test_pathways_enumerates_both_complete_routes(small_cell, small_cell_hg):
    """Two routes lead to T:
         Route B: {R7}                      (direct: N2 -> T)
         Route A: {R3, R4, R5, R6}          (N1+N2 -> M4 -> M5 -> T; N2 -> M3a -> T)
       A correct enumerate_pathways returns exactly these two pathways.
    """
    c = small_cell.chems
    r = small_cell.rxns

    cascade = traceback(small_cell_hg, c["T"], shell_cutoff=1)
    pathways = enumerate_pathways(cascade, small_cell_hg)

    reaction_sets = _pathway_reaction_sets(pathways)

    assert frozenset({r["R7"]}) in reaction_sets
    route_a = frozenset({r["R3"], r["R4"], r["R5"], r["R6"]})
    assert route_a in reaction_sets, (
        f"Route A {{R3, R4, R5, R6}} missing; got reaction sets {reaction_sets}"
    )
    assert len(pathways) == 2


def test_parser_roundtrip_tsv(small_cell, tmp_path):
    """Dump the fixture to TSV, parse back, re-run synthesize, and confirm the
       id->shell map is identical. Exercises parse_chemicals, parse_reactions,
       and parse_metabolite_list end-to-end against the small cell."""
    write_small_cell_tsvs(small_cell, tmp_path)

    chems_by_id = parse_chemicals(tmp_path / "enzymemap_chems.tsv")
    reactions = parse_reactions(tmp_path / "enzymemap_reactions.tsv", chems_by_id)
    natives = parse_metabolite_list(tmp_path / "minimal_metabolites.txt", chems_by_id)
    universals = parse_metabolite_list(
        tmp_path / "ubiquitous_metabolites.txt", chems_by_id
    )

    hg_tsv = synthesize(reactions, natives, universals, verbose=False)

    hg_py = synthesize(
        small_cell.all_reactions,
        small_cell.natives,
        small_cell.universals,
        verbose=False,
    )

    shells_py = {c.id: s for c, s in hg_py.chemical_to_shell.items()}
    shells_tsv = {c.id: s for c, s in hg_tsv.chemical_to_shell.items()}
    assert shells_tsv == shells_py
