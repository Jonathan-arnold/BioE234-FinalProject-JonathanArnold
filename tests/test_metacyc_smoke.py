"""End-to-end smoke test on the real MetaCyc corpus.

Pins shell numbers for well-known targets against the pre-computed
data/metacyc_L2_reachables.txt. If any of these drift, something in the
synthesize pipeline or shell-0 seeding has changed materially.

Skips gracefully when data/ files are absent (they are gitignored for
non-authors of the course corpus).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from synthesis_helper.parser import (
    parse_chemicals,
    parse_metabolite_list,
    parse_reactions,
)
from synthesis_helper.pathways import enumerate_pathways
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback

DATA_DIR = Path(__file__).parent.parent / "data"
REQUIRED_FILES = (
    "good_chems.txt",
    "good_reactions.txt",
    "minimal_metabolites.txt",
    "ubiquitous_metabolites.txt",
)


def _data_available() -> bool:
    return all((DATA_DIR / f).exists() for f in REQUIRED_FILES)


@pytest.fixture(scope="module")
def metacyc_hg():
    if not _data_available():
        pytest.skip(f"MetaCyc corpus files missing under {DATA_DIR}")

    chems = parse_chemicals(DATA_DIR / "good_chems.txt")
    rxns = parse_reactions(DATA_DIR / "good_reactions.txt", chems)
    natives = parse_metabolite_list(DATA_DIR / "minimal_metabolites.txt", chems)
    universals = parse_metabolite_list(DATA_DIR / "ubiquitous_metabolites.txt", chems)
    hg = synthesize(rxns, natives, universals, verbose=False)
    return chems, hg


def test_corpus_load_counts(metacyc_hg):
    """Parser and synthesize output snapshot for the committed corpus.

    These numbers are a regression canary — if they move, investigate before
    updating. A silent drift likely means the parser matching layer, the
    BFS termination, or the input files have changed unintentionally.
    """
    chems, hg = metacyc_hg
    assert len(chems) == 9361
    assert len(hg.chemical_to_shell) == 4310
    assert len(hg.reaction_to_shell) == 12138
    assert max(hg.chemical_to_shell.values()) == 25


def test_kanamycin_c_at_shell_8(metacyc_hg):
    """Lecture-cited example: kanamycin c reachable at shell 8 from E. coli
       minimals + ubiquitous cofactors. See CLAUDE.md §Example Query.
    """
    chems, hg = metacyc_hg
    target = chems[153337]
    assert target.name == "kanamycin c"
    assert hg.chemical_to_shell.get(target) == 8


def test_2prime_dehydrokanamycin_a_at_shell_10(metacyc_hg):
    """2'-dehydrokanamycin a sits one oxidation step past kanamycin c → shell 10.
       Pins the pre-computed metacyc_L2_reachables.txt snapshot.
    """
    chems, hg = metacyc_hg
    target = chems[148443]
    assert target.name == "2'-dehydrokanamycin a"
    assert hg.chemical_to_shell.get(target) == 10


def test_traceback_and_pathways_on_real_target(metacyc_hg):
    """Full pipeline on a real target must not crash and must return at least
       one pathway.

       Caps chosen conservatively: enumerate_pathways has a known exponential
       blowup in max_pathways when cascades branch (paths_for materializes
       max_pathways sub-paths per non-native substrate, cartesian-producting
       max_pathways**K combinations per merge point). Keep both caps small
       here so the smoke test stays sub-second; correctness at higher caps is
       covered by the small_cell fixture.
    """
    chems, hg = metacyc_hg
    target = chems[153337]  # kanamycin c

    cascade = traceback(hg, target, shell_cutoff=-1)
    assert cascade.target == target
    assert len(cascade.reactions) > 0

    pathways = enumerate_pathways(cascade, hg, max_pathways=10)
    assert len(pathways) > 0
    assert all(p.target == target for p in pathways)
    assert all(len(p.reactions) > 0 for p in pathways)
