"""Scenario (d-traceback): Cascade for a multi-producer target must include ALL routes."""

from __future__ import annotations

import pytest

from synthesis_helper.traceback import traceback


def test_cascade_for_T_contains_both_producers(small_cell, small_cell_hg):
    """Cascade(T) must include R7 (direct N2 -> T) AND R6 (M5+M3a -> T) plus
       R6's upstream chain: R4 (makes M5), R3 (makes M4 for R4), R5 (makes M3a).
    """
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(small_cell_hg, c["T"], shell_cutoff=1)

    assert cascade.target == c["T"]
    expected = {r["R7"], r["R6"], r["R4"], r["R3"], r["R5"]}
    assert cascade.reactions == expected, (
        f"Expected reactions {sorted(x.id for x in expected)}, "
        f"got {sorted(x.id for x in cascade.reactions)}"
    )


def test_cascade_stops_at_shell_zero(small_cell, small_cell_hg):
    """Cascade(T) must NOT include R1, R2, R2b — they produce chemicals unrelated to T."""
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(small_cell_hg, c["T"], shell_cutoff=1)
    for unrelated in ("R1", "R2", "R2b", "R8", "R9"):
        assert r[unrelated] not in cascade.reactions, (
            f"{unrelated} should not be in Cascade(T)"
        )


def test_cascade_for_intermediate_subset(small_cell, small_cell_hg):
    """Cascade(M5) should include only the chain producing M5 and its substrates:
       R4 (direct), R3 (makes M4). It must NOT include R5, R6, R7 — those are
       T-specific or M3a-specific.
    """
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(small_cell_hg, c["M5"])
    assert cascade.reactions == {r["R4"], r["R3"]}


def test_traceback_on_unreachable_raises(small_cell, small_cell_hg):
    """(e-traceback) Target not in the hypergraph raises ValueError."""
    c = small_cell.chems
    with pytest.raises(ValueError):
        traceback(small_cell_hg, c["M7"])


def test_traceback_shell_cutoff_minus_one_keeps_only_strictly_shorter_producers(
    small_cell, small_cell_hg
):
    """T sits at shell 1. Its producers are R7 (sub N2, shell 0) and R6
       (subs M5 shell 2, M3a shell 1). With shell_cutoff=-1 a substrate must
       satisfy shell <= 1 + (-1) = 0, so only R7 qualifies; R6 is pruned
       because M5 (and M3a) sit at or above T's own shell.
    """
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(small_cell_hg, c["T"], shell_cutoff=-1)

    assert cascade.reactions == {r["R7"]}, (
        f"shell_cutoff=-1 should retain only R7; got "
        f"{sorted(x.id for x in cascade.reactions)}"
    )


def test_traceback_shell_cutoff_one_recovers_full_cascade(small_cell, small_cell_hg):
    """shell_cutoff=1 lets a producer's substrates sit one shell above the
       chemical itself, which admits R6 (max sub-shell 2, T's shell 1). The
       resulting cascade matches the unrestricted ground truth.
    """
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(small_cell_hg, c["T"], shell_cutoff=1)
    assert cascade.reactions == {r["R7"], r["R6"], r["R4"], r["R3"], r["R5"]}


def test_traceback_default_is_unlimited(small_cell, small_cell_hg):
    """The default shell_cutoff is None (unlimited — no shell-based pruning)."""
    c = small_cell.chems
    default_cascade = traceback(small_cell_hg, c["T"])
    explicit = traceback(small_cell_hg, c["T"], shell_cutoff=None)
    assert default_cascade.reactions == explicit.reactions
    # And unlimited admits at least as many reactions as the tightest prune.
    pruned = traceback(small_cell_hg, c["T"], shell_cutoff=-1)
    assert pruned.reactions <= default_cascade.reactions


def test_traceback_shell_cutoff_below_minus_one_raises(small_cell, small_cell_hg):
    """shell_cutoff must be >= -1; -2 is invalid."""
    c = small_cell.chems
    with pytest.raises(ValueError):
        traceback(small_cell_hg, c["T"], shell_cutoff=-2)


def test_traceback_max_producers_keeps_lowest_shell_first(small_cell, small_cell_hg):
    """T has two producers: R7 (shell 1) and R6 (shell 3). With shell_cutoff
       wide enough to admit both, max_producers_per_chemical=1 should keep
       only the lower-shell producer (R7).
    """
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(
        small_cell_hg,
        c["T"],
        shell_cutoff=2,
        max_producers_per_chemical=1,
    )
    assert cascade.reactions == {r["R7"]}


def test_traceback_max_producers_filter_then_cap(small_cell, small_cell_hg):
    """The cap applies AFTER the shell filter — shell_cutoff=-1 already
       drops R6, so the cap has nothing to slice.
    """
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(
        small_cell_hg,
        c["T"],
        shell_cutoff=-1,
        max_producers_per_chemical=10,
    )
    assert cascade.reactions == {r["R7"]}


def test_traceback_max_producers_none_disables_cap(small_cell, small_cell_hg):
    """Explicit None for the cap matches the cap-absent behavior."""
    c = small_cell.chems
    capped = traceback(
        small_cell_hg, c["T"], shell_cutoff=2, max_producers_per_chemical=None
    )
    uncapped = traceback(small_cell_hg, c["T"], shell_cutoff=2)
    assert capped.reactions == uncapped.reactions


def test_traceback_max_producers_zero_raises(small_cell, small_cell_hg):
    """max_producers_per_chemical must be None or >= 1."""
    c = small_cell.chems
    with pytest.raises(ValueError):
        traceback(small_cell_hg, c["T"], max_producers_per_chemical=0)
