"""Tests for the synthesis_helper.mcp.server tools + resources.

We drive the small_cell fixture through the real lazy-loading state
(writing TSVs into a tmp data dir) and then call the tool functions
directly. The FastMCP @tool decorator keeps them as plain callables,
so exercising them this way covers all business logic without async
plumbing.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from synthesis_helper.mcp import state
from synthesis_helper.mcp.server import (
    chemical_resource,
    compare_pathways,
    describe_pathway,
    enumerate_pathways_for,
    find_by_structure,
    find_chemical,
    get_cascade,
    get_shell,
    list_reachables,
    pathway_to_composition,
    reachables_resource,
    reaction_resource,
    resynthesize_with_fed,
    stats_resource,
)
from tests.fixtures.small_cell import build_small_cell, write_small_cell_tsvs


@pytest.fixture(scope="module")
def small_cell_data_dir(tmp_path_factory) -> Path:
    """Write small_cell TSVs to a fresh tmp dir once per module."""
    sc = build_small_cell()
    d = tmp_path_factory.mktemp("mcp_small_cell_data")
    write_small_cell_tsvs(sc, d)
    return d


@pytest.fixture(autouse=True)
def bootstrap_from_small_cell(small_cell_data_dir: Path):
    """Reset + point state at the small_cell tmp dir before each test."""
    state.reset()
    state.configure(data_dir=small_cell_data_dir)
    yield
    state.reset()


# --- tools -------------------------------------------------------------------


def test_find_chemical_by_exact_name():
    hits = find_chemical("T")
    assert any(h["name"] == "T" and h["id"] == 300 for h in hits)


def test_find_chemical_case_insensitive_substring():
    hits = find_chemical("m_dee")
    names = {h["name"] for h in hits}
    assert "M_deep" in names


def test_get_shell_native():
    assert get_shell("N1")["shell"] == 0


def test_get_shell_target_shortcut():
    """T is reachable via R7 (direct N2 → T) at shell 1 even though R6 is shell 3."""
    assert get_shell("T")["shell"] == 1


def test_get_shell_unreachable_returns_null_shell():
    """M7 is never reachable (only producer needs X_missing)."""
    assert get_shell("M7")["shell"] is None


def test_list_reachables_shell_filter():
    out = list_reachables(max_shell=1)
    assert all(c["shell"] <= 1 for c in out)
    names = {c["name"] for c in out}
    assert {"N1", "N2", "U1"}.issubset(names)
    assert "M_deep" not in names  # shell 3


def test_list_reachables_name_filter():
    out = list_reachables(name_contains="M3")
    names = {c["name"] for c in out}
    assert names <= {"M3a", "M3b"}
    assert names  # at least one


def test_get_cascade_for_target():
    # shell_cutoff=1 admits R6's longer route (subs span shells 1-2 above T's
    # own shell-1); the default (-1) prunes it and would leave only R7.
    cascade = get_cascade("T", shell_cutoff=1)
    assert cascade["target"]["name"] == "T"
    # Every leaf must be shell 0
    assert all(leaf["shell"] == 0 for leaf in cascade["leaf_metabolites"])
    rxn_ids = {r["id"] for r in cascade["reactions"]}
    # Both direct (R7 id=8) and indirect (R6 id=7) routes appear in cascade
    assert 7 in rxn_ids and 8 in rxn_ids


def test_get_cascade_max_producers_cap_caps_branching():
    # T has two producers (R6 id=7 shell-3, R7 id=8 shell-1). With the shell
    # filter loose enough to admit both, max_producers_per_chemical=1 should
    # keep only the lower-shell producer (R7).
    cascade = get_cascade("T", shell_cutoff=2, max_producers_per_chemical=1)
    rxn_ids = {r["id"] for r in cascade["reactions"]}
    assert rxn_ids == {8}


def test_enumerate_pathways_for_target():
    pathways = enumerate_pathways_for("T", max_pathways=5)
    assert pathways, "expected at least one pathway to T"
    # every pathway's final reaction must produce T
    for p in pathways:
        last = p["reactions"][-1]
        assert 300 in last["product_ids"]


def test_describe_pathway_markdown_contains_arrow():
    md = describe_pathway("T", pathway_index=0)
    assert "→" in md
    assert "Pathway 0" in md


def test_ec_names_loaded_from_data_dir(tmp_path):
    """state._bootstrap picks up optional ec_names.tsv."""
    from synthesis_helper.mcp.server import _resolve_or_raise  # noqa: F401
    sc = build_small_cell()
    d = tmp_path / "cell_with_ec"
    d.mkdir()
    write_small_cell_tsvs(sc, d)
    (d / "ec_names.tsv").write_text("ecnum\tname\n1.1.1.1\tAlcohol dehydrogenase\n")

    state.reset()
    state.configure(data_dir=d)
    assert state.get_ec_names() == {"1.1.1.1": "Alcohol dehydrogenase"}
    state.reset()


def test_resynthesize_with_fed_enables_new_chemicals():
    baseline_hg = state.get_hypergraph()
    assert not any(c.name == "M6" for c in baseline_hg.chemical_to_shell)

    out = resynthesize_with_fed(fed_chemical_refs=["F1"])
    assert out["delta_chemicals"] >= 2  # F1 itself + M6
    new_names = {c["name"] for c in out["new_reachables"]}
    assert {"F1", "M6"}.issubset(new_names)


def test_resynthesize_unresolved_refs_are_reported():
    out = resynthesize_with_fed(fed_chemical_refs=["F1", "NO_SUCH_THING"])
    assert "NO_SUCH_THING" in out["unresolved"]


# --- resources ---------------------------------------------------------------


def test_stats_resource_shape():
    s = stats_resource()
    assert s["total_chemicals"] == 15  # small_cell has 15 chems (incl X_missing)
    assert s["reachable_chemicals"] >= 8  # most are reachable
    assert s["shells"] >= 3


def test_reachables_resource_tsv_header():
    tsv = reachables_resource()
    first = tsv.splitlines()[0]
    assert first == "id\tname\tinchi\tshell"
    # N1 (native, shell 0) should appear
    assert any(line.startswith("101\tN1\t") for line in tsv.splitlines())


def test_chemical_resource_produced_and_consumed():
    res = chemical_resource("300")  # T
    assert res["name"] == "T"
    # T is produced by R6 (id=7) and R7 (id=8)
    assert {7, 8}.issubset(set(res["produced_by_reaction_ids"]))


def test_reaction_resource_expands_chemicals():
    res = reaction_resource("8")  # R7: N2 → T
    assert 102 in res["substrate_ids"]
    assert 300 in res["product_ids"]
    sub_names = {s["name"] for s in res["substrates"]}
    prod_names = {p["name"] for p in res["products"]}
    assert sub_names == {"N2"}
    assert prod_names == {"T"}


def test_unknown_chemical_resource_raises():
    with pytest.raises(ValueError):
        chemical_resource("99999999")


def test_unknown_reaction_resource_raises():
    with pytest.raises(ValueError):
        reaction_resource("99999999")


# ---------------------------------------------------------------------------
# find_by_structure  (SMILES + InChI → Morgan FP + Tanimoto)
# ---------------------------------------------------------------------------


def test_find_by_structure_exact_match():
    hits = find_by_structure("CCO", limit=5, similarity_threshold=0.0)
    assert hits, "expected at least one structural hit"
    top = hits[0]
    assert top["name"] == "N1"
    assert top["tanimoto"] >= 0.99


def test_find_by_structure_includes_shell():
    hits = find_by_structure("CCO", limit=3, similarity_threshold=0.0)
    assert "shell" in hits[0]
    assert hits[0]["shell"] == 0  # N1 is a native


def test_find_by_structure_ranks_similar_structures():
    # CCO=ethanol (N1), CC=O=acetaldehyde (M1), CC(=O)O=acetate (M2).
    # All three carry real SMILES in the fixture. Querying with acetaldehyde
    # should surface the other two alcohols with reasonable scores.
    hits = find_by_structure("CC=O", limit=5, similarity_threshold=0.0)
    names = [h["name"] for h in hits]
    assert "M1" in names
    # Expect the related structures to also appear
    assert set(names) & {"N1", "M2"}, f"expected overlap with N1/M2 in {names}"


def test_find_by_structure_threshold_filters():
    # Benzene is nothing like the three 2-carbon molecules in the fixture
    # — at a high threshold every hit should be dropped.
    hits = find_by_structure("c1ccccc1", limit=10, similarity_threshold=0.9)
    assert hits == []


def test_find_by_structure_invalid_smiles_raises():
    with pytest.raises(ValueError, match="Invalid SMILES"):
        find_by_structure("not_a_smiles_$$$", limit=5)


def test_find_by_structure_exposes_fingerprint_coverage_in_stats():
    # Trigger index build.
    find_by_structure("CCO", limit=1, similarity_threshold=0.0)
    s = stats_resource()
    fp = s.get("fingerprint_coverage")
    assert fp is not None
    assert fp["indexed"] >= 3   # at least N1, M1, M2
    assert fp["fraction"] >= 0.0
    assert "build_ms" in fp


# ---------------------------------------------------------------------------
# pathway_to_composition  (annotated enzyme list)
# ---------------------------------------------------------------------------


def test_pathway_to_composition_shape():
    out = pathway_to_composition("T", pathway_index=0)
    assert out["target"]["name"] == "T"
    assert out["step_count"] == len(out["enzymes"])
    for e in out["enzymes"]:
        assert {"step", "ecnum", "enzyme_name", "is_orphan", "is_p450",
                "is_heme", "ec_class", "substrate_ids", "product_ids"} <= set(e)


def test_pathway_to_composition_flags_p450_and_orphan_consistently():
    # There are multiple pathways to T (via R7, and via R3→R4→R5+R6). One
    # path hits R7 (EC 2.7.1.1 — hexokinase, not orphan) and another hits
    # R6 (EC 1.14.13.1 — P450, deliberately absent from ec_names → orphan).
    # Walk every pathway and confirm the per-step invariants hold.
    for i in range(3):
        try:
            out = pathway_to_composition("T", pathway_index=i, shell_cutoff=1)
        except ValueError:
            break
        for e in out["enzymes"]:
            ec = e["ecnum"]
            if ec.startswith("1.14."):
                assert e["is_p450"] is True
                assert e["is_heme"] is True  # heme is strict superset of P450
                assert e["is_orphan"] is True  # no entry in fixture's ec_names
            else:
                assert e["is_p450"] is False
            # Heme ⊇ P450, so is_heme ⇒ is_p450 OR (peroxidase / dioxygenase)
            if e["is_heme"]:
                assert (
                    ec.startswith("1.14.")
                    or ec.startswith("1.11.1.")
                    or ec.startswith("1.13.11.")
                )
            if ec == "2.7.1.1":
                assert e["is_orphan"] is False
                assert e["is_heme"] is False
                assert e["enzyme_name"] == "hexokinase"
                assert e["ec_class"] == "Transferase"


def test_pathway_to_composition_out_of_range_raises():
    with pytest.raises(ValueError, match="out of range"):
        pathway_to_composition("T", pathway_index=9999)


def test_pathway_to_composition_unreachable_raises():
    with pytest.raises(ValueError, match="not reachable"):
        pathway_to_composition("M7")


# ---------------------------------------------------------------------------
# compare_pathways  (multi-dim scorecard)
# ---------------------------------------------------------------------------


def test_compare_pathways_basic_shape():
    out = compare_pathways("T", n=5)
    assert out["target"]["name"] == "T"
    assert isinstance(out["comparison"], list)
    # Every row must carry the expected cofactor keys, even when 0.
    expected_keys = {"NADPH", "NADH", "NADP+", "NAD+", "ATP", "CoA", "SAM"}
    assert set(out["cofactor_keys"]) == expected_keys
    # Every row must carry all flag fields (count or hint) even when 0.
    flag_fields = {
        "n_p450", "n_heme", "n_orphan_steps", "n_toxic_intermediate",
        "toxic_intermediates", "nadh_hint", "nadph_hint", "atp_hint",
    }
    for row in out["comparison"]:
        assert set(row["cofactor_uses"]) == expected_keys
        assert flag_fields <= set(row)


def test_compare_pathways_sort_order():
    # Rows are sorted by (step_count, n_heme, n_toxic_intermediate,
    # n_orphan_steps) — all ascending.
    out = compare_pathways("T", n=5)
    rows = out["comparison"]
    keys = [
        (r["step_count"], r["n_heme"], r["n_toxic_intermediate"], r["n_orphan_steps"])
        for r in rows
    ]
    assert keys == sorted(keys)


def test_compare_pathways_flags_p450_pathway():
    # At least one of the T pathways includes R6 (EC 1.14.13.1). Find it
    # and assert the P450 + heme + orphan flags propagate up to the scorecard.
    out = compare_pathways("T", n=5, shell_cutoff=1)
    has_p450_row = any(r["n_p450"] >= 1 for r in out["comparison"])
    assert has_p450_row, "expected at least one T pathway to hit R6 (EC 1.14.)"
    # n_heme must be >= n_p450 per-row (heme is a strict superset).
    for r in out["comparison"]:
        assert r["n_heme"] >= r["n_p450"]
    has_orphan_row = any(r["n_orphan_steps"] >= 1 for r in out["comparison"])
    assert has_orphan_row


def test_compare_pathways_pathway_index_preserved():
    # After sort, each row must still carry its ORIGINAL enumeration index
    # so the LLM can drill down via pathway_to_composition(target, pidx).
    out = compare_pathways("T", n=5)
    indices = sorted(r["pathway_index"] for r in out["comparison"])
    assert indices == list(range(len(indices)))
    # Every original index is unique
    assert len(indices) == len(set(indices))


def test_compare_pathways_unreachable_raises():
    with pytest.raises(ValueError, match="not reachable"):
        compare_pathways("M7")
