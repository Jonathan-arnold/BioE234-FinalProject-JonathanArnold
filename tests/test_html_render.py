"""Tests for synthesis_helper.mcp.html_render."""

from __future__ import annotations

import json
import re
from collections import Counter

import pytest

from synthesis_helper.mcp import html_render as hr
from synthesis_helper.models import HyperGraph
from synthesis_helper.pathways import enumerate_pathways
from synthesis_helper.traceback import traceback as build_cascade
from tests.fixtures.small_cell import SmallCell


_PAYLOAD_RE = re.compile(
    r'<script type="application/json" id="cy-data">(.*?)</script>',
    re.DOTALL,
)


def _extract_payload(html: str) -> dict:
    m = _PAYLOAD_RE.search(html)
    assert m is not None, "payload script tag not found"
    return json.loads(m.group(1))


def _content_nodes(payload: dict) -> list[dict]:
    return [n for n in payload["nodes"] if n["data"]["kind"] != "group"]


def _group_nodes(payload: dict) -> list[dict]:
    return [n for n in payload["nodes"] if n["data"]["kind"] == "group"]


# ---------------------------------------------------------------------------
# Pathway
# ---------------------------------------------------------------------------


@pytest.fixture
def target_and_pathways(small_cell: SmallCell, small_cell_hg: HyperGraph):
    T = small_cell.chems["T"]
    cascade = build_cascade(small_cell_hg, T, shell_cutoff=1)
    pathways = enumerate_pathways(cascade, small_cell_hg, max_pathways=10)
    assert pathways, "small_cell fixture should produce at least one pathway to T"
    return T, cascade, pathways


def test_pathway_html_starts_with_doctype(target_and_pathways, small_cell_hg):
    _, _, pathways = target_and_pathways
    html = hr.render_pathway_html(pathways[0], small_cell_hg)
    assert html.lstrip().lower().startswith("<!doctype html>")


def test_pathway_payload_node_counts(target_and_pathways, small_cell_hg):
    _, _, pathways = target_and_pathways
    pathway = pathways[0]
    html = hr.render_pathway_html(pathway, small_cell_hg)
    payload = _extract_payload(html)

    chem_nodes = [n for n in payload["nodes"] if n["data"]["kind"] == "chemical"]
    rxn_nodes = [n for n in payload["nodes"] if n["data"]["kind"] == "reaction"]

    # Every reaction in the pathway becomes one reaction node.
    assert len(rxn_nodes) == len(pathway.reactions)
    # Every chemical touched by a reaction (substrate or product) or the
    # target itself becomes one chemical node.
    expected_chems = set()
    for r in pathway.reactions:
        expected_chems.update(r.substrates)
        expected_chems.update(r.products)
    expected_chems.add(pathway.target)
    assert len(chem_nodes) == len(expected_chems)


def test_pathway_target_is_flagged(target_and_pathways, small_cell_hg):
    T, _, pathways = target_and_pathways
    html = hr.render_pathway_html(pathways[0], small_cell_hg)
    payload = _extract_payload(html)
    targets = [n for n in _content_nodes(payload) if n["data"].get("is_target")]
    assert len(targets) == 1
    assert targets[0]["data"]["chem_id"] == T.id


def test_pathway_themes_cover_every_content_node(target_and_pathways, small_cell_hg):
    _, _, pathways = target_and_pathways
    html = hr.render_pathway_html(pathways[0], small_cell_hg)
    payload = _extract_payload(html)
    for n in _content_nodes(payload):
        assert set(n["data"]["themes"].keys()) == set(hr.PATHWAY_THEMES), (
            f"node {n['data']['id']} missing themes"
        )


def test_pathway_group_parents_match_theme_assignments(target_and_pathways, small_cell_hg):
    _, _, pathways = target_and_pathways
    html = hr.render_pathway_html(pathways[0], small_cell_hg)
    payload = _extract_payload(html)

    # For each theme, the set of group node ids must equal the set of
    # group ids referenced by content nodes under that theme. No orphans,
    # no dangling references.
    groups_by_theme: dict[str, set[str]] = {}
    for g in _group_nodes(payload):
        groups_by_theme.setdefault(g["data"]["theme"], set()).add(g["data"]["id"])

    for theme in hr.PATHWAY_THEMES:
        referenced: set[str] = set()
        for n in _content_nodes(payload):
            gid = n["data"]["themes"].get(theme)
            if gid is not None:
                referenced.add(gid)
        assert groups_by_theme.get(theme, set()) == referenced, (
            f"theme {theme}: group parents do not match content assignments"
        )


def test_pathway_meta_defaults(target_and_pathways, small_cell_hg):
    _, _, pathways = target_and_pathways
    html = hr.render_pathway_html(pathways[0], small_cell_hg)
    payload = _extract_payload(html)
    meta = payload["meta"]
    assert meta["view_kind"] == "pathway"
    assert meta["default_theme"] == "shell"
    assert meta["themes"] == list(hr.PATHWAY_THEMES)


def test_pathway_reaction_flags_populated(small_cell, target_and_pathways, small_cell_hg):
    # Every reaction node must carry the four engineering flags. R6 in the
    # fixture (EC 1.14.13.1) has no ec_names entry → is_p450 ∧ is_heme ∧
    # is_orphan all true. R7 (EC 2.7.1.1 / hexokinase) has an entry → none
    # of the flags fire. Verify both cases appear in at least one rendered
    # pathway.
    _, _, pathways = target_and_pathways
    ec_names = {"2.7.1.1": "hexokinase", "1.1.1.1": "alcohol dehydrogenase"}
    required = {"is_p450", "is_heme", "is_orphan", "has_toxic_neighbor"}

    seen_p450_orphan = False
    seen_clean = False
    for pw in pathways:
        html = hr.render_pathway_html(pw, small_cell_hg, ec_to_name=ec_names)
        payload = _extract_payload(html)
        for n in payload["nodes"]:
            d = n["data"]
            if d["kind"] != "reaction":
                continue
            assert required <= set(d), f"rxn {d['rxn_id']} missing flags"
            if d["ecnum"] == "1.14.13.1":
                assert d["is_p450"] and d["is_heme"] and d["is_orphan"]
                seen_p450_orphan = True
            if d["ecnum"] == "2.7.1.1":
                assert not d["is_p450"] and not d["is_heme"] and not d["is_orphan"]
                seen_clean = True

    assert seen_p450_orphan, "expected at least one pathway to route through R6"
    assert seen_clean, "expected at least one pathway to route through R7"


def test_pathway_toxic_chem_flagging(small_cell, target_and_pathways, small_cell_hg):
    # Synthetic test: pass a toxic_chem_ids set that names an arbitrary
    # intermediate (M1 — acetaldehyde-flavored node in the fixture) and
    # verify it propagates to both the chemical node (is_toxic) and any
    # reaction whose substrate/product set touches it (has_toxic_neighbor).
    _, _, pathways = target_and_pathways
    M1 = small_cell.chems["M1"]
    toxic = {M1.id}

    saw_toxic_chem = False
    saw_touching_rxn = False
    for pw in pathways:
        html = hr.render_pathway_html(
            pw, small_cell_hg, toxic_chem_ids=toxic,
        )
        payload = _extract_payload(html)
        for n in payload["nodes"]:
            d = n["data"]
            if d["kind"] == "chemical" and d["chem_id"] == M1.id:
                assert d["is_toxic"] is True
                saw_toxic_chem = True
            if d["kind"] == "reaction":
                touches = M1.id in d["substrate_ids"] or M1.id in d["product_ids"]
                assert d["has_toxic_neighbor"] is touches

                if touches:
                    saw_touching_rxn = True

    # M1 is produced by R1 and consumed by R2 / R9, at least one must appear.
    if saw_toxic_chem:
        assert saw_touching_rxn, "M1 should be touched by at least one rxn"


# ---------------------------------------------------------------------------
# Cascade
# ---------------------------------------------------------------------------


def test_cascade_html_basic(target_and_pathways, small_cell_hg):
    _, cascade, _ = target_and_pathways
    html = hr.render_cascade_html(cascade, small_cell_hg)
    assert html.lstrip().lower().startswith("<!doctype html>")
    payload = _extract_payload(html)
    assert payload["meta"]["view_kind"] == "cascade"
    assert payload["meta"]["default_theme"] == "shell"
    assert payload["meta"]["themes"] == list(hr.CASCADE_THEMES)


def test_cascade_themes_cover_every_content_node(target_and_pathways, small_cell_hg):
    _, cascade, _ = target_and_pathways
    html = hr.render_cascade_html(cascade, small_cell_hg)
    payload = _extract_payload(html)
    for n in _content_nodes(payload):
        assert set(n["data"]["themes"].keys()) == set(hr.CASCADE_THEMES)


def test_cascade_group_parents_match_theme_assignments(target_and_pathways, small_cell_hg):
    _, cascade, _ = target_and_pathways
    html = hr.render_cascade_html(cascade, small_cell_hg)
    payload = _extract_payload(html)

    groups_by_theme: dict[str, set[str]] = {}
    for g in _group_nodes(payload):
        groups_by_theme.setdefault(g["data"]["theme"], set()).add(g["data"]["id"])

    for theme in hr.CASCADE_THEMES:
        referenced: set[str] = set()
        for n in _content_nodes(payload):
            gid = n["data"]["themes"].get(theme)
            if gid is not None:
                referenced.add(gid)
        assert groups_by_theme.get(theme, set()) == referenced, (
            f"cascade theme {theme}: group parents mismatch"
        )


def test_cascade_producer_depth_has_zero_for_target(target_and_pathways, small_cell_hg):
    T, cascade, _ = target_and_pathways
    html = hr.render_cascade_html(cascade, small_cell_hg)
    payload = _extract_payload(html)
    target_node = next(n for n in _content_nodes(payload) if n["data"].get("is_target"))
    assert target_node["data"]["themes"]["producer_depth"] == "pd_0"


def test_cascade_max_reactions_enforced(target_and_pathways, small_cell_hg):
    _, cascade, _ = target_and_pathways
    with pytest.raises(ValueError, match="exceeds max_reactions"):
        hr.render_cascade_html(cascade, small_cell_hg, max_reactions=0)


# ---------------------------------------------------------------------------
# Reaction node payload shape (for the double-click copy flow)
# ---------------------------------------------------------------------------


def test_reaction_nodes_carry_substrate_and_product_names(target_and_pathways, small_cell_hg):
    _, _, pathways = target_and_pathways
    pathway = pathways[0]
    html = hr.render_pathway_html(pathway, small_cell_hg)
    payload = _extract_payload(html)
    for n in payload["nodes"]:
        d = n["data"]
        if d["kind"] != "reaction":
            continue
        assert len(d["substrate_names"]) == len(d["substrate_ids"])
        assert len(d["product_names"]) == len(d["product_ids"])
        assert all(isinstance(s, str) and s for s in d["substrate_names"])
        assert all(isinstance(p, str) and p for p in d["product_names"])


# ---------------------------------------------------------------------------
# Theme-label helpers
# ---------------------------------------------------------------------------


def test_label_helpers():
    assert hr._label_for_group("shell", "shell_0") == "Shell 0 · Native"
    assert hr._label_for_group("shell", "shell_5") == "Shell 5"
    assert hr._label_for_group("ec_class", "ec_1") == "1.x Oxidoreductase"
    assert hr._label_for_group("ec_class", "ec_other") == "No EC / other"
    assert hr._label_for_group("role", "role_target") == "Target"
    assert hr._label_for_group("producer_depth", "pd_0") == "Target layer"
    assert hr._label_for_group("producer_depth", "pd_2") == "Producer depth 2"
