"""Interactive HTML pathway / cascade viewer (Cytoscape.js + Flat Design).

Emits a self-contained HTML file with an inline Cytoscape.js bundle. The
rendered page supports:

- 5 layout algorithms (preset / breadthfirst / concentric / cose / grid)
- Per-view group-by themes (pathway: shell/step/ec_class/role;
  cascade: shell/ec_class/role/producer_depth), swapped client-side
  with zero server round-trip
- Draggable nodes
- Double-click a node to copy its structured info to the clipboard
"""

from __future__ import annotations

import json
from collections import deque
from pathlib import Path
from typing import Iterable, Mapping

from synthesis_helper.models import (
    Cascade,
    Chemical,
    HyperGraph,
    Pathway,
    Reaction,
)

_ASSETS = Path(__file__).parent / "assets"
_CYTOSCAPE_JS = _ASSETS / "cytoscape.min.js"

PATHWAY_THEMES: tuple[str, ...] = ("shell", "ec_class", "role")
CASCADE_THEMES: tuple[str, ...] = ("shell", "ec_class", "role", "producer_depth")


# ---------------------------------------------------------------------------
# Theme computation
# ---------------------------------------------------------------------------

_EC_LABELS: dict[str, str] = {
    "ec_1": "1.x Oxidoreductase",
    "ec_2": "2.x Transferase",
    "ec_3": "3.x Hydrolase",
    "ec_4": "4.x Lyase",
    "ec_5": "5.x Isomerase",
    "ec_6": "6.x Ligase",
    "ec_other": "No EC / other",
}

_ROLE_LABELS: dict[str, str] = {
    "role_native": "Native",
    "role_intermediate": "Intermediate",
    "role_target": "Target",
}


# Engineering-relevant flag detection. Mirrors the definitions used by
# compare_pathways / annotate_pathway so the visualizer's overlay agrees with
# the numeric scorecards.
def _is_p450_ec(ecnum: str) -> bool:
    return ecnum.startswith("1.14.")


def _is_heme_ec(ecnum: str) -> bool:
    return (
        ecnum.startswith("1.14.")
        or ecnum.startswith("1.11.1.")
        or ecnum.startswith("1.13.11.")
    )


def _is_orphan_ec(ecnum: str, ec_to_name: Mapping[str, str] | None) -> bool:
    if not ecnum:
        return True
    if ec_to_name is None:
        return False
    return ecnum not in ec_to_name


def _ec_class_group(ecnum: str) -> str:
    if not ecnum:
        return "ec_other"
    first = ecnum.split(".", 1)[0]
    if first in {"1", "2", "3", "4", "5", "6"}:
        return f"ec_{first}"
    return "ec_other"


def _chem_role(chem: Chemical, hg: HyperGraph, target: Chemical) -> str:
    if chem == target:
        return "role_target"
    if hg.chemical_to_shell.get(chem) == 0:
        return "role_native"
    return "role_intermediate"


def _label_for_group(theme: str, gid: str) -> str:
    if theme == "shell":
        n = gid[len("shell_") :]
        if n.isdigit():
            return "Shell 0 · Native" if n == "0" else f"Shell {n}"
        return "Shell ?"
    if theme == "ec_class":
        return _EC_LABELS.get(gid, gid)
    if theme == "role":
        return _ROLE_LABELS.get(gid, gid)
    if theme == "producer_depth":
        n = gid[len("pd_") :]
        if n.isdigit():
            return "Target layer" if n == "0" else f"Producer depth {n}"
        return "Producer depth ?"
    return gid


def _compute_themes_for_pathway(
    pathway: Pathway, hg: HyperGraph
) -> tuple[
    dict[Chemical, dict[str, str | None]],
    dict[Reaction, dict[str, str | None]],
    list[tuple[str, str, str]],
]:
    target = pathway.target
    chem_themes: dict[Chemical, dict[str, str | None]] = {}
    rxn_themes: dict[Reaction, dict[str, str | None]] = {}

    for chem in pathway.metabolites:
        chem_themes[chem] = {
            "shell": f"shell_{hg.chemical_to_shell.get(chem, 0)}",
            "ec_class": None,
            "role": _chem_role(chem, hg, target),
        }

    for rxn in pathway.reactions:
        rxn_themes[rxn] = {
            "shell": f"shell_{hg.reaction_to_shell.get(rxn, 0)}",
            "ec_class": _ec_class_group(rxn.ecnum),
            "role": None,
        }

    group_specs = _collect_group_specs(chem_themes, rxn_themes, PATHWAY_THEMES)
    return chem_themes, rxn_themes, group_specs


def _compute_themes_for_cascade(
    cascade: Cascade, hg: HyperGraph
) -> tuple[
    dict[Chemical, dict[str, str | None]],
    dict[Reaction, dict[str, str | None]],
    list[tuple[str, str, str]],
]:
    target = cascade.target
    chem_themes: dict[Chemical, dict[str, str | None]] = {}
    rxn_themes: dict[Reaction, dict[str, str | None]] = {}

    # Producer-depth BFS from target backward.
    producers_of: dict[Chemical, list[Reaction]] = {}
    for r in cascade.reactions:
        for p in r.products:
            producers_of.setdefault(p, []).append(r)

    chem_depth: dict[Chemical, int] = {target: 0}
    rxn_depth: dict[Reaction, int] = {}
    queue: deque[Chemical] = deque([target])
    while queue:
        chem = queue.popleft()
        d = chem_depth[chem]
        for r in producers_of.get(chem, []):
            if r in rxn_depth:
                continue
            rxn_depth[r] = d
            for s in r.substrates:
                if s not in chem_depth:
                    chem_depth[s] = d + 1
                    queue.append(s)

    all_chems: set[Chemical] = {target}
    for r in cascade.reactions:
        all_chems.update(r.substrates)
        all_chems.update(r.products)

    for chem in all_chems:
        d = chem_depth.get(chem)
        chem_themes[chem] = {
            "shell": f"shell_{hg.chemical_to_shell.get(chem, 0)}",
            "ec_class": None,
            "role": _chem_role(chem, hg, target),
            "producer_depth": f"pd_{d}" if d is not None else None,
        }

    for rxn in cascade.reactions:
        d = rxn_depth.get(rxn)
        rxn_themes[rxn] = {
            "shell": f"shell_{hg.reaction_to_shell.get(rxn, 0)}",
            "ec_class": _ec_class_group(rxn.ecnum),
            "role": None,
            "producer_depth": f"pd_{d}" if d is not None else None,
        }

    group_specs = _collect_group_specs(chem_themes, rxn_themes, CASCADE_THEMES)
    return chem_themes, rxn_themes, group_specs


def _collect_group_specs(
    chem_themes: Mapping[Chemical, dict[str, str | None]],
    rxn_themes: Mapping[Reaction, dict[str, str | None]],
    themes: Iterable[str],
) -> list[tuple[str, str, str]]:
    """Returns unique (group_id, theme, label) triples, sorted for stable output."""
    seen: set[tuple[str, str]] = set()
    for store in (chem_themes, rxn_themes):
        for td in store.values():
            for theme in themes:
                gid = td.get(theme)
                if gid is None:
                    continue
                seen.add((theme, gid))
    return sorted(
        ((gid, theme, _label_for_group(theme, gid)) for (theme, gid) in seen),
        key=lambda t: (t[1], t[0]),
    )


# ---------------------------------------------------------------------------
# Payload builder
# ---------------------------------------------------------------------------

def _real_name(name: str | None) -> str:
    """Return `name` if it's a real label, else empty string.
    Source data (MetaCyc) uses the literal sentinel "undefined" for anonymous
    chemicals — treat that as missing."""
    if name and name != "undefined":
        return name
    return ""


def _real_smiles(smiles: str | None) -> str:
    """Empty-out SMILES sentinels ("null", "")."""
    if smiles and smiles != "null":
        return smiles
    return ""


def _chem_label_short(chem: Chemical) -> str:
    real = _real_name(chem.name)
    return real if real else f"#{chem.id}"


def _rxn_label_short(rxn: Reaction, ec_to_name: Mapping[str, str] | None) -> str:
    if rxn.ecnum and ec_to_name and rxn.ecnum in ec_to_name:
        return ec_to_name[rxn.ecnum]
    if rxn.ecnum:
        return f"EC {rxn.ecnum}"
    return f"rxn #{rxn.id}"


def _build_cytoscape_payload(
    reactions: Iterable[Reaction],
    target: Chemical,
    hg: HyperGraph,
    ec_to_name: Mapping[str, str] | None,
    chem_themes: Mapping[Chemical, dict[str, str | None]],
    rxn_themes: Mapping[Reaction, dict[str, str | None]],
    group_specs: list[tuple[str, str, str]],
    view_kind: str,
    themes_in_view: tuple[str, ...],
    chem_id_to_name: Mapping[int, str],
    currency_chemicals: set[Chemical] | None = None,
    toxic_chem_ids: set[int] | None = None,
) -> dict:
    currency_set = currency_chemicals or set()
    toxic_ids = toxic_chem_ids or set()
    reactions_list = list(reactions)

    chems_used: dict[Chemical, int | None] = {}
    for r in reactions_list:
        for c in r.substrates:
            chems_used.setdefault(c, hg.chemical_to_shell.get(c))
        for c in r.products:
            chems_used.setdefault(c, hg.chemical_to_shell.get(c))
    chems_used.setdefault(target, hg.chemical_to_shell.get(target))

    nodes: list[dict] = []

    for gid, theme, label in group_specs:
        nodes.append(
            {"data": {"id": gid, "kind": "group", "theme": theme, "label": label}}
        )

    for chem in sorted(chems_used, key=lambda c: c.id):
        shell = chems_used[chem]
        is_target = chem == target
        themes_dict = chem_themes.get(
            chem, {t: None for t in themes_in_view}
        )
        nodes.append(
            {
                "data": {
                    "id": f"c{chem.id}",
                    "kind": "chemical",
                    "label": _chem_label_short(chem),
                    "shell": shell,
                    "is_target": is_target,
                    "is_currency": chem in currency_set,
                    "is_toxic": chem.id in toxic_ids,
                    "name": _real_name(chem.name),
                    "chem_id": chem.id,
                    "inchi": chem.inchi or "",
                    "smiles": _real_smiles(chem.smiles),
                    "themes": dict(themes_dict),
                }
            }
        )

    for rxn in sorted(reactions_list, key=lambda r: r.id):
        enzyme_name = ""
        if rxn.ecnum and ec_to_name:
            enzyme_name = ec_to_name.get(rxn.ecnum, "")
        themes_dict = rxn_themes.get(rxn, {t: None for t in themes_in_view})
        subs_sorted = sorted(rxn.substrates, key=lambda x: x.id)
        prods_sorted = sorted(rxn.products, key=lambda x: x.id)
        touching_ids = {s.id for s in subs_sorted} | {p.id for p in prods_sorted}
        ec = rxn.ecnum or ""
        nodes.append(
            {
                "data": {
                    "id": f"r{rxn.id}",
                    "kind": "reaction",
                    "label": _rxn_label_short(rxn, ec_to_name),
                    "shell": hg.reaction_to_shell.get(rxn),
                    "rxn_id": rxn.id,
                    "ecnum": ec,
                    "enzyme_name": enzyme_name,
                    "is_p450": _is_p450_ec(ec),
                    "is_heme": _is_heme_ec(ec),
                    "is_orphan": _is_orphan_ec(ec, ec_to_name),
                    "has_toxic_neighbor": bool(touching_ids & toxic_ids),
                    "substrate_ids": [s.id for s in subs_sorted],
                    "product_ids": [p.id for p in prods_sorted],
                    "substrate_names": [
                        chem_id_to_name.get(s.id, f"#{s.id}") for s in subs_sorted
                    ],
                    "product_names": [
                        chem_id_to_name.get(p.id, f"#{p.id}") for p in prods_sorted
                    ],
                    "themes": dict(themes_dict),
                }
            }
        )

    edges: list[dict] = []
    for rxn in reactions_list:
        for s in sorted(rxn.substrates, key=lambda x: x.id):
            edges.append(
                {
                    "data": {
                        "id": f"e_c{s.id}_r{rxn.id}",
                        "source": f"c{s.id}",
                        "target": f"r{rxn.id}",
                        "kind": "sub",
                    }
                }
            )
        for p in sorted(rxn.products, key=lambda x: x.id):
            edges.append(
                {
                    "data": {
                        "id": f"e_r{rxn.id}_c{p.id}",
                        "source": f"r{rxn.id}",
                        "target": f"c{p.id}",
                        "kind": "prod",
                    }
                }
            )

    return {
        "nodes": nodes,
        "edges": edges,
        "meta": {
            "view_kind": view_kind,
            "target_id": f"c{target.id}",
            "target_name": _real_name(target.name),
            "target_chem_id": target.id,
            "target_shell": hg.chemical_to_shell.get(target),
            "themes": list(themes_in_view),
            "default_theme": "shell",
            "counts": {
                "chemicals": len(chems_used),
                "reactions": len(reactions_list),
                "groups": len(group_specs),
            },
        },
    }


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
  :root {
    --bg: #FAFAFA;
    --toolbar-bg: #FFFFFF;
    --toolbar-divider: #E0E0E0;
    --button-active: #E3F2FD;
    --text: #212121;
    --subtle: #616161;
    --ok: #66BB6A;
  }
  * { box-sizing: border-box; }
  html, body { margin: 0; padding: 0; height: 100%; background: var(--bg);
               font-family: -apple-system, "Helvetica Neue", Arial, sans-serif;
               color: var(--text); font-size: 14px; }
  #toolbar {
    position: fixed; top: 0; left: 0; right: 0; height: 160px;
    background: var(--toolbar-bg); border-bottom: 1px solid var(--toolbar-divider);
    display: flex; align-items: center; padding: 0 16px; gap: 16px; z-index: 10;
  }
  #title { flex: 0 0 auto; min-width: 240px; }
  #title h1 { margin: 0; font-size: 15px; font-weight: 600; letter-spacing: 0.2px; }
  #title .sub { font-size: 12px; color: var(--subtle); margin-top: 4px; }
  #pickers { flex: 1 1 auto; display: flex; flex-direction: column; gap: 6px; min-width: 0; }
  .picker-row { display: flex; align-items: center; gap: 8px; }
  .picker-row .lbl { font-size: 12px; color: var(--subtle); width: 72px; flex: 0 0 72px; }
  .picker-row .buttons { display: flex; gap: 6px; flex-wrap: wrap; }
  button.pick {
    background: transparent; border: 1px solid transparent;
    padding: 4px 10px; font-size: 13px; font-family: inherit;
    color: var(--text); cursor: pointer; border-radius: 4px;
    transition: background-color 80ms ease;
  }
  button.pick:hover { background: #F5F5F5; }
  button.pick.active { background: var(--button-active); border-color: #BBDEFB; }
  #actions { flex: 0 0 auto; display: flex; gap: 8px; align-items: center; }
  #actions input {
    font-family: inherit; font-size: 13px; padding: 4px 8px;
    border: 1px solid var(--toolbar-divider); border-radius: 4px; outline: none;
    width: 180px;
  }
  #actions input:focus { border-color: #42A5F5; }
  #cy { position: fixed; top: 160px; left: 0; right: 0; bottom: 0; background: var(--bg); }
  button.pick .count {
    margin-left: 6px; font-size: 11px; color: var(--subtle);
    font-variant-numeric: tabular-nums;
  }
  button.pick.active .count { color: var(--text); }
  button.pick:disabled {
    color: #BDBDBD; cursor: not-allowed; background: transparent;
  }
  #toast {
    position: fixed; bottom: 16px; right: 16px;
    background: #FFFFFF; border: 1px solid var(--ok);
    padding: 10px 20px; font-size: 13px; color: var(--text);
    border-radius: 4px; opacity: 0; pointer-events: none;
    transition: opacity 150ms ease;
  }
  #toast.show { opacity: 1; }
  #tooltip {
    position: fixed; background: #FFFFFF; border: 1px solid var(--toolbar-divider);
    padding: 8px 10px; font-size: 12px; border-radius: 4px;
    pointer-events: none; z-index: 20; max-width: 340px; line-height: 1.45;
    white-space: pre-wrap; display: none;
  }
</style>
<script>__LIB__</script>
<script type="application/json" id="cy-data">__PAYLOAD__</script>
</head>
<body>
<div id="toolbar">
  <div id="title">
    <h1 id="title-main"></h1>
    <div class="sub" id="title-sub"></div>
  </div>
  <div id="pickers">
    <div class="picker-row">
      <div class="lbl">Layout</div>
      <div class="buttons" id="layout-buttons"></div>
    </div>
    <div class="picker-row">
      <div class="lbl">Group by</div>
      <div class="buttons" id="theme-buttons"></div>
    </div>
    <div class="picker-row">
      <div class="lbl">Cofactors</div>
      <div class="buttons" id="cofactor-buttons"></div>
    </div>
    <div class="picker-row">
      <div class="lbl">Flags</div>
      <div class="buttons" id="flag-buttons"></div>
    </div>
  </div>
  <div id="actions">
    <button class="pick" id="fit-btn">Fit</button>
    <button class="pick" id="reset-btn">Reset</button>
    <input id="search" placeholder="Search name or id…">
  </div>
</div>
<div id="cy"></div>
<div id="toast"></div>
<div id="tooltip"></div>

<script>
(function () {
  const PAYLOAD = JSON.parse(document.getElementById('cy-data').textContent);
  const META = PAYLOAD.meta;

  const MATERIAL_BLUE = ['#E3F2FD','#BBDEFB','#90CAF9','#64B5F6','#42A5F5','#1E88E5','#1976D2','#1565C0','#0D47A1'];
  const DEEP_PURPLE   = ['#EDE7F6','#D1C4E9','#B39DDB','#9575CD','#7E57C2','#673AB7','#5E35B1','#512DA8','#4527A0'];
  const AMBER         = ['#FFF8E1','#FFECB3','#FFE082','#FFD54F','#FFCA28','#FFC107','#FFB300','#FFA000','#FF8F00'];
  function rampColor(pal, idx) {
    if (idx < 0 || !Number.isFinite(idx)) return pal[0];
    return pal[Math.min(idx, pal.length - 1)];
  }
  function numericPart(gid) {
    const m = /_(-?\d+)$/.exec(gid || '');
    return m ? parseInt(m[1], 10) : -1;
  }
  const PALETTE = {
    shell:          (gid) => rampColor(MATERIAL_BLUE, numericPart(gid)),
    step:           (gid) => rampColor(DEEP_PURPLE,   numericPart(gid)),
    producer_depth: (gid) => rampColor(AMBER,         numericPart(gid)),
    ec_class: (gid) => ({
      'ec_1': '#EF9A9A', 'ec_2': '#CE93D8', 'ec_3': '#9FA8DA',
      'ec_4': '#80DEEA', 'ec_5': '#80CBC4', 'ec_6': '#A5D6A7',
      'ec_other': '#BDBDBD',
    }[gid] || '#BDBDBD'),
    role: (gid) => ({
      'role_native': '#66BB6A',
      'role_intermediate': '#BBDEFB',
      'role_target': '#EF5350',
    }[gid] || '#ECEFF1'),
  };
  const NEUTRAL = { chemical: '#ECEFF1', reaction: '#FFF176' };
  const THEME_LABELS = {
    shell: 'Shell',
    ec_class: 'EC class',
    role: 'Role',
    producer_depth: 'Producer depth',
  };
  const THEME_RANK = {
    role: { 'role_native': 0, 'role_intermediate': 1, 'role_target': 2 },
    ec_class: { 'ec_1': 1, 'ec_2': 2, 'ec_3': 3, 'ec_4': 4, 'ec_5': 5, 'ec_6': 6, 'ec_other': 7 },
  };
  function groupRank(theme, gid) {
    if (!gid) return -1;
    if (THEME_RANK[theme] && THEME_RANK[theme][gid] !== undefined) return THEME_RANK[theme][gid];
    const n = numericPart(gid);
    return n >= 0 ? n : 0;
  }

  document.getElementById('title-main').textContent =
    META.view_kind === 'pathway'
      ? 'Pathway → ' + (META.target_name || '#' + META.target_chem_id)
      : 'Cascade → ' + (META.target_name || '#' + META.target_chem_id);
  document.getElementById('title-sub').textContent =
    META.counts.chemicals + ' chemicals · ' +
    META.counts.reactions + ' reactions · target shell ' + META.target_shell;

  // Keep originals pristine; every rebuild clones from these.
  const ORIGINAL_NODES = PAYLOAD.nodes;
  const ORIGINAL_EDGES = PAYLOAD.edges;

  let currentLayoutName = 'preset';
  let currentThemeName = META.default_theme;
  let currentCofactorMode = 'inline';
  let currentFocusId = null;   // id of the node currently in click-focus mode, or null
  let tapTimer = null;         // pending single-tap focus toggle (cancelled by dblclick)

  function currencyChemIdSet() {
    const s = new Set();
    for (const n of ORIGINAL_NODES) {
      if (n.data.kind === 'chemical' && n.data.is_currency && !n.data.is_target) {
        s.add(n.data.id);
      }
    }
    return s;
  }
  const HAS_CURRENCY = currencyChemIdSet().size > 0;

  // Expand PAYLOAD → an array of fresh element specs for the given cofactor
  // mode. "shared" is the raw graph; "hidden" drops currency chemicals and
  // their incident edges; "inline" clones a per-edge duplicate of each
  // currency chemical so no single cofactor node acts as a hub.
  function buildElementsForMode(mode) {
    const currencyIds = currencyChemIdSet();
    const out = [];
    if (mode === 'shared' || currencyIds.size === 0) {
      for (const n of ORIGINAL_NODES) out.push({ data: Object.assign({}, n.data) });
      for (const e of ORIGINAL_EDGES) out.push({ data: Object.assign({}, e.data) });
      return out;
    }
    if (mode === 'hidden') {
      for (const n of ORIGINAL_NODES) {
        if (n.data.kind === 'chemical' && currencyIds.has(n.data.id)) continue;
        out.push({ data: Object.assign({}, n.data) });
      }
      for (const e of ORIGINAL_EDGES) {
        if (currencyIds.has(e.data.source) || currencyIds.has(e.data.target)) continue;
        out.push({ data: Object.assign({}, e.data) });
      }
      return out;
    }
    // inline mode
    const rxnById = {};
    const chemById = {};
    for (const n of ORIGINAL_NODES) {
      if (n.data.kind === 'reaction') rxnById[n.data.id] = n;
      else if (n.data.kind === 'chemical') chemById[n.data.id] = n;
    }
    for (const n of ORIGINAL_NODES) {
      if (n.data.kind === 'chemical' && currencyIds.has(n.data.id)) continue;
      out.push({ data: Object.assign({}, n.data) });
    }
    for (const e of ORIGINAL_EDGES) {
      const srcCurrency = currencyIds.has(e.data.source);
      const tgtCurrency = currencyIds.has(e.data.target);
      if (!srcCurrency && !tgtCurrency) {
        out.push({ data: Object.assign({}, e.data) });
        continue;
      }
      const chemId = srcCurrency ? e.data.source : e.data.target;
      const rxnId  = srcCurrency ? e.data.target : e.data.source;
      const chemNode = chemById[chemId];
      const rxnNode  = rxnById[rxnId];
      if (!chemNode || !rxnNode) continue;
      const dupId = chemId + '__' + e.data.id;
      const dupData = Object.assign({}, chemNode.data);
      dupData.id = dupId;
      dupData.is_duplicate = true;
      // Layout-parent follows the reaction (so dup sits next to its rxn);
      // color still uses the chemical's own theme mapping (cofactors keep
      // their identity color wherever they appear).
      const parentMap = Object.assign({}, rxnNode.data.themes);
      for (const t in chemNode.data.themes) {
        if (parentMap[t] == null) parentMap[t] = chemNode.data.themes[t];
      }
      dupData.themes_parent = parentMap;
      out.push({ data: dupData });
      const newEdge = Object.assign({}, e.data);
      if (srcCurrency) {
        newEdge.source = dupId;
        newEdge.id = 'e_' + dupId + '_to_' + e.data.target;
      } else {
        newEdge.target = dupId;
        newEdge.id = 'e_' + e.data.source + '_to_' + dupId;
      }
      out.push({ data: newEdge });
    }
    return out;
  }

  // Pre-populate each element's theme-derived parent, display, and color
  // BEFORE handing them to Cytoscape, so the first paint already matches
  // the current theme (no flicker).
  function stampElementsForTheme(elements, themeName) {
    for (const el of elements) {
      const d = el.data;
      if (d.kind === 'group') {
        d.display = (d.theme === themeName) ? 'element' : 'none';
        continue;
      }
      const parentSource = d.themes_parent || d.themes;
      const parentGid = parentSource ? parentSource[themeName] : null;
      if (parentGid) d.parent = parentGid;
      else delete d.parent;
      const colorGid = d.themes ? d.themes[themeName] : null;
      d.color = d.is_target
        ? '#EF5350'
        : (colorGid && PALETTE[themeName] ? PALETTE[themeName](colorGid) : NEUTRAL[d.kind]);
    }
  }

  const initialElements = buildElementsForMode(currentCofactorMode);
  stampElementsForTheme(initialElements, currentThemeName);

  const cy = cytoscape({
    container: document.getElementById('cy'),
    elements: initialElements,
    wheelSensitivity: 0.2,
    style: [
      { selector: 'node[kind = "chemical"]',
        style: {
          'shape': 'round-rectangle',
          'background-color': 'data(color)',
          'border-width': 1.5,
          'border-color': '#78909C',
          'label': 'data(label)',
          'color': '#212121',
          'font-size': 12,
          'font-family': '-apple-system, "Helvetica Neue", Arial, sans-serif',
          'text-valign': 'center',
          'text-halign': 'center',
          'text-wrap': 'wrap',
          'text-max-width': 160,
          'padding': '6px',
          'width': 'label',
          'height': 'label',
        }},
      { selector: 'node[kind = "chemical"][?is_target]',
        style: {
          'border-width': 2.5,
          'border-color': '#C62828',
          'color': '#FFFFFF',
          'font-weight': 600,
        }},
      { selector: 'node[kind = "chemical"][?is_duplicate]',
        style: {
          'background-opacity': 0.6,
          'border-style': 'dashed',
          'border-color': '#B0BEC5',
          'border-width': 1,
          'font-size': 10,
        }},
      { selector: 'node[kind = "reaction"]',
        style: {
          'shape': 'hexagon',
          'background-color': 'data(color)',
          'border-width': 1.5,
          'border-color': '#F9A825',
          'label': 'data(label)',
          'color': '#212121',
          'font-size': 11,
          'font-family': '-apple-system, "Helvetica Neue", Arial, sans-serif',
          'text-valign': 'center',
          'text-halign': 'center',
          'text-wrap': 'wrap',
          'text-max-width': 140,
          'padding': '6px',
          'width': 'label',
          'height': 'label',
        }},
      { selector: 'node[kind = "group"]',
        style: {
          'shape': 'round-rectangle',
          'background-opacity': 0,
          'border-width': 1,
          'border-style': 'dashed',
          'border-color': '#BDBDBD',
          'label': 'data(label)',
          'color': '#616161',
          'font-size': 11,
          'font-weight': 600,
          'text-valign': 'top',
          'text-halign': 'left',
          'text-margin-y': -2,
          'padding': '14px',
        }},
      { selector: 'node[display = "none"]', style: { 'display': 'none' } },
      // Flag overlays. Each rule only fires when the user has toggled the
      // corresponding flag on — gated by the node-level `flag_*_on` data
      // attribute the JS sets in refreshFlags(). Priority (overridden later
      // in this array = higher): orphan → toxic → heme → p450.
      { selector: 'node[kind = "reaction"][?flag_orphan_on]',
        style: {
          'background-color': '#BDBDBD',
          'border-width': 2,
          'border-color': '#616161',
          'border-style': 'dashed',
          'color': '#212121',
        }},
      { selector: 'node[?flag_toxic_on]',
        style: {
          'border-width': 3,
          'border-color': '#EF6C00',
          'border-style': 'solid',
        }},
      { selector: 'node[kind = "reaction"][?flag_heme_on]',
        style: {
          'border-width': 3,
          'border-color': '#E53935',
          'border-style': 'solid',
        }},
      { selector: 'node[kind = "reaction"][?flag_p450_on]',
        style: {
          'border-width': 4,
          'border-color': '#B71C1C',
          'border-style': 'solid',
        }},
      { selector: 'node:selected',
        style: {
          'border-width': 3,
          'border-color': '#1E88E5',
        }},
      { selector: 'edge',
        style: {
          'width': 1.8,
          'line-color': '#9E9E9E',
          'curve-style': 'bezier',
          'target-arrow-shape': 'triangle',
          'target-arrow-color': '#9E9E9E',
          'arrow-scale': 0.9,
          'opacity': 0.85,
        }},
      { selector: 'edge[kind = "sub"]',
        style: {
          'target-arrow-shape': 'none',
        }},
    ],
    layout: { name: 'preset', fit: true, padding: 24 },
  });

  function computePresetPositions() {
    const positions = {};
    const byBucket = {};
    cy.nodes().not('[kind = "group"]').forEach(n => {
      const d = n.data();
      const source = d.themes_parent || d.themes;
      const gid = source ? source[currentThemeName] : null;
      const rank = groupRank(currentThemeName, gid);
      if (!byBucket[rank]) byBucket[rank] = [];
      byBucket[rank].push(n.id());
    });
    const ranks = Object.keys(byBucket).map(Number).sort((a, b) => a - b);
    const Y_SPACING = 140;
    const X_SPACING = 200;
    ranks.forEach((r, yi) => {
      const ids = byBucket[r].slice().sort();
      const width = Math.max(1, ids.length);
      const offset = -((width - 1) / 2) * X_SPACING;
      ids.forEach((id, xi) => {
        positions[id] = { x: offset + xi * X_SPACING, y: yi * Y_SPACING };
      });
    });
    return positions;
  }

  function layoutConfig(name) {
    const base = { fit: true, animate: false, padding: 24 };
    if (name === 'preset') {
      return Object.assign({}, base, { name: 'preset', positions: computePresetPositions() });
    }
    if (name === 'breadthfirst') {
      const roots = cy.nodes('[kind = "chemical"]')
        .filter(n => n.data('shell') === 0)
        .map(n => n.id());
      return Object.assign({}, base, {
        name: 'breadthfirst', directed: true, spacingFactor: 1.2,
        roots: roots.length ? roots : undefined,
      });
    }
    if (name === 'concentric') {
      return Object.assign({}, base, {
        name: 'concentric',
        concentric: (n) => (n.data('is_target') ? 999 : 100 - (n.data('shell') || 0)),
        levelWidth: () => 1,
        minNodeSpacing: 40,
      });
    }
    if (name === 'cose') {
      return Object.assign({}, base, { name: 'cose', idealEdgeLength: 90,
        nodeRepulsion: 8000, numIter: 1000 });
    }
    if (name === 'grid') {
      return Object.assign({}, base, { name: 'grid', avoidOverlap: true });
    }
    return Object.assign({}, base, { name });
  }

  function runLayout(name) {
    currentLayoutName = name;
    cy.layout(layoutConfig(name)).run();
    refreshLayoutButtons();
  }

  function applyTheme(themeName) {
    currentThemeName = themeName;
    cy.batch(() => {
      // 1) Re-parent content nodes to the new theme's groups first. This is
      //    important: if we hid groups first, any node parented under a
      //    display:none group would itself disappear.
      cy.nodes().not('[kind = "group"]').forEach(n => {
        const d = n.data();
        const parentSource = d.themes_parent || d.themes;
        const parentGid = parentSource ? parentSource[themeName] : null;
        n.move({ parent: parentGid || null });
        if (!d.is_target) {
          const colorGid = d.themes ? d.themes[themeName] : null;
          const color = colorGid && PALETTE[themeName]
            ? PALETTE[themeName](colorGid)
            : NEUTRAL[d.kind];
          n.data('color', color);
        }
      });
      // 2) Toggle group-node visibility.
      cy.nodes('[kind = "group"]').forEach(n => {
        n.style('display', n.data('theme') === themeName ? 'element' : 'none');
      });
    });
    runLayout(currentLayoutName);
    refreshThemeButtons();
  }

  function applyCofactorMode(mode) {
    currentCofactorMode = mode;
    // Any focus state refers to nodes that are about to be removed.
    if (tapTimer) { clearTimeout(tapTimer); tapTimer = null; }
    currentFocusId = null;
    const elements = buildElementsForMode(mode);
    stampElementsForTheme(elements, currentThemeName);
    cy.elements().remove();
    cy.add(elements);
    runLayout(currentLayoutName);
    refreshCofactorButtons();
    // Newly-added nodes need their flag_*_on data populated to match the
    // currently active flag set; otherwise the overlay silently drops after
    // a cofactor-mode switch.
    refreshFlags();
  }

  // Toolbar wiring
  const LAYOUTS = ['preset', 'breadthfirst', 'concentric', 'cose', 'grid'];
  const LAYOUT_LABELS = {
    preset: 'Layered', breadthfirst: 'Breadth-first',
    concentric: 'Concentric', cose: 'Force', grid: 'Grid',
  };
  const layoutBtnContainer = document.getElementById('layout-buttons');
  LAYOUTS.forEach(name => {
    const btn = document.createElement('button');
    btn.className = 'pick';
    btn.textContent = LAYOUT_LABELS[name];
    btn.dataset.layout = name;
    btn.addEventListener('click', () => runLayout(name));
    layoutBtnContainer.appendChild(btn);
  });
  function refreshLayoutButtons() {
    Array.from(layoutBtnContainer.children).forEach(b => {
      b.classList.toggle('active', b.dataset.layout === currentLayoutName);
    });
  }

  const themeBtnContainer = document.getElementById('theme-buttons');
  META.themes.forEach(theme => {
    const btn = document.createElement('button');
    btn.className = 'pick';
    btn.textContent = THEME_LABELS[theme] || theme;
    btn.dataset.theme = theme;
    btn.addEventListener('click', () => applyTheme(theme));
    themeBtnContainer.appendChild(btn);
  });
  function refreshThemeButtons() {
    Array.from(themeBtnContainer.children).forEach(b => {
      b.classList.toggle('active', b.dataset.theme === currentThemeName);
    });
  }

  const COFACTOR_MODES = ['inline', 'shared', 'hidden'];
  const COFACTOR_LABELS = {
    inline: 'Inline', shared: 'Shared hub', hidden: 'Hide',
  };
  const cofactorBtnContainer = document.getElementById('cofactor-buttons');
  COFACTOR_MODES.forEach(mode => {
    const btn = document.createElement('button');
    btn.className = 'pick';
    btn.textContent = COFACTOR_LABELS[mode];
    btn.dataset.cofactor = mode;
    btn.disabled = !HAS_CURRENCY;
    btn.addEventListener('click', () => applyCofactorMode(mode));
    cofactorBtnContainer.appendChild(btn);
  });
  if (!HAS_CURRENCY) {
    const note = document.createElement('span');
    note.style.fontSize = '11px';
    note.style.color = 'var(--subtle)';
    note.textContent = '(no cofactors detected)';
    cofactorBtnContainer.appendChild(note);
  }
  function refreshCofactorButtons() {
    Array.from(cofactorBtnContainer.children).forEach(b => {
      if (b.dataset && b.dataset.cofactor) {
        b.classList.toggle('active', b.dataset.cofactor === currentCofactorMode);
      }
    });
  }

  // ---- Flag overlay (P450 / Heme / Orphan / Toxic) ----
  // Each flag toggles independently and multiply-applies. The p450/heme/orphan
  // flags live on reaction nodes; toxic lives on chemicals (is_toxic) and
  // spreads to any reaction touching one (has_toxic_neighbor).
  const FLAGS = [
    { key: 'p450',   label: 'P450',   matches: (d) => d.kind === 'reaction' && !!d.is_p450 },
    { key: 'heme',   label: 'Heme',   matches: (d) => d.kind === 'reaction' && !!d.is_heme },
    { key: 'orphan', label: 'Orphan', matches: (d) => d.kind === 'reaction' && !!d.is_orphan },
    { key: 'toxic',  label: 'Toxic',
      matches: (d) => (d.kind === 'chemical' && !!d.is_toxic)
                   || (d.kind === 'reaction' && !!d.has_toxic_neighbor) },
  ];
  const activeFlags = new Set();
  const flagBtnContainer = document.getElementById('flag-buttons');

  function countFlag(flag) {
    let n = 0;
    for (const node of ORIGINAL_NODES) {
      if (node.data.kind === 'group') continue;
      if (flag.matches(node.data)) n += 1;
    }
    return n;
  }

  function refreshFlags() {
    // For each node currently in the graph, set its flag_<key>_on data to
    // true iff the node matches and the flag is active. Cytoscape style
    // selectors (declared in the stylesheet) pick these up automatically.
    cy.batch(() => {
      cy.nodes().not('[kind = "group"]').forEach(n => {
        const d = n.data();
        for (const flag of FLAGS) {
          const on = activeFlags.has(flag.key) && flag.matches(d);
          n.data('flag_' + flag.key + '_on', on);
        }
      });
    });
  }

  function refreshFlagButtons() {
    Array.from(flagBtnContainer.children).forEach(b => {
      if (!b.dataset || !b.dataset.flag) return;
      b.classList.toggle('active', activeFlags.has(b.dataset.flag));
    });
  }

  FLAGS.forEach(flag => {
    const n = countFlag(flag);
    const btn = document.createElement('button');
    btn.className = 'pick';
    btn.dataset.flag = flag.key;
    btn.disabled = (n === 0);
    btn.innerHTML = flag.label + '<span class="count">' + n + '</span>';
    btn.addEventListener('click', () => {
      if (btn.disabled) return;
      if (activeFlags.has(flag.key)) activeFlags.delete(flag.key);
      else activeFlags.add(flag.key);
      refreshFlags();
      refreshFlagButtons();
    });
    flagBtnContainer.appendChild(btn);
  });

  document.getElementById('fit-btn').addEventListener('click', () => cy.fit(null, 24));
  document.getElementById('reset-btn').addEventListener('click', () => {
    clearFocus();
    currentThemeName = META.default_theme;
    activeFlags.clear();
    applyCofactorMode('inline');
    refreshFlagButtons();
    runLayout('preset');
    cy.fit(null, 24);
  });

  // Search
  const searchInput = document.getElementById('search');
  searchInput.addEventListener('input', (e) => {
    // Typing in search supersedes any click-focus: clear it to keep state consistent.
    if (currentFocusId) clearFocus();
    const q = e.target.value.trim().toLowerCase();
    if (!q) {
      cy.elements().style('opacity', null);
      return;
    }
    const matches = cy.nodes().filter(n => {
      const d = n.data();
      if (d.kind === 'group') return false;
      return (d.label || '').toLowerCase().includes(q)
          || (d.name || '').toLowerCase().includes(q)
          || String(d.chem_id || d.rxn_id || '').includes(q);
    });
    cy.elements().style('opacity', 0.18);
    matches.style('opacity', 1);
    matches.connectedEdges().style('opacity', 0.7);
    matches.neighborhood('node').style('opacity', 0.9);
  });

  // Copy-on-dblclick
  function flagsLine(d) {
    const tags = [];
    if (d.is_p450) tags.push('P450');
    if (d.is_heme && !d.is_p450) tags.push('heme');
    if (d.is_orphan) tags.push('orphan');
    if (d.has_toxic_neighbor) tags.push('touches toxic');
    if (d.is_toxic) tags.push('toxic intermediate');
    return tags.length ? '\nFlags: [' + tags.join(', ') + ']' : '';
  }
  function formatForCopy(d) {
    if (d.kind === 'chemical') {
      const nameSuffix = d.name ? ' "' + d.name + '"' : '';
      const lines = ['Chemical #' + d.chem_id + nameSuffix + ' — shell=' + d.shell];
      if (d.inchi) lines.push('InChI=' + d.inchi);
      if (d.smiles) lines.push('SMILES=' + d.smiles);
      return lines.join('\n') + flagsLine(d);
    }
    if (d.kind === 'reaction') {
      const ec = d.ecnum ? ('EC ' + d.ecnum) : 'no EC';
      const enz = d.enzyme_name ? (' "' + d.enzyme_name + '"') : '';
      const subs = d.substrate_names.map((n, i) => '#' + d.substrate_ids[i] + ' ' + n).join(' + ');
      const prods = d.product_names.map((n, i) => '#' + d.product_ids[i] + ' ' + n).join(' + ');
      return 'Reaction #' + d.rxn_id + ' — ' + ec + enz + ' — shell=' + d.shell +
             '\nSubstrates: ' + subs +
             '\nProducts:   ' + prods +
             flagsLine(d);
    }
    return '';
  }

  // ---- Click-to-focus: dim everything except the clicked node + 1-hop neighborhood ----
  function applyFocus(nodeId) {
    const focal = cy.getElementById(nodeId);
    if (!focal || focal.empty() || focal.data('kind') === 'group') {
      clearFocus();
      return;
    }
    currentFocusId = nodeId;
    const hood = focal.closedNeighborhood(); // focal + neighbors + connecting edges
    const groups = cy.nodes('[kind = "group"]');
    cy.batch(() => {
      cy.elements().style('opacity', 0.12);
      groups.style('opacity', 1);           // keep group frames as context
      hood.style('opacity', 1);
      hood.edges().style({
        'line-color': '#1E88E5',
        'target-arrow-color': '#1E88E5',
        'width': 2.5,
      });
      focal.style({ 'border-width': 3, 'border-color': '#1E88E5' });
    });
  }

  function clearFocus() {
    currentFocusId = null;
    cy.batch(() => {
      cy.elements().style('opacity', null);
      cy.edges().removeStyle('line-color target-arrow-color width');
      cy.nodes().removeStyle('border-width border-color');
    });
  }

  cy.on('tap', 'node', (evt) => {
    const d = evt.target.data();
    if (d.kind === 'group') return;
    const nodeId = evt.target.id();
    // Defer the single-tap handler so a follow-up dblclick can cancel it
    // (dblclick should copy without a focus flicker).
    if (tapTimer) clearTimeout(tapTimer);
    tapTimer = setTimeout(() => {
      tapTimer = null;
      if (currentFocusId === nodeId) clearFocus();
      else applyFocus(nodeId);
    }, 220);
  });

  cy.on('tap', (evt) => {
    // Tap on empty background clears focus.
    if (evt.target === cy) {
      if (tapTimer) { clearTimeout(tapTimer); tapTimer = null; }
      clearFocus();
    }
  });

  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape') clearFocus();
  });

  cy.on('dblclick', 'node', (evt) => {
    // Cancel the pending single-tap focus toggle so dblclick-to-copy stays clean.
    if (tapTimer) { clearTimeout(tapTimer); tapTimer = null; }
    const d = evt.target.data();
    if (d.kind === 'group') return;
    const text = formatForCopy(d);
    if (!text) return;
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(
        () => showToast('Copied — paste back into Claude'),
        () => showToast('Copy failed; select + ⌘C manually')
      );
    } else {
      showToast('Clipboard API unavailable');
    }
  });

  const toastEl = document.getElementById('toast');
  let toastTimer = null;
  function showToast(msg) {
    toastEl.textContent = msg;
    toastEl.classList.add('show');
    if (toastTimer) clearTimeout(toastTimer);
    toastTimer = setTimeout(() => toastEl.classList.remove('show'), 2000);
  }

  const tooltipEl = document.getElementById('tooltip');
  cy.on('mouseover', 'node', (evt) => {
    const d = evt.target.data();
    if (d.kind === 'group') return;
    tooltipEl.textContent = formatForCopy(d);
    tooltipEl.style.display = 'block';
  });
  cy.on('mousemove', 'node', (evt) => {
    tooltipEl.style.left = (evt.originalEvent.clientX + 14) + 'px';
    tooltipEl.style.top  = (evt.originalEvent.clientY + 14) + 'px';
  });
  cy.on('mouseout', 'node', () => { tooltipEl.style.display = 'none'; });

  refreshLayoutButtons();
  refreshThemeButtons();
  refreshCofactorButtons();
  refreshFlagButtons();
  refreshFlags();
  runLayout('preset');
})();
</script>
</body>
</html>
"""


def _escape_html(s: str) -> str:
    return (
        s.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _render_html(payload: dict, title: str) -> str:
    try:
        lib_js = _CYTOSCAPE_JS.read_text(encoding="utf-8")
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"Cytoscape.js asset missing at {_CYTOSCAPE_JS}. "
            "Run `python scripts/fetch_cytoscape.py` to download it."
        ) from e
    payload_json = json.dumps(payload, separators=(",", ":"), ensure_ascii=False)
    # Prevent the JSON body from breaking out of the <script> tag.
    payload_json = payload_json.replace("</", "<\\/")
    return (
        _HTML_TEMPLATE
        .replace("__TITLE__", _escape_html(title))
        .replace("__LIB__", lib_js)
        .replace("__PAYLOAD__", payload_json)
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def render_pathway_html(
    pathway: Pathway,
    hg: HyperGraph,
    ec_to_name: Mapping[str, str] | None = None,
    chem_id_to_name: Mapping[int, str] | None = None,
    currency_chemicals: set[Chemical] | None = None,
    toxic_chem_ids: set[int] | None = None,
) -> str:
    chem_themes, rxn_themes, group_specs = _compute_themes_for_pathway(pathway, hg)
    chem_map = chem_id_to_name or {c.id: _chem_label_short(c) for c in pathway.metabolites}
    payload = _build_cytoscape_payload(
        reactions=pathway.reactions,
        target=pathway.target,
        hg=hg,
        ec_to_name=ec_to_name,
        chem_themes=chem_themes,
        rxn_themes=rxn_themes,
        group_specs=group_specs,
        view_kind="pathway",
        themes_in_view=PATHWAY_THEMES,
        chem_id_to_name=chem_map,
        currency_chemicals=currency_chemicals,
        toxic_chem_ids=toxic_chem_ids,
    )
    title = f"Pathway → {pathway.target.name or f'#{pathway.target.id}'}"
    return _render_html(payload, title)


def render_cascade_html(
    cascade: Cascade,
    hg: HyperGraph,
    ec_to_name: Mapping[str, str] | None = None,
    max_reactions: int = 200,
    chem_id_to_name: Mapping[int, str] | None = None,
    currency_chemicals: set[Chemical] | None = None,
    toxic_chem_ids: set[int] | None = None,
) -> str:
    n = len(cascade.reactions)
    if n > max_reactions:
        raise ValueError(
            f"Cascade has {n} reactions, exceeds max_reactions={max_reactions}. "
            "Lower shell_cutoff (toward -1), or raise max_reactions."
        )
    chem_themes, rxn_themes, group_specs = _compute_themes_for_cascade(cascade, hg)
    all_chems: set[Chemical] = {cascade.target}
    for r in cascade.reactions:
        all_chems.update(r.substrates)
        all_chems.update(r.products)
    chem_map = chem_id_to_name or {c.id: _chem_label_short(c) for c in all_chems}
    payload = _build_cytoscape_payload(
        reactions=cascade.reactions,
        target=cascade.target,
        hg=hg,
        ec_to_name=ec_to_name,
        chem_themes=chem_themes,
        rxn_themes=rxn_themes,
        group_specs=group_specs,
        view_kind="cascade",
        themes_in_view=CASCADE_THEMES,
        chem_id_to_name=chem_map,
        currency_chemicals=currency_chemicals,
        toxic_chem_ids=toxic_chem_ids,
    )
    title = f"Cascade → {cascade.target.name or f'#{cascade.target.id}'}"
    return _render_html(payload, title)
