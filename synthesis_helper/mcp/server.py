"""FastMCP server exposing synthesize / traceback / pathways to Claude Code."""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

from mcp.server.fastmcp import FastMCP

from synthesis_helper.composition import annotate_pathway
from synthesis_helper.models import Chemical
from synthesis_helper.pathways import enumerate_pathways
from synthesis_helper.traceback import traceback as build_cascade
from synthesis_helper.mcp import lookup, state
from synthesis_helper.mcp.serializers import (
    annotated_enzyme_to_dto,
    cascade_to_dto,
    chemical_to_dto,
    pathway_comparison_row_to_dto,
    pathway_to_dto,
    reaction_to_dto,
    similarity_hit_to_dto,
)
from synthesis_helper.mcp.html_render import render_cascade_html, render_pathway_html


mcp = FastMCP("synthesis-helper")


_MAX_PATHWAYS_CAP = int(os.environ.get("SYNTHESIS_MAX_PATHWAYS", "1000"))


def _resolve_or_raise(chemical_ref: str | int) -> Chemical:
    chems = state.get_chemicals()
    hit = lookup.resolve_one(chemical_ref, chems)
    if hit is None:
        raise ValueError(
            f"No chemical matches reference {chemical_ref!r}. "
            "Try an id, an exact name, or an InChI string."
        )
    return hit


# -----------------------------------------------------------------------------
# Tools
# -----------------------------------------------------------------------------


@mcp.tool()
def find_chemical(query: str, limit: int = 10) -> list[dict[str, Any]]:
    """Search chemicals by id, name (substring), or InChI.

    Returns up to *limit* matches as dicts with id, name, inchi, smiles, shell.
    Shell is null when the chemical is not reachable from the baseline cell.
    """
    hg = state.get_hypergraph()
    chems = state.get_chemicals()
    hits = lookup.resolve(query, chems, limit=limit)
    return [chemical_to_dto(c, hg) for c in hits]


@mcp.tool()
def find_by_structure(
    smiles: str,
    limit: int = 10,
    similarity_threshold: float = 0.4,
) -> list[dict[str, Any]]:
    """**Look up chemicals by structural similarity (SMILES input).**
    Call this whenever the user gives you a SMILES string, a drawn structure,
    or asks "find chemicals like this" / "what is this molecule in MetaCyc".
    Use instead of ``find_chemical`` when the query is a structure rather
    than a name, id, or InChI — ``find_chemical`` cannot match by shape.

    Returns top-N chemicals ranked by Tanimoto similarity on Morgan
    fingerprints (radius 2, 2048 bits). Each hit carries the usual chemical
    fields (id, name, inchi, smiles, shell) plus a ``tanimoto`` score in
    ``[0.0, 1.0]``. Hits below ``similarity_threshold`` are dropped; if
    fewer than ``limit`` hits pass, returns fewer (never an error).

    Raises ``ValueError`` on invalid SMILES.

    Note: the first call in a process builds a corpus-wide fingerprint index
    (~3-5 s for 9k chemicals); subsequent calls are fast. Coverage is close
    to 100% because chemicals without SMILES fall back to InChI parsing.
    """
    from synthesis_helper.mcp import similarity

    query_fp = similarity.fingerprint_from_smiles(smiles)
    index = state.get_fingerprint_index()
    chems = state.get_chemicals()
    hg = state.get_hypergraph()
    hits = similarity.tanimoto_search(query_fp, index, similarity_threshold, limit)
    return [similarity_hit_to_dto(chems[cid], score, hg) for cid, score in hits if cid in chems]


@mcp.tool()
def get_shell(chemical_ref: str | int) -> dict[str, Any]:
    """Return the shell (min reaction steps from native metabolites) of a chemical.

    shell=0 means native/universal; shell=null means not reachable.
    """
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    return {
        "id": chem.id,
        "name": chem.name,
        "shell": hg.chemical_to_shell.get(chem),
    }


@mcp.tool()
def list_reachables(
    max_shell: int | None = None,
    name_contains: str | None = None,
    limit: int = 100,
) -> list[dict[str, Any]]:
    """List reachable chemicals in the baseline hypergraph, filtered.

    Sorted by (shell, id). Default cap is 100 entries — bump *limit* for more.
    """
    hg = state.get_hypergraph()
    needle = name_contains.lower() if name_contains else None
    out = []
    for chem, shell in sorted(
        hg.chemical_to_shell.items(), key=lambda kv: (kv[1], kv[0].id)
    ):
        if max_shell is not None and shell > max_shell:
            continue
        if needle and needle not in chem.name.lower():
            continue
        out.append(chemical_to_dto(chem, hg))
        if len(out) >= limit:
            break
    return out


@mcp.tool()
def get_cascade(
    chemical_ref: str | int,
    shell_cutoff: int | None = None,
) -> dict[str, Any]:
    """Build the full cascade (all reactions reaching a target) from the baseline hypergraph.

    The cascade is the unflattened tree: every reaction that can produce the target,
    recursively, with shell-0 metabolites as leaves. Returned DTO includes
    reaction list, all metabolites touched, and the subset of leaf (shell-0)
    metabolites.

    *shell_cutoff* controls which producing reactions enter the cascade.
    The default ``None`` disables pruning: every producing reaction is
    admitted. Pass an integer ``n >= -1`` to admit a reaction only if all
    its substrates lie in shell <= (chemical's shell + n); ``n = -1`` is
    the tightest prune and keeps only producers whose substrates are
    strictly closer to shell 0 than the chemical itself.
    """
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    if chem not in hg.chemical_to_shell:
        raise ValueError(
            f"{chem.name!r} (id={chem.id}) is not reachable from the baseline cell."
        )
    cascade = build_cascade(hg, chem, shell_cutoff=shell_cutoff)
    return cascade_to_dto(cascade, hg)


@mcp.tool()
def enumerate_pathways_for(
    chemical_ref: str | int,
    max_pathways: int = 10,
    shell_cutoff: int | None = None,
) -> list[dict[str, Any]]:
    """Return up to *max_pathways* individual routes from native metabolites to target.

    Each pathway is a topologically ordered list of reactions whose combined
    output produces the target, starting entirely from shell-0 inputs.
    """
    capped = min(max_pathways, _MAX_PATHWAYS_CAP)
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    if chem not in hg.chemical_to_shell:
        raise ValueError(
            f"{chem.name!r} (id={chem.id}) is not reachable from the baseline cell."
        )
    cascade = build_cascade(hg, chem, shell_cutoff=shell_cutoff)
    pathways = enumerate_pathways(cascade, hg, max_pathways=capped)
    return [pathway_to_dto(p, i, hg) for i, p in enumerate(pathways)]


@mcp.tool()
def describe_pathway(
    chemical_ref: str | int,
    pathway_index: int,
    shell_cutoff: int | None = None,
) -> str:
    """Render a human-readable Markdown trace of a single pathway.

    Re-enumerates pathways (baseline hypergraph only) and picks index.
    """
    hg = state.get_hypergraph()
    chems = state.get_chemicals()
    chem = _resolve_or_raise(chemical_ref)
    cascade = build_cascade(hg, chem, shell_cutoff=shell_cutoff)
    pathways = enumerate_pathways(cascade, hg, max_pathways=pathway_index + 1)
    if pathway_index >= len(pathways):
        raise ValueError(
            f"Only {len(pathways)} pathway(s) found; index {pathway_index} out of range."
        )
    p = pathways[pathway_index]

    lines = [
        f"# Pathway {pathway_index} → {chem.name} (id={chem.id}, shell={hg.chemical_to_shell.get(chem)})",
        f"*{len(p.reactions)} reaction step(s)*",
        "",
    ]
    for step, rxn in enumerate(p.reactions, 1):
        subs = " + ".join(
            f"{chems[s.id].name}(#{s.id})" if s.id in chems else f"#{s.id}"
            for s in sorted(rxn.substrates, key=lambda c: c.id)
        )
        prods = " + ".join(
            f"{chems[p_.id].name}(#{p_.id})" if p_.id in chems else f"#{p_.id}"
            for p_ in sorted(rxn.products, key=lambda c: c.id)
        )
        ec = f"EC {rxn.ecnum}" if rxn.ecnum else f"rxn #{rxn.id}"
        lines.append(f"{step}. {subs} → {prods}  *({ec})*")
    return "\n".join(lines)


@mcp.tool()
def pathway_to_composition(
    chemical_ref: str | int,
    pathway_index: int = 0,
    shell_cutoff: int | None = None,
) -> dict[str, Any]:
    """**Machine-readable enzyme list for ONE pathway, with engineering flags.**
    Call this whenever the user asks for the enzyme list, composition, "which
    enzymes / which genes to express", or any per-step enzyme detail for a
    single pathway. Use instead of ``describe_pathway`` when the user wants
    structured data (for filtering, ranking, copying into a gene order) rather
    than a Markdown narrative.

    Returns ordered per-step rows: ``{step, ecnum, enzyme_name, is_orphan,
    is_p450, is_heme, ec_class, substrate_ids, product_ids}``.

    - ``is_orphan``: the reaction has no EC number OR the EC is missing from
      the EC → enzyme name table. These steps need curation before cloning.
    - ``is_p450``: EC starts with ``1.14.`` — cytochrome P450 / oxygenase
      family. Notoriously hard to express in *E. coli* (needs electron
      partners, often membrane-bound). Useful flag for "please give me a
      pathway without any P450s".
    - ``is_heme``: broader superset of ``is_p450`` — also flags peroxidases
      (EC ``1.11.1.*``) and dioxygenases (``1.13.11.*``). Any of these
      competes with the host for heme + heme-biosynthesis supply.
    - ``ec_class``: Oxidoreductase / Transferase / … / Unknown.

    Use ``compare_pathways`` first when the user wants to pick *which*
    pathway to compose; use this tool for the chosen one. Raises
    ``ValueError`` if ``pathway_index`` is out of range.
    """
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    if chem not in hg.chemical_to_shell:
        raise ValueError(
            f"{chem.name!r} (id={chem.id}) is not reachable from the baseline cell."
        )
    cascade = build_cascade(hg, chem, shell_cutoff=shell_cutoff)
    pathways = enumerate_pathways(cascade, hg, max_pathways=pathway_index + 1)
    if pathway_index >= len(pathways):
        raise ValueError(
            f"Only {len(pathways)} pathway(s) found; index {pathway_index} out of range."
        )
    pathway = pathways[pathway_index]
    enzymes = annotate_pathway(pathway, state.get_ec_names())
    return {
        "target": chemical_to_dto(chem, hg),
        "pathway_index": pathway_index,
        "step_count": len(enzymes),
        "enzymes": [annotated_enzyme_to_dto(e) for e in enzymes],
    }


@mcp.tool()
def compare_pathways(
    chemical_ref: str | int,
    n: int = 5,
    shell_cutoff: int | None = None,
) -> dict[str, Any]:
    """**Multi-dimensional scorecard across up to N pathways to one target.**
    Call this whenever the user asks to compare, rank, triage, or pick between
    pathways — "which pathway is most feasible", "avoid P450s", "simplest
    route". Use instead of calling ``pathway_to_composition`` N times.

    Returns a dict with a ``comparison`` list, sorted by
    ``(step_count, n_heme, n_toxic_intermediate, n_orphan_steps)`` so row 0
    is the first-pass most-viable candidate. Each row carries:

    - ``pathway_index`` — the original enumeration index (hand back to
      ``pathway_to_composition``, ``describe_pathway``, or
      ``open_pathway_interactive`` for drill-down)
    - ``step_count``
    - ``unique_ecs`` (set of EC strings), ``unique_ec_classes`` (EC-class
      diversity — higher usually means a more varied cascade)
    - ``n_orphan_steps`` (reactions with no / unmapped EC)
    - ``n_p450`` (EC starts with ``1.14.``) — narrow P450 flag.
    - ``n_heme`` — superset of ``n_p450``: also counts peroxidases
      (``1.11.1.*``) and dioxygenases (``1.13.11.*``). Any of these competes
      with the host for heme + heme-biosynthesis (ALA, hemA) supply, so this
      is a better "hard to express aerobically" signal than P450 alone.
    - ``n_toxic_intermediate`` — count of pathway reactions whose substrates
      or products touch a curated blacklist of reactive / toxic metabolites
      (reactive aldehydes, catechol / hydroquinone, acrylate, menadione, …).
      ``toxic_intermediates`` (list of canonical names) lets the caller see
      which ones were hit.
    - ``cofactor_uses`` — dict keyed by ATP / NADH / NADPH / NAD+ / NADP+ /
      CoA / SAM; value = reaction count consuming or producing that cofactor
      in this pathway (stoichiometry-agnostic: one reaction → one count).
    - ``nadh_hint`` / ``nadph_hint`` / ``atp_hint`` — crude net balance for
      each cofactor: +1 per reaction consuming (on substrate side), -1 per
      reaction regenerating (on product side). Positive = pathway drains
      that pool. NADH and NADPH pools are biologically separate (NADPH is
      biosynthetic, scarcer) so they're reported independently.

    ``cofactor_keys`` at the top level lists every key present in each row's
    ``cofactor_uses`` dict so callers can build tabular views reliably.
    """
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    if chem not in hg.chemical_to_shell:
        raise ValueError(
            f"{chem.name!r} (id={chem.id}) is not reachable from the baseline cell."
        )
    cap = max(1, min(n, _MAX_PATHWAYS_CAP))
    cascade = build_cascade(hg, chem, shell_cutoff=shell_cutoff)
    pathways = enumerate_pathways(cascade, hg, max_pathways=cap)

    ec_names = state.get_ec_names()
    cofactor_ids = state.get_named_cofactor_ids()
    toxic_map = state.get_toxic_intermediate_map()
    toxic_ids = frozenset(toxic_map.keys())
    nadh_ids = cofactor_ids["NADH"]
    nadph_ids = cofactor_ids["NADPH"]
    atp_ids = cofactor_ids["ATP"]

    def _metrics_for(pw) -> dict[str, Any]:
        ecs = {r.ecnum for r in pw.reactions if r.ecnum}
        ec_classes = {r.ecnum.split(".", 1)[0] for r in pw.reactions
                      if r.ecnum and r.ecnum[:1].isdigit()}
        n_orphan = sum(1 for r in pw.reactions if (not r.ecnum) or r.ecnum not in ec_names)
        n_p450 = sum(1 for r in pw.reactions if r.ecnum.startswith("1.14."))
        n_heme = sum(
            1 for r in pw.reactions
            if r.ecnum.startswith("1.14.")
            or r.ecnum.startswith("1.11.1.")
            or r.ecnum.startswith("1.13.11.")
        )
        cofactor_uses: dict[str, int] = {key: 0 for key in state.COFACTOR_KEYS}
        nadh_hint = 0
        nadph_hint = 0
        atp_hint = 0
        n_toxic = 0
        toxic_seen: set[str] = set()
        for rxn in pw.reactions:
            sub_ids = {s.id for s in rxn.substrates}
            prod_ids = {p.id for p in rxn.products}
            touching = sub_ids | prod_ids
            for key in state.COFACTOR_KEYS:
                if touching & cofactor_ids[key]:
                    cofactor_uses[key] += 1
            # Net direction per pool: +1 if consumed (substrate), -1 if
            # regenerated (product). Stoichiometry-agnostic.
            if sub_ids & nadh_ids:
                nadh_hint += 1
            if prod_ids & nadh_ids:
                nadh_hint -= 1
            if sub_ids & nadph_ids:
                nadph_hint += 1
            if prod_ids & nadph_ids:
                nadph_hint -= 1
            if sub_ids & atp_ids:
                atp_hint += 1
            if prod_ids & atp_ids:
                atp_hint -= 1
            hit_toxic = touching & toxic_ids
            if hit_toxic:
                n_toxic += 1
                for cid in hit_toxic:
                    toxic_seen.add(toxic_map[cid])
        return {
            "unique_ecs": sorted(ecs),
            "unique_ec_classes": sorted(ec_classes),
            "n_orphan_steps": n_orphan,
            "n_p450": n_p450,
            "n_heme": n_heme,
            "n_toxic_intermediate": n_toxic,
            "toxic_intermediates": sorted(toxic_seen),
            "cofactor_uses": cofactor_uses,
            "nadh_hint": nadh_hint,
            "nadph_hint": nadph_hint,
            "atp_hint": atp_hint,
        }

    rows: list[tuple[int, Any, dict[str, Any]]] = [
        (i, pw, _metrics_for(pw)) for i, pw in enumerate(pathways)
    ]
    rows.sort(key=lambda t: (
        len(t[1].reactions),
        t[2]["n_heme"],
        t[2]["n_toxic_intermediate"],
        t[2]["n_orphan_steps"],
    ))
    comparison = [
        pathway_comparison_row_to_dto(i, pw, metrics, hg) for (i, pw, metrics) in rows
    ]
    return {
        "target": chemical_to_dto(chem, hg),
        "n_pathways_considered": len(pathways),
        "cofactor_keys": list(state.COFACTOR_KEYS),
        "comparison": comparison,
    }


# -----------------------------------------------------------------------------
# Interactive HTML viewer helpers
# -----------------------------------------------------------------------------


_FILENAME_SAFE = re.compile(r"[^\w\-]+")


def _slug(s: str, max_len: int = 40) -> str:
    cleaned = _FILENAME_SAFE.sub("_", s).strip("_")[:max_len]
    return cleaned or "x"


def _open_in_viewer(path: Path) -> tuple[bool, str | None]:
    """Launch the OS default app for *path*. Returns (ok, error_message)."""
    try:
        if sys.platform == "darwin":
            cmd = ["open", str(path)]
        elif sys.platform.startswith("linux"):
            cmd = ["xdg-open", str(path)]
        elif sys.platform == "win32":
            cmd = ["cmd", "/c", "start", "", str(path)]
        else:
            return False, f"unsupported platform: {sys.platform}"
        subprocess.run(cmd, check=True, timeout=10)
        return True, None
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return False, str(e)


def _write_html(prefix: str, html: str) -> Path:
    fd, tmp = tempfile.mkstemp(prefix=prefix, suffix=".html", dir=tempfile.gettempdir())
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write(html)
    return Path(tmp)


_DO_NOT_VERIFY = (
    "The HTML file and browser tab are on the user's local machine, which is "
    "NOT visible to Bash / Read tools in your sandbox — do not attempt to "
    "verify. Do NOT re-call this tool for the same request; every call "
    "opens another browser tab."
)


def _viewer_result(
    *,
    path: Path,
    opened: bool,
    viewer_error: str | None,
    html_bytes: int,
    extra: dict[str, Any],
) -> dict[str, Any]:
    """Build a consistent return shape for the interactive-viewer tools.

    ``status`` + ``message`` lead the dict so the LLM has an unambiguous
    success/failure signal before scanning the rest of the payload — this
    prevents the "retry to verify" loop that causes duplicate browser tabs.
    """
    if opened:
        status = "opened"
        message = (
            f"Pathway viewer launched in the user's default browser "
            f"(file: {path}). {_DO_NOT_VERIFY}"
        )
    else:
        status = "write_only"
        message = (
            f"HTML written to {path} but the system viewer failed to launch: "
            f"{viewer_error}. The file is on the user's machine; they can "
            f"open it manually. Do not re-call this tool automatically."
        )
    result: dict[str, Any] = {
        "status": status,
        "message": message,
        "path": str(path),
        "opened_in_viewer": opened,
        "viewer_error": viewer_error,
        "html_bytes": html_bytes,
    }
    result.update(extra)
    return result


@mcp.tool()
def open_pathway_interactive(
    chemical_ref: str | int,
    pathway_index: int = 0,
    shell_cutoff: int | None = None,
) -> dict[str, Any]:
    """**THE visualization tool for a single pathway.** Call this whenever the
    user asks to see, draw, view, graph, diagram, or visualize one pathway to
    a target chemical — there is no PNG / text-graph alternative; this is it.

    **SIDE EFFECT — READ BEFORE CALLING TWICE:** This tool WRITES an HTML
    file to the user's local ``/tmp`` and LAUNCHES their default browser,
    opening a new tab. Both actions happen on the user's machine (where
    the MCP server runs), NOT in the sandbox where your Bash / Read tools
    operate. The file path in the return value CANNOT be verified with
    Bash or Read — any such attempt will fail even though the tool
    succeeded. **Trust ``status == "opened"`` and do not re-call the tool
    to "check" or retry — every call opens a new browser tab.**

    Renders pathway #``pathway_index`` (0 = shortest / first enumerated) as a
    self-contained interactive HTML page (Cytoscape.js inlined, works offline)
    and launches the user's default browser. The page gives the user:

    - 5 layout algorithms (Layered / Breadth-first / Concentric / Force / Grid)
    - Group-by themes: shell, ec_class, role
    - Cofactor mode toggle (Inline / Shared / Hide) — defaults to Inline,
      which duplicates H2O / NAD(P)H / ATP / … next to each reaction so the
      carbon backbone stands out instead of a hairball of cofactor edges
    - Click any node to dim everything except its 1-hop neighborhood (ESC or
      background-click clears)
    - Double-click any node to copy its structured details to the clipboard
      for pasting back into Claude

    Returns a dict with ``status``, ``message``, the file path, step count,
    etc. Use ``open_cascade_interactive`` instead when the user wants the
    *full tree* of ways to make the target (branching producer graph), not
    a single linear pathway.
    """
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    if chem not in hg.chemical_to_shell:
        raise ValueError(
            f"{chem.name!r} (id={chem.id}) is not reachable from the baseline cell."
        )
    cascade = build_cascade(hg, chem, shell_cutoff=shell_cutoff)
    pathways = enumerate_pathways(cascade, hg, max_pathways=pathway_index + 1)
    if pathway_index >= len(pathways):
        raise ValueError(
            f"Only {len(pathways)} pathway(s) found; index {pathway_index} out of range."
        )
    chem_id_to_name = {
        c.id: (c.name if c.name and c.name != "undefined" else f"#{c.id}")
        for c in state.get_chemicals().values()
    }
    html = render_pathway_html(
        pathways[pathway_index], hg, state.get_ec_names(),
        chem_id_to_name=chem_id_to_name,
        currency_chemicals=state.get_currency_chemicals(),
        toxic_chem_ids=set(state.get_toxic_intermediate_map().keys()),
    )
    path = _write_html(
        f"pathway_{_slug(chem.name or f'chem{chem.id}')}_{pathway_index}_", html
    )
    opened, err = _open_in_viewer(path)
    return _viewer_result(
        path=path,
        opened=opened,
        viewer_error=err,
        html_bytes=len(html),
        extra={
            "chemical": chem.name,
            "chemical_id": chem.id,
            "pathway_index": pathway_index,
            "steps": len(pathways[pathway_index].reactions),
        },
    )


@mcp.tool()
def open_cascade_interactive(
    chemical_ref: str | int,
    shell_cutoff: int | None = None,
    max_reactions: int = 200,
) -> dict[str, Any]:
    """**THE visualization tool for a full cascade / producer tree.** Call this
    whenever the user asks to see, draw, view, graph, diagram, or visualize
    *all* the ways to synthesize a target (not one linear route) — there is
    no PNG / text-graph alternative; this is it.

    **SIDE EFFECT — READ BEFORE CALLING TWICE:** This tool WRITES an HTML
    file to the user's local ``/tmp`` and LAUNCHES their default browser,
    opening a new tab. Both actions happen on the user's machine (where
    the MCP server runs), NOT in the sandbox where your Bash / Read tools
    operate. The file path in the return value CANNOT be verified with
    Bash or Read — any such attempt will fail even though the tool
    succeeded. **Trust ``status == "opened"`` and do not re-call the tool
    to "check" or retry — every call opens a new browser tab.**

    Renders the complete cascade (every reaction that can produce the target,
    recursively back to native metabolites) as a self-contained interactive
    HTML page (Cytoscape.js inlined, works offline) and launches the user's
    default browser. Same page interactions as ``open_pathway_interactive``
    (layouts, cofactor toggle, click-to-focus, double-click-to-copy); the
    group-by themes are {shell, ec_class, role, producer_depth}.

    Large cascades can be unwieldy — if the target has many producer
    reactions, set ``shell_cutoff`` to an integer (default ``None`` is
    unlimited; lower toward ``-1`` to prune harder). ``max_reactions``
    is a safety ceiling (default 200); the call raises if the cascade
    exceeds it so the browser doesn't choke on a thousand-node graph.

    Use ``open_pathway_interactive`` instead when the user wants one specific
    linear pathway to the target.
    """
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    if chem not in hg.chemical_to_shell:
        raise ValueError(
            f"{chem.name!r} (id={chem.id}) is not reachable from the baseline cell."
        )
    cascade = build_cascade(hg, chem, shell_cutoff=shell_cutoff)
    chem_id_to_name = {
        c.id: (c.name if c.name and c.name != "undefined" else f"#{c.id}")
        for c in state.get_chemicals().values()
    }
    html = render_cascade_html(
        cascade, hg, state.get_ec_names(),
        max_reactions=max_reactions, chem_id_to_name=chem_id_to_name,
        currency_chemicals=state.get_currency_chemicals(),
        toxic_chem_ids=set(state.get_toxic_intermediate_map().keys()),
    )
    path = _write_html(f"cascade_{_slug(chem.name or f'chem{chem.id}')}_", html)
    opened, err = _open_in_viewer(path)
    return _viewer_result(
        path=path,
        opened=opened,
        viewer_error=err,
        html_bytes=len(html),
        extra={
            "chemical": chem.name,
            "chemical_id": chem.id,
            "reactions": len(cascade.reactions),
        },
    )


@mcp.tool()
def resynthesize_with_fed(
    fed_chemical_refs: list[str | int],
    include_new_reachables: bool = True,
    new_reachables_limit: int = 50,
) -> dict[str, Any]:
    """Re-run synthesize() with additional fed chemicals added to shell 0.

    Results are cached per fed-set so repeated calls with the same inputs
    are fast. Returns summary stats plus (optionally) the list of chemicals
    that became reachable *only* with the fed additions, up to *new_reachables_limit*.
    """
    chems = state.get_chemicals()
    fed: list[Chemical] = []
    unresolved: list[Any] = []
    for ref in fed_chemical_refs:
        hit = lookup.resolve_one(ref, chems)
        if hit is None:
            unresolved.append(ref)
        else:
            fed.append(hit)

    fed_ids = tuple(sorted({c.id for c in fed}))
    baseline = state.get_hypergraph()
    fed_hg = state.get_or_build_fed_hypergraph(fed_ids)

    shells = max(fed_hg.reaction_to_shell.values(), default=0)
    new_reachables: list[dict[str, Any]] = []
    if include_new_reachables:
        diff = [
            c for c in fed_hg.chemical_to_shell
            if c not in baseline.chemical_to_shell
        ]
        diff.sort(key=lambda c: (fed_hg.chemical_to_shell[c], c.id))
        new_reachables = [chemical_to_dto(c, fed_hg) for c in diff[:new_reachables_limit]]

    return {
        "fed_chemicals": [chemical_to_dto(c, fed_hg) for c in fed],
        "unresolved": [str(u) for u in unresolved],
        "shells": shells,
        "reachable_chemicals": len(fed_hg.chemical_to_shell),
        "reachable_reactions": len(fed_hg.reaction_to_shell),
        "baseline_reachable_chemicals": len(baseline.chemical_to_shell),
        "delta_chemicals": len(fed_hg.chemical_to_shell) - len(baseline.chemical_to_shell),
        "delta_reactions": len(fed_hg.reaction_to_shell) - len(baseline.reaction_to_shell),
        "new_reachables": new_reachables,
        "new_reachables_truncated": (
            include_new_reachables
            and (len(fed_hg.chemical_to_shell) - len(baseline.chemical_to_shell))
            > new_reachables_limit
        ),
    }


# -----------------------------------------------------------------------------
# Resources
# -----------------------------------------------------------------------------


@mcp.resource("synthesis://stats")
def stats_resource() -> dict[str, Any]:
    """Summary stats of the loaded hypergraph."""
    return state.get_stats()


@mcp.resource("synthesis://reachables")
def reachables_resource() -> str:
    """TSV dump of all reachables (id, name, inchi, shell)."""
    hg = state.get_hypergraph()
    lines = ["id\tname\tinchi\tshell"]
    for chem, shell in sorted(
        hg.chemical_to_shell.items(), key=lambda kv: (kv[1], kv[0].id)
    ):
        lines.append(f"{chem.id}\t{chem.name}\t{chem.inchi}\t{shell}")
    return "\n".join(lines) + "\n"


@mcp.resource("synthesis://chemical/{chem_id}")
def chemical_resource(chem_id: str) -> dict[str, Any]:
    """Chemical detail + which reactions produce / consume it."""
    hg = state.get_hypergraph()
    chems = state.get_chemicals()
    try:
        cid = int(chem_id)
    except ValueError:
        raise ValueError(f"chemical_id must be integer, got {chem_id!r}")
    if cid not in chems:
        raise ValueError(f"Unknown chemical id {cid}")
    chem = chems[cid]

    produced_by = [r.id for r in hg.reaction_to_shell if chem in r.products]
    consumed_by = [r.id for r in hg.reaction_to_shell if chem in r.substrates]
    return {
        **chemical_to_dto(chem, hg),
        "produced_by_reaction_ids": sorted(produced_by),
        "consumed_by_reaction_ids": sorted(consumed_by),
    }


@mcp.resource("synthesis://reaction/{rxn_id}")
def reaction_resource(rxn_id: str) -> dict[str, Any]:
    """Reaction detail with expanded substrate / product names."""
    hg = state.get_hypergraph()
    try:
        rid = int(rxn_id)
    except ValueError:
        raise ValueError(f"reaction_id must be integer, got {rxn_id!r}")
    for r in state.get_reactions():
        if r.id == rid:
            return {
                **reaction_to_dto(r, hg),
                "substrates": [chemical_to_dto(s, hg) for s in sorted(r.substrates, key=lambda c: c.id)],
                "products": [chemical_to_dto(p, hg) for p in sorted(r.products, key=lambda c: c.id)],
            }
    raise ValueError(f"Unknown reaction id {rid}")


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------


def main() -> None:
    """Run the server on stdio (what Claude Code connects to)."""
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
