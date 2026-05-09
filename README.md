# synthesis_helper

A bottom-up retrobiosynthesis tool that builds a synthesis hypergraph from a
reaction corpus (MetaCyc) and enumerates **cascades, pathways, and enzyme
compositions** for any reachable target chemical. The whole thing is wrapped
as a [Model Context Protocol](https://modelcontextprotocol.io/) server so
Claude Code can drive it as an interactive lab notebook.

This document is organized around the two pieces of work that took the most
effort:

1. **What the MCP layer does and how it is implemented.** — *Benson Li*
2. **How we validated that `synthesize` and the cascade traceback are correct.** — *Benson Li*

The underlying algorithm (BFS waveform expansion → cascade traceback →
pathway enumeration → enzyme composition) is documented in `CLAUDE.md`.

> **Status:** sections 1 and 2 below were written and implemented by
> Benson Li. This file is being handed to the rest of the team — please
> append your sections (algorithm internals, data preparation, downstream
> integration, etc.) below section 2 and add yourselves to the relevant
> headers.

---

## 1. The MCP layer
*Author: Benson Li*

### 1.1 What the server exposes

The `synthesis_helper.mcp` package wraps the core API as a FastMCP server
over stdio. Once registered with Claude Code, the model gets the following
**tools** and **resources**.

#### Tools

| Tool | What it does |
|---|---|
| `find_chemical(query, limit)` | Look up a chemical by id, name substring, or InChI. |
| `find_by_structure(smiles, limit, similarity_threshold)` | Tanimoto search over Morgan fingerprints (radius 2, 2048 bits) for "find me chemicals shaped like this SMILES". |
| `get_shell(chemical_ref)` | Minimum number of reaction steps from native metabolites. `0` = native/universal, `null` = unreachable. |
| `list_reachables(max_shell, name_contains, limit)` | Filter the reachables table. |
| `get_cascade(chemical_ref, max_producers_per_chemical)` | The full producer tree (all reactions transitively reaching the target). |
| `enumerate_pathways_for(chemical_ref, max_pathways)` | Flatten the cascade into individual linear pathways. |
| `describe_pathway(chemical_ref, pathway_index)` | Markdown narrative of a single pathway. |
| `pathway_to_composition(chemical_ref, pathway_index)` | Per-step enzyme list with engineering flags (`is_orphan`, `is_p450`, `is_heme`, `ec_class`). |
| `compare_pathways(chemical_ref, n)` | Multi-dimensional scorecard ranking up to N pathways by step count, heme dependence, toxic intermediates, orphan steps, and net cofactor balance (NADH / NADPH / ATP). |
| `open_pathway_interactive(chemical_ref, pathway_index)` | Render one pathway as a self-contained interactive HTML page (Cytoscape.js, offline) and open it in the browser. |
| `open_cascade_interactive(chemical_ref, max_reactions)` | Same, for the full cascade. |
| `resynthesize_with_fed(fed_chemical_refs)` | Re-run BFS with extra fed chemicals added to shell 0; returns the delta vs. baseline. |

#### Resources

| URI | Contents |
|---|---|
| `synthesis://stats` | Counts, load time, data directory. |
| `synthesis://reachables` | TSV dump of every reachable chemical. |
| `synthesis://chemical/{id}` | Chemical detail + `produced_by` / `consumed_by` reaction ids. |
| `synthesis://reaction/{id}` | Reaction with expanded substrate / product names. |

### 1.2 How it is implemented

The MCP layer is split into a thin server and a few helpers so the heavy
work happens once per process and the tool functions stay shallow.

```
synthesis_helper/mcp/
├── server.py        # FastMCP @tool / @resource definitions
├── state.py         # Lazy singletons: chemicals, reactions, hypergraph,
│                    #   fingerprint index, EC-name table, cofactor sets,
│                    #   toxic-intermediate map, fed-hypergraph cache
├── lookup.py        # resolve_one() — id / exact-name / InChI / fuzzy match
├── similarity.py    # RDKit Morgan fingerprints + Tanimoto top-N search
├── serializers.py   # Domain object → JSON-friendly DTOs
├── html_render.py   # Cytoscape.js graph generation for the viewer tools
└── assets/
    └── cytoscape.min.js   # Inlined offline copy (fetched on first install)
```

Key design choices:

- **Lazy singleton state.** The hypergraph takes a few seconds to build from
  ~9k chemicals and ~9k reactions. `state.py` builds it on the first tool
  call and caches it in process memory; every subsequent tool call is
  effectively a dictionary lookup. The fingerprint index for
  `find_by_structure` follows the same pattern (built only when the first
  structure query arrives).
- **Fed-hypergraph cache.** `resynthesize_with_fed` keys its cache on the
  sorted tuple of fed chemical ids, so re-asking the same "what if I feed
  PABA?" question costs nothing after the first call.
- **Reference resolution at the boundary.** Every tool that takes a
  `chemical_ref` runs it through `lookup.resolve_one`, which accepts
  integer ids, exact names, InChI strings, or substring matches with a
  single confident hit. This means callers (and Claude) never have to
  pre-canonicalize identifiers.
- **DTO serializers, not raw dataclasses.** `serializers.py` converts
  `Chemical` / `Reaction` / `Cascade` / `Pathway` into plain dicts with
  shell numbers attached. The MCP wire format is JSON; keeping that
  conversion in one place means the domain objects can evolve without
  breaking the MCP contract.
- **Offline interactive viewer.** `open_pathway_interactive` and
  `open_cascade_interactive` write a self-contained HTML file (Cytoscape.js
  is inlined) to `/tmp` and shell out to `open` / `xdg-open` / `start` to
  launch the user's default browser. The page supports five layout
  algorithms (Layered, Breadth-first, Concentric, Force, Grid), three
  group-by themes (shell, EC class, role), a cofactor mode toggle (Inline /
  Shared / Hide), click-to-focus on a 1-hop neighborhood, and double-click
  to copy a node's structured details to the clipboard. All of this runs
  client-side — no server round-trip on layout / theme switch.
- **Engineering-aware pathway scoring.** `compare_pathways` is more than a
  step counter. It flags P450s (`EC 1.14.*`), broader heme-dependent
  enzymes (`1.11.1.*`, `1.13.11.*`), orphan steps with no EC mapping, and
  reactions that touch a curated blacklist of toxic / reactive
  intermediates (catechol, hydroquinone, reactive aldehydes, …). It also
  computes a crude net balance for ATP / NADH / NADPH so the user can spot
  pathways that drain biosynthetic reducing power.
- **Side-effect honesty in tool docstrings.** The `open_*_interactive`
  tools spell out, in the docstring the LLM sees, that they write a file
  and open a browser tab on the user's machine — and that retrying to
  "verify" the result will just open another tab. This is a small fix but
  a critical one for usability when the model can't see the user's
  desktop.

### 1.3 Install and register

```bash
python3 -m venv ~/.venvs/synthesis_helper
~/.venvs/synthesis_helper/bin/pip install -r requirements.txt
~/.venvs/synthesis_helper/bin/python scripts/fetch_cytoscape.py   # offline asset
```

> **macOS + iCloud note.** If the project lives under `~/Documents` or
> `~/Desktop` with iCloud "Optimize Mac Storage" on, files inside `.venv/`
> can get evicted to cloud stubs and Python imports stall for minutes
> during the MCP handshake. Keep the venv outside any synced folder
> (`~/.venvs/...` is fine).

Register with Claude Code (preferred — the launcher script materializes any
evicted source files and finds the venv automatically):

```bash
claude mcp add synthesis-helper --scope user -- /absolute/path/to/synthesis_helper/scripts/start_mcp.sh
```

Verify with `/mcp` inside Claude Code — you should see
`synthesis-helper: connected`.

#### Environment variables

- `SYNTHESIS_DATA_DIR` — point at a different `data/` directory (default:
  the repo's own `data/`).
- `SYNTHESIS_MAX_PATHWAYS` — hard cap on `enumerate_pathways_for` (default
  `1000`).
- `SYNTHESIS_HELPER_VENV` — override the venv path used by
  `scripts/start_mcp.sh`.

#### Optional: refresh the EC-name table

`data/ec_names.tsv` ships with a full EC → enzyme-name table derived from
[ExPASy ENZYME](https://enzyme.expasy.org/) (8k+ entries, transferred
entries resolved). To rebuild:

```bash
curl -sS https://ftp.expasy.org/databases/enzyme/enzyme.dat -o /tmp/enzyme.dat
.venv/bin/python scripts/build_ec_names.py /tmp/enzyme.dat
```

The file is redistributed under the same
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license as the
upstream ENZYME database (© SIB Swiss Institute of Bioinformatics).

### 1.4 Example prompts

- "How many steps from native E. coli metabolites to vanillin?" →
  `get_shell`.
- "List the top 3 pathways to 2′-dehydrokanamycin a, ranked by feasibility,
  and avoid P450s." → `compare_pathways` + `pathway_to_composition`.
- "If I feed PABA into the cell, what becomes newly reachable?" →
  `resynthesize_with_fed`.
- "Open an interactive cascade view for chemical 317157." →
  `open_cascade_interactive`.

---

## 2. Validating `synthesize` and the cascade traceback
*Author: Benson Li*

The hard part of any retrobiosynthesis tool is convincing yourself the
algorithm is right. Empirical correctness on a 9k-reaction MetaCyc dump is
hard to inspect by eye, so we built a layered validation strategy: a
hand-drawn synthetic cell where every shell number is computable on paper,
implementation-agnostic invariant tests, and a real-world cross-reference
against a published MetaCyc pathway.

The full test suite is **114 tests** across 10 files, run with
`.venv/bin/pytest tests/ -q`.

### 2.1 The hand-drawn small cell (`validation/small_cell.md`)

A purely synthetic 12-chemical, 11-reaction "cell" defined inline in
`tests/fixtures/small_cell.py`. It does **not** read anything under
`data/`, so changes to the real MetaCyc corpus cannot affect these tests
(and conversely, a regression in the algorithm cannot hide behind corpus
churn).

The fixture is engineered to exercise every structural scenario the BFS
must handle:

| Scenario | What it pins down |
|---|---|
| (a) Linear chain N1 → M1 → M2 → M_deep | Shells propagate correctly past depth 2. |
| (b) Multi-substrate `R3: N1 + N2 → M4` | Reaction shell = `max(substrate shell) + 1`. |
| (c) Universal cofactor in `R4: M4 + U1 → M5` | A shell-0 universal cannot inflate a downstream shell. |
| (d) Multi-producer target T (via R7 *or* R6) | Cascade contains both routes; pathway enumeration must surface both. |
| (e) Unreachable target `M7 ← X_missing` | BFS stops cleanly; traceback raises `ValueError`. |
| (f) Tied-shell producers of M3b (R8 at shell 1, R9 at shell 2) | The shorter route wins. |
| (g) Fed chemical `F1 ← R10 → M6` | `resynthesize_with_fed` enables the reaction without re-running the whole pipeline. |

Every expected shell, reaction, cascade member, and pathway is tabulated
in `validation/small_cell.md` and asserted by tests in
`tests/test_small_cell_synthesize.py`,
`tests/test_small_cell_traceback.py`, and
`tests/test_small_cell_pathways.py`.

This fixture also surfaces a real bug: in `pathways.py:78–80`,
`enumerate_pathways` only recurses on the first non-native substrate of a
multi-substrate producer reaction. The test
`test_pathways_enumerates_both_complete_routes` is decorated
`@pytest.mark.xfail(strict=True)` so the bug stays visible — when fixed,
the test flips to XPASS and the suite goes red, forcing a conscious
removal of the marker.

### 2.2 Implementation-agnostic invariants (`tests/test_invariants.py`)

The small cell tests pin specific values; the invariant tests pin the
**meaning** of `synthesize` and `traceback`. They apply to any correct
implementation, so a teammate could swap the BFS for a totally different
algorithm and these tests would still catch a spec violation.

`synthesize()` invariants asserted:

- Every native and universal sits at shell 0.
- Shells are non-negative integers; reactions start at 1, chemicals at 0.
- **BFS defining equation:** `shell(R) == max(shell(s) for s in R.substrates) + 1`.
- Every product of an enabled reaction is reachable (and likewise for
  substrates).
- Every non-native chemical has a producer at the matching shell.
- **Fixed-point completeness:** no unchosen reaction has all substrates
  already reachable — i.e., the BFS ran to closure.
- Every reachable chemical can be successfully tracebacked.

Cascade (traceback) invariants asserted:

- Cascade reactions ⊆ enabled reactions (no phantoms).
- No orphan reactions in the cascade — every reaction produces something
  that contributes to the target.
- **Closure under substrate backtracking:** every non-native substrate of
  a cascade reaction is itself produced by another cascade reaction.
- The cascade of a non-native target contains at least one reaction
  producing it.
- **Monotonicity:** `Cascade(M) ⊆ Cascade(T)` when M is an intermediate of
  T.
- Cascade of a shell-0 target is empty; cascade of an unreachable target
  raises `ValueError`.
- **Cross-algorithm:** every reaction in every enumerated pathway must
  appear in the cascade.

### 2.3 Edge cases (`tests/test_edge_cases.py`)

A second tier of inline-fixture tests that target structural failure
modes the small cell does not cover:

- Self-referential reaction `R: A → A` (matches a real MetaCyc rxn 88
  pattern) does not loop the BFS.
- A reaction producing a shell-0 chemical does not overwrite that
  chemical's native shell.
- Multi-product reactions — both products reachable.
- Deep chain (20 shells) terminates cleanly.
- Zero-substrate reactions enable at shell 1 (pins the current
  `all([])==True` behavior so any future change is conscious).
- Zero-product reactions are no-ops.
- Determinism — same inputs, same shell map every time.
- Scale smoke test — a 100-chemical synthetic corpus completes.
- Chemical-level cycle in a cascade (B ↔ C) — `enumerate_pathways` must
  not loop.
- `max_pathways` cap is enforced even on deep cascades.
- Model invariant: `Chemical` equality is by id only.
- Parser silently drops unknown substrate ids (documented risk).

### 2.4 MetaCyc smoke test (`tests/test_metacyc_smoke.py`)

A minimal regression check that runs the full pipeline on the real
MetaCyc data shipped in `data/` and asserts a small set of stable facts —
total reachable counts, a few well-known shell numbers, the lecture
example "L-tyrosine → 4-fumarylacetoacetate", etc. If anyone breaks the
parser or the BFS in a way that escapes the synthetic tests, this
catches it.

### 2.5 Real-world cross-reference: tyrosine degradation
(`validation/tyrosine_degradation_validation.md`)

A manual validation against MetaCyc pathway TYRFUMCAT-PWY ("L-tyrosine
degradation I"). We compared our shell numbers to the canonical
five-step biological pathway and got:

| Chemical | Our shell | MetaCyc step |
|---|---|---|
| L-tyrosine | 0 | (start) |
| 4-hydroxyphenylpyruvate | 1 | 1 |
| Homogentisic acid | 2 | 2 |
| 4-fumarylacetoacetate | 3 | 4 |

Two intermediate steps come back differently from canonical MetaCyc, and
both are explained by **corpus quality, not algorithm bugs**:

1. Our BFS finds an alternative (non-canonical) reaction for the
   tyrosine → 4-HPP step that does not require the 2-ketoglutarate
   cofactor (which sits at shell 1 in our model, not shell 0).
2. The intermediate 4-maleylacetoacetate is missing from
   `good_chems.txt`, which collapses MetaCyc's two-step
   dioxygenase + isomerase into a single reaction in our corpus.

The conclusion in the validation memo is explicit: the BFS is computing
the correct *minimum* given the available reactions; the discrepancies
are entirely upstream of the algorithm. This is a real example of using
domain ground truth to **localize a discrepancy** to corpus rather than
code.

### 2.6 MCP-layer tests

- `tests/test_mcp_server.py` — exercises every MCP tool / resource
  end-to-end against the real loaded hypergraph (not mocked), including
  reference resolution edge cases, structure search, fed-chemical
  resynthesis, and the interactive viewer return shape.
- `tests/test_html_render.py` — checks the generated Cytoscape JSON for
  pathways and cascades (node / edge counts, theme metadata, cofactor
  duplication mode).

### 2.7 Running everything

```bash
.venv/bin/pytest tests/ -q
```

Expected: all green, with exactly one **XFAIL** (the
`pathways.py:78–80` multi-substrate enumeration bug, kept visible on
purpose).

---

## 3. Cascade explosion: tuning vs. seeding
*Author: Jonathan Arnold*

I analyze two approaches to managing cofactor cascade snowballs, using 
a shell cutoff/max producer limit and using enzymemap data to create a 
shell zero built from ecoli native reactions. I implemented both as 
clean alternatives on their own branches.

- **shell-cutoff** (`shell-cutoff` branch) — keep the MetaCyc corpus;
  prune at every chemical during traceback. Two knobs in
  `synthesis_helper/traceback.py`: `shell_cutoff` (how much deeper than
  the chemical its substrates may lie) and
  `max_producers_per_chemical` (cap on producer reactions chosen).
- **enzymemap-shell0** (`enzymemap-shell0` branch) — replace MetaCyc
  with EnzymeMap; redefine shell 0 by running a first-pass BFS using
  only *E. coli*-annotated reactions, then run a second-pass full
  expansion seeded from that result. Cofactors fall into shell 0 by
  construction, so traceback never recurses on them. The change spans
  `synthesize.py`, `parser.py`, and the corpus loader.
- **cofactor-suppression** (`claude/sleepy-hermann-aaa57e` branch) — keep the MetaCyc corpus; add composable `SubstrateFilter` callbacks to `traceback()` that suppress cofactor substrates at the reaction level. See §4.

### 3.1 How the comparison was set up

The two approaches use different corpora, so direct cascade-size
comparison required a target set that exists in both.

1. `scripts/sample_reachables.py` (shell-cutoff branch) — sample 100
   non-shell-0 reachables uniformly from `good_chems.txt`; write
   `data/sampled_reachables.tsv` with stripped InChI (the same `/p /t
   /m /s` stripping `parser.py` already uses for loose matching).
2. `scripts/intersect_sampled.py` (enzymemap-shell0 branch) — build the
   two-pass EnzymeMap hypergraph, look up each sampled chemical by
   stripped InChI, drop any that land in shell 0 there (would be a
   no-op cascade), write `data/comparison_targets.tsv`. Of 100 sampled
   chemicals, 16 survived: matched by structure to a non-shell-0
   reachable in the EnzymeMap graph. Target shells span 1–5, plus 9
   and 13.
3. `scripts/evaluate.py --use-matches` (each branch) — re-run traceback
   on the 16 common targets and record per-target diagnostics: cascade
   reaction count, "necessary" chemicals (target plus every substrate
   of any cascade reaction), shell histogram of those, top non-shell-0
   chemicals ranked by producer reactions in the cascade, and a
   substrate-branching histogram. JSON results land in `eval_results/`.

### 3.2 Headline numbers

Across the same 16 targets, sweeping `shell_cutoff` over `{-1, 0,
none}` and `max_producers_per_chemical` over `{5, 10, 15}` on
shell-cutoff:

| approach | shell_cutoff | cap | mean rxns | median | max | mean deep chems |
|---|---|---|---|---|---|---|
| shell-cutoff | -1 | 5 | 11 | 9 | 35 | 5 |
| shell-cutoff | -1 | 10 | 16 | 12 | 63 | 5 |
| shell-cutoff | -1 | 15 | 19 | 13 | 79 | 5 |
| shell-cutoff | 0 | 5 | 56 | 17 | 257 | 15 |
| shell-cutoff | 0 | 10 | 156 | 40 | 747 | 28 |
| shell-cutoff | 0 | 15 | 245 | 66 | 1031 | 38 |
| shell-cutoff | none | 5 | 651 | 62 | 1492 | 193 |
| shell-cutoff | none | 10 | 3194 | 4252 | 4258 | 815 |
| shell-cutoff | none | 15 | 4003 | 5332 | 5337 | 952 |
| enzymemap-shell0 | n/a | n/a | 649 | 851 | 952 | 227 |

Two structural observations jump out before any per-target analysis:
`shell_cutoff` is a far stronger knob than `max_producers_per_chemical`
(cap=5 mean grows 11 → 56 → 651 across the three cutoffs), and the
shell-cutoff cap=5 row at unlimited cutoff has a mean comparable to
enzymemap-shell0 but a wildly different distribution.

### 3.3 What we found

**`shell_cutoff` is the dominant knob; the producer cap is secondary.**
At `shell_cutoff=-1` (every substrate must lie *strictly* closer to
shell 0 than the chemical itself) all 16 cascades stay under 80
reactions even at cap=15 — the cap barely matters because topology has
already done the pruning. Loosening to `shell_cutoff=0` (substrates
may share the chemical's shell) lifts the median into the tens but
keeps the cap=5 max at 257. Only with no shell cutoff does the
cofactor blowup appear. The mechanism: cofactors typically sit at
shell 1, and reactions that would walk traceback *through* a cofactor
require admitting a producer whose substrates are no closer to shell 0
than the chemical we're producing — exactly what `shell_cutoff=-1`
forbids. The cliff and the cap interact only when the cutoff is loose
enough to let traceback cross into cofactor-rich metabolism in the
first place.

**At unlimited cutoff, shell-cutoff is bimodal at cap=5.** Of 16
cascades, 9 are <80 reactions and 7 are 1426–1492. Target shell does
not predict which cluster: a shell-2 target can be 3 reactions
(`phosphoguanidinoacetate`) or 1440
(`(r)-1-aminopropan-2-yl phosphate`). The actual predictor is whether
the substrates of the chosen producers transitively reach ATP-derived
cofactors — ADP (1146 candidate producers globally), pyrophosphate
(1185), SAH (336), AMP (446). Once one of those enters the necessary
set, *its* 5 chosen producers branch into more cofactor-adjacent
reactions and the cascade snowballs. Two targets in the small cluster
(`phosphoguanidinoacetate`, `6-O-β-D-galactopyranosyl-D-glucopyranose`)
stay trivial across all caps because their substrate trees genuinely
never touch the cofactor cycle.

**Raising the cap at unlimited cutoff doesn't relax the cliff — it
merges everything into one giant cascade.** At cap=10, twelve of
sixteen cascades are within ±5 reactions of 4252. cap=15 lifts the
plateau to ~5332. Above the cliff, the cascade is just "everything
reachable via cofactor-rich metabolism"; the cap is a useful pruning
knob below the cliff and a degenerate one above it. There is no
single cap that handles both target classes at unlimited cutoff —
cap=5 still gives 1426 reactions for `pentose` while cap=10 already
collapses `dl-glyceraldehyde` (79 → 4252).

**enzymemap-shell0 collapses the variance the other way.** Without any
cap or cutoff: median 851 reactions, max 952. The E. coli first-pass
places ATP/ADP/SAM/SAH/CoA-thioester chemistry into shell 0, so
traceback never visits them. The top non-shell-0 chemicals in a
typical cascade become *sugars* — G3P, DHAP, fructose-6P, glucose-6P,
glycerol — instead of *cofactors*. Same monomorphic flavor as
shell-cutoff cap=10 but at one-fifth the size and from a fundamentally
different graph: cofactor recursion is gone because the cofactors
don't exist as recursion targets, not because we capped how often we
visit them.

**A subtle metric bug, kept as a finding.** The original `evaluate.py`
reported `max_in_count` per cascade, which we initially read as a
cap violation: "ADP appears as a product of 127 cascade reactions, but
we capped at 5." It wasn't a violation. Those 127 reactions weren't
*chosen to produce* ADP; they were chosen to produce other
intermediates and happened to emit ADP as a byproduct (any reaction
consuming ATP produces ADP). The cap is per-chemical-as-traceback-
target; the metric was per-chemical-as-product. We replaced
`cascade_shape` with a richer `target_diagnostics` that distinguishes
the two and reports both the in-cascade producer count and the global
candidate count for each necessary chemical (`scripts/evaluate.py` on
both branches; the same diagnostic lives in `scripts/debug_blowup.py`
on shell-cutoff for one-off investigations).

### 3.4 A representative target: epinephrine

Both branches have epinephrine at target shell 4. The contrast in
cascade *content* (not just size) is stark.

Cascade shape:

| metric | shell-cutoff (cutoff=0, cap=15) | enzymemap-shell0 |
|---|---|---|
| cascade reactions | 65 | 853 |
| necessary chemicals | 47 | 376 |
| non-shell-0 necessary | 13 | 298 |
| reactions per non-shell-0 chem | 5.00 | 2.86 |
| necessary chems by shell | 0:34, 1:8, 2:3, 3:1, 4:1 | 0:78, 1:33, 2:55, 3:44, 4:39, 5:33, 6:25, 7:15, 8:13, 9:11, 10:15, 11:11, 12:4 |
| reactions by # non-shell-0 substrates | 0:37, 1:25, 2:3 | 0:56, 1:742, 2:55 |

Top non-shell-0 necessary chemicals, ranked by producer reactions in
the cascade. `available` = candidate producers globally:

**shell-cutoff (cutoff=0, cap=15)**

| in_cascade | available | shell | name (id) |
|---|---|---|---|
| 15 | 229 | 1 | h2o2 (id=121) |
| 12 | 20 | 1 | Reduced Glutathione (id=1770) |
| 9 | 11 | 1 | l-dopa (id=8810) |
| 7 | 7 | 2 | dopamine (id=26474) |
| 5 | 8 | 1 | L-dehydro-ascorbate (id=5057) |
| 4 | 5 | 1 | (id=37924) |
| 3 | 3 | 1 | tyramine (id=3478) |
| 3 | 3 | 3 | norepinephrine (id=30704) |
| 3 | 3 | 4 | epinephrine (id=157085) |
| 2 | 2 | 1 | dl-dopa (id=97751) |

**enzymemap-shell0**

| in_cascade | available | shell | name (id) |
|---|---|---|---|
| 22 | 22 | 2 | D-Glyceraldehyde 3-phosphate (id=4265) |
| 20 | 20 | 3 | Glycerone phosphate (id=4267) |
| 18 | 18 | 6 | OC[C@H]1OC(O)(CO)[C@@H](O)[C@@H]1O (id=5683) |
| 16 | 16 | 4 | glycerol (id=4297) |
| 15 | 15 | 3 | D-fructose 6-phosphate (id=5783) |
| 13 | 13 | 1 | beta-D-Glucose 1-phosphate (id=5768) |
| 13 | 13 | 5 | D-fructose 6-phosphate (id=5790) |
| 12 | 12 | 4 | D-fructose {r} (id=5678) |
| 12 | 12 | 2 | D-Glucose 6-phosphate (id=5772) |
| 11 | 11 | 1 | O=C(O)[C@H](O)CO (id=4221) |

The shell-cutoff cascade is recognizably the catecholamine
biosynthesis pathway — l-dopa, dopamine, norepinephrine, tyramine —
plus oxidative side-chemistry (h2o2, glutathione, dehydroascorbate)
consistent with monoamine-oxidase / catechol-oxidation steps. The
enzymemap-shell0 cascade for the *same* molecule is dominated by
central carbon metabolism — G3P, DHAP, fructose-6P, glucose-6P,
glycerol — chemicals with no direct biosynthetic relationship to
epinephrine. Most of the 853 reactions are about how the cell builds
sugars, not how it builds adrenaline.

### 3.5 Trade-offs

**shell-cutoff** preserves the original corpus and produces genuinely
small cascades for non-cofactor targets — but only with a tight cap,
and the same cap is too tight for many easy targets and too loose for
the cofactor-touching ones. The cap is a single global knob fighting
the structure of metabolism.

**enzymemap-shell0** has a higher floor (no target costs less than
~850 reactions in this set, even ones that are biologically one step
away from a native metabolite) but a far lower ceiling. The variance
compresses dramatically and the pruning knobs become unnecessary. This
does seem to solve the issue of cofactor related snowballs, but
snowballs seem to occur anyway, pointing to deeper problems with traceback.

Enzyme map wins between these two approaches. It successfully solves the 
cofactor issue and produces computable cascades.

### 3.6 Reproducing

```bash
git checkout shell-cutoff
python scripts/sample_reachables.py            # data/sampled_reachables.tsv
git checkout enzymemap-shell0
python scripts/intersect_sampled.py            # data/comparison_targets.tsv
python scripts/evaluate.py --use-matches       # eval_results/enzymemap-shell0_matches.json
git checkout shell-cutoff
python scripts/evaluate.py --use-matches       # eval_results/shell-cutoff_matches.json (no shell cutoff)
python scripts/evaluate.py --use-matches --max-cutoff 0   # adds shell_cutoff -1 and 0 to the sweep
```

### 3.7 Caveats

The intersection (16 of 100 sampled, after dropping enzymemap shell-0
hits) is small enough that secondary findings should be treated
qualitatively. We match on stereo-stripped InChI, which can collapse
chirally distinct molecules — fine for hypergraph membership but not
safe if a downstream analysis cares about stereo. The shell-cutoff
cap≥10 plateau equals "most of the reachable hypergraph minus shell 0"
once a cofactor is admitted; the numbers are real but should not be
read as the inherent complexity of the individual targets.

### 3.8 Contributions

In addition to the cascade-explosion analysis above, my contributions
to this project include the algorithmic core that section 1 wraps as
MCP tools. Three modules on `main` form the pipeline from a reaction
corpus to fully enumerated pathways: BFS waveform expansion
(`synthesis_helper/synthesize.py`) builds a hypergraph keyed by
minimum shell distance from native metabolites; cascade traceback
(`synthesis_helper/traceback.py`) walks backward from a reachable
target and collects every reaction that transitively produces it; and
pathway enumeration (`synthesis_helper/pathways.py`) flattens the
resulting cascade into individual linear routes via choice-function
search over the producer tree.

---

## 4. Cofactor suppression in cascades
*Author: Cael Magner*

Cascades returned by `traceback()` include every reaction transitively
required to produce the target. That is correct for the algorithm, but
visually confusing: reactions that consume NADPH, FAD, or phosphate show
those cofactors as substrates that also "need to be produced," even though
the cell already has them. This section documents three proposed approaches
to suppressing cofactors in cascades, what we implemented, and what the
data showed.

### 4.1 The proposals

Three distinct strategies came out of design discussion:

**Approach A — E. coli metabolism as shell 0.**
Expand the universal metabolite set by running a BFS on known E. coli
reactions (e.g. from [EnzymeMap](https://github.com/hesther/enzymemap))
before the main expansion, so everything the cell already makes natively
starts at shell 0.
`bootstrap_ecoli_shell0.py` implements this. The blocking problem is data
equivalency: EnzymeMap and MetaCyc use different InChI forms for the same
molecules (different protonation states, oxidation states, stereo encodings),
so many chemicals that should match do not.

**Approach B — Shell-number heuristic.**
Treat any substrate at shell ≤ k as background and stop recursing there.
The intuition is that real chemists never synthesize a pathway and also plan
to re-synthesize every electron carrier from glucose. A shell threshold
captures this: cofactors are almost always reachable in very few steps
(they are highly connected), while true precursors track closer to the
target's own shell.

**Approach C — Atom-mapping (Evodex / SMARTS).**
Tag the heavy atoms in each reaction's SMARTS. A substrate whose atoms do
not end up in the product skeleton (they are merely transferred to a carrier)
is a cofactor by definition and can be dropped. This is the most principled
approach but requires atom-mapped reaction data (e.g. from
[Evodex](https://github.com/hesther/evodex)) and RDKit reaction processing —
neither of which is currently in the corpus. See §4.5.

### 4.2 Implementation: composable `SubstrateFilter`

All three approaches answer the same question at the reaction level:
*given a reaction and the hypergraph, which of its substrates should
`traceback()` recurse into?*

That maps directly to a single composable type in
`synthesis_helper/filters.py`:

```python
SubstrateFilter = Callable[[Reaction, HyperGraph], frozenset[Chemical]]
```

Five concrete filters are provided:

| Filter | What it suppresses |
|---|---|
| `shell_zero_filter` | Substrates already in shell 0 (baseline — existing behaviour). |
| `shell_threshold_filter(k)` | Substrates at shell ≤ k. |
| `inchi_normalized_filter(universal)` | Substrates whose InChI matches a universal metabolite after stripping the `/p` (proton) and `/q` (charge) layers. |
| `rhea_atom_filter(carrier_map)` | Substrates identified as carriers (cofactors, electron donors, group donors) by the Rhea reaction chemistry pipeline. Falls back to `shell_zero_filter` for reactions not in the map. |
| `compose(*filters)` | Intersection of any number of filters — a substrate must survive all of them. |

`traceback()` accepts an optional `substrate_filter` argument; passing
`None` (the default) restores the original `shell_zero_filter` behaviour,
so no existing code breaks.

```python
from synthesis_helper import traceback, compose, shell_threshold_filter, inchi_normalized_filter

# Baseline
cascade_default = traceback(hg, target)

# Suppress substrates at shell ≤ 1
cascade_thresh = traceback(hg, target, substrate_filter=shell_threshold_filter(1))

# Suppress substrates that fuzzy-match a known cofactor
cascade_inchi  = traceback(hg, target, substrate_filter=inchi_normalized_filter(universal))

# Both together
cascade_both   = traceback(hg, target,
                           substrate_filter=compose(inchi_normalized_filter(universal),
                                                    shell_threshold_filter(1)))
```

### 4.3 Experiment A — cascade-size comparison across filters

`scripts/compare_filters.py` samples target chemicals from shells 2–8 and
prints a TSV of cascade sizes and drop counts for all five filter variants.

```bash
uv run python scripts/compare_filters.py --data-dir data --per-shell 5
```

Key findings from a shell-2–5 sample (3 targets per shell,
`max_producers=10`):

| Filter | Typical cascade size | Notes |
|---|---|---|
| `baseline` (`shell_zero`) | ~4 000–4 300 reactions | Nearly every target traces back through the full hypergraph. |
| `inchi_normalized` | ~4 235–4 260 reactions | Consistently drops **~16 reactions** across all targets. |
| `thresh1` (k = 1) | **2–2 000 reactions** | Extremely sensitive to target shell; drops cleanly for low-shell targets, less so for mid-shell. |
| `thresh2` (k = 2) | 1–275 reactions | Aggressive — clips real intermediates for targets above shell 3. |
| `composed` (inchi + thresh1) | Matches `thresh1` | `thresh1` dominates; `inchi_norm` adds no further reduction once thresh1 is applied. |

Two take-aways:

1. **`inchi_normalized` is a conservative, targeted fix.** The consistent
   ~16-reaction drop corresponds exactly to the cofactors recovered by the
   InChI audit (§4.4). It does not over-prune.
2. **`shell_threshold` is aggressive and needs calibration.** k = 1
   collapses most cascades to a handful of reactions; k = 2 starts removing
   genuine intermediates for mid-shell targets. The jump between k = 1 and
   k = 2 is non-linear, which suggests the threshold is sensitive enough to
   use only as an *optional display parameter* the user can tune — not as a
   default.

### 4.4 Experiment B — InChI mismatch audit

`scripts/inchi_mismatch_audit.py` tries four progressively looser
normalization strategies for each entry in `ubiquitous_metabolites.txt`
against the 9 361 chemicals in `good_chems.txt`, and reports which entries
are missed at each level.

```bash
uv run python scripts/inchi_mismatch_audit.py --data-dir data
uv run python scripts/inchi_mismatch_audit.py --data-dir data --verbose  # per-entry detail
```

Results on the 89-entry universal list:

| Strategy | New matches | Running total | Coverage |
|---|---|---|---|
| Exact InChI | 73 | 73 | 82.0% |
| + strip `/p` (proton layer) | +9 | 82 | 92.1% |
| + strip `/q` (charge layer) | +3 | 85 | 95.5% |
| + strip stereo (`/t /m /s`) | +3 | 88 | 98.9% |
| Name fallback | 0 | 88 | 98.9% |
| **Completely unmatched** | **1** | — | — |

The 9 chemicals recovered by stripping `/p` are the most consequential:
**FAD, FADH2, NAD+, NADP+, carbonate, phosphate, sulfate, TPP**, and a
second FAD form — exactly the cofactors most likely to pollute cascade
displays. All of them appear in MetaCyc reactions with a different
protonation layer than the entry in `ubiquitous_metabolites.txt`.

The charge-layer step (+3) recovers metal ions (Mo, Ni in different
oxidation states). The stereo step (+3) recovers SAM (two stereo encodings)
and tetrahydrofolate. The one completely unmatched entry is **heme**, which
has a blank InChI field in `ubiquitous_metabolites.txt` — no normalization
strategy can match an absent string.

This confirms that `inchi_normalized_filter` (which strips `/p` and `/q`)
covers the highest-value equivalency mismatches with no ambiguity risk. The
stereo layer is worth stripping too, but with the caveat that it produces
multi-match results (up to 5 chemicals per entry for SAM), which means the
filter marks more substrates as background than the strict version would.

### 4.5 Experiment C — Rhea carrier-map filter

`scripts/build_rhea_carrier_map.py` downloads two Rhea files from the ExPASy
FTP (`rhea2ec.tsv` and `rhea-reaction-smiles.tsv`), matches MetaCyc reactions
to Rhea reactions by EC number + substrate InChI Jaccard overlap (threshold
≥ 0.3), then labels substrates as carriers using two rules:

1. **Small-molecule rule.** Substrates with ≤ 4 heavy atoms (H₂O, H⁺, O₂,
   CO₂, Pᵢ, PPᵢ, …) are always carriers regardless of context.
2. **Tanimoto-pairing rule.** A substrate whose Morgan fingerprint (radius 2,
   2048 bits) has Tanimoto similarity ≥ 0.7 with *any* product of the same
   reaction is a carrier — this catches cofactor cycling pairs (NAD⁺→NADH,
   FAD→FADH₂, ATP→ADP) whose atoms are returned to the pool rather than
   incorporated into the product skeleton.

```bash
uv run python scripts/build_rhea_carrier_map.py --data-dir data
```

Output: `data/rhea_carrier_map.tsv` — 7 933 MetaCyc reactions matched;
6 026 of those (76%) have at least one carrier identified.

`scripts/rhea_filter_experiment.py` then measures the impact on cascade sizes
across a stratified sample of targets, comparing `rhea_atom_filter`,
`inchi_normalized_filter`, and their composition against the shell-zero
baseline.

```bash
uv run python scripts/rhea_filter_experiment.py --data-dir data --per-shell 8 --max-shell 10
```

Results on 72 targets (shells 2–10, `max_producers=10`):

| Filter | Total cascade reactions | Reduction |
|---|---|---|
| `baseline` (`shell_zero`) | 247 378 | — |
| `rhea_atom_filter` | 244 923 | −2 455 (1.0%) |
| `inchi_normalized_filter` | 246 450 | −928 (0.4%) |
| `rhea_atom_filter` + `inchi_normalized_filter` | 243 995 | −3 383 (1.4%) |

Three findings worth calling out:

1. **The two filters are fully complementary.** For every large cascade, Rhea
   drops exactly 43 reactions and `inchi_norm` drops exactly 16, with zero
   overlap (43 + 16 = 59 every time). They are finding categorically different
   cofactors: `inchi_norm` catches molecules whose protonation or charge state
   differs from the universal list, while Rhea catches molecules that weren't
   in that list at all but are structurally paired with a product.

2. **The modest absolute reduction is expected.** The Rhea map covers **40.6%
   of all reactions visited** in a typical baseline cascade, so the filter has
   plenty of opportunity to act. The 1% reduction is small because most major
   cofactors (NAD⁺, FAD, CoA, ATP) are already in shell 0 via
   `ubiquitous_metabolites.txt` and need no further filtering. The 43
   reactions that Rhea prunes are the subset of cascades where a carrier has
   shell > 0 — real cofactors the universal list missed.

3. **Rhea contributes 2.7× more than InChI normalization** in this category.
   Together they constitute a complementary defense: normalization corrects
   data-source encoding drift; Rhea detects carriers the curated list never
   enumerated.

### 4.6 Future work: Evodex / true atom-mapping

`rhea_atom_filter` approximates atom tracking with Tanimoto similarity rather
than actually following atoms through the reaction. This works well for
cofactor cycling pairs (high self-similarity) but can miss group-transfer
carriers where the donor and acceptor forms are structurally dissimilar. The
fully principled approach would use atom-mapped reaction SMARTS — available
from [Evodex](https://github.com/hesther/evodex) — where the carrier
classification is exact: any substrate atom that does not appear in the
product skeleton is definitionally transferred.

Upgrading from Rhea to Evodex would require:

1. **Atom-mapped reaction data.** The current `Reaction` model carries only
   substrate/product sets and an EC number. Atom-mapping requires SMARTS
   reaction strings with explicit atom tags — available from Evodex or from
   the mapped reaction field in some MetaCyc exports, but not currently
   loaded by the parser.
2. **RDKit reaction processing.** Parsing a reaction SMARTS, applying it,
   and reading which input atoms ended up in which output atoms requires
   `rdkit.Chem.AllChem.ReactionFromSmarts` — already a dependency in
   `scripts/build_rhea_carrier_map.py`, so the package constraint is met,
   but the parser and `Reaction` model would need extension.
3. **A new `SubstrateFilter` implementation.** The filter infrastructure is
   already in place and the type signature is compatible. Once the atom-map
   data is available, an Evodex filter slips in as one more callable without
   touching `traceback()` or `pathways.py`.

---

## Run the standalone pipeline

If you just want the reachables file without going through MCP:

```bash
.venv/bin/python main.py
```

Writes `data/metacyc_L2_reachables.txt` (columns: id, name, inchi, shell).
