## synthesis_helper

A bottom-up retrobiosynthesis tool that builds a synthesis hypergraph from a
reaction corpus (MetaCyc) and enumerates **cascades, pathways, and enzyme
compositions** for any reachable target chemical. The whole thing is wrapped
as a [Model Context Protocol](https://modelcontextprotocol.io/) server so
Claude Code can drive it as an interactive lab notebook.

# Final Project Report
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
- **cofactor-suppression** (`claude/sleepy-hermann-aaa57e` branch) — this is owned and analyzed below by Cael Magner.

### 1 How the comparison was set up

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

### 2 Headline numbers

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

### 3 What we found

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

### 4 A representative target: epinephrine

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

### 5 Trade-offs

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

### 6 Reproducing

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

### 7 Caveats

The intersection (16 of 100 sampled, after dropping enzymemap shell-0
hits) is small enough that secondary findings should be treated
qualitatively. We match on stereo-stripped InChI, which can collapse
chirally distinct molecules — fine for hypergraph membership but not
safe if a downstream analysis cares about stereo. The shell-cutoff
cap≥10 plateau equals "most of the reachable hypergraph minus shell 0"
once a cofactor is admitted; the numbers are real but should not be
read as the inherent complexity of the individual targets.

### 8 Contributions

In addition to the cascade-explosion analysis above, my contributions
to this project include the algorithmic core of the code. Three modules 
on `main` form the pipeline from a reaction
corpus to fully enumerated pathways: BFS waveform expansion
(`synthesis_helper/synthesize.py`) builds a hypergraph keyed by
minimum shell distance from native metabolites; cascade traceback
(`synthesis_helper/traceback.py`) walks backward from a reachable
target and collects every reaction that transitively produces it; and
pathway enumeration (`synthesis_helper/pathways.py`) flattens the
resulting cascade into individual linear routes via choice-function
search over the producer tree.
