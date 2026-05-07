# Synthesis Helper — Cascade Explosion Experiments

This repo compares three approaches to controlling the combinatorial explosion that occurs when building cascades for a reachable. Each approach lives on its own branch as a pure implementation; `main` holds the working hybrid that integrates with the MCP wrappers.

## The three approaches

1. **shell-cutoff** — When recursing through reactions that produce a chemical, include a reaction only if all its substrates lie in shell ≤ (chemical's shell + n), for a chosen cutoff `n >= -1`. Tighter `n` prunes more aggressively (e.g., `n=-1` on a shell-3 chemical admits only reactions whose substrates are all in shell ≤ 2).

2. **enzymemap-shell0** — Replace MetaCyc with EnzyMap. Seed shell 0 with ubiquitous metabolites, run the full hypergraph expansion using only reactions annotated as natively occurring in *E. coli*, then treat every reachable from that expansion as the new shell 0 for downstream cascade work.

3. **evodex-heavy-atom** — Replace MetaCyc with EVODEX partial reactions. A reaction producing a reachable is included in its cascade only if it contributes heavy atoms to that reachable.
