"""Sample reachables from the shell-cutoff hypergraph for cross-branch comparison.

Builds the shell-cutoff hypergraph (MetaCyc good_chems / good_reactions seeded
with native + universal metabolites), samples N reachables uniformly from all
non-shell-0 chemicals, and writes them to data/sampled_reachables.tsv with
their stripped InChI so the enzymemap-shell0 branch can match by structure.

Usage:
    python scripts/sample_reachables.py
"""

from __future__ import annotations

import random
from pathlib import Path

from synthesis_helper.parser import (
    _strip_stereo,
    parse_chemicals,
    parse_metabolite_list,
    parse_reactions,
)
from synthesis_helper.synthesize import synthesize

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
OUT_PATH = DATA_DIR / "sampled_reachables.tsv"
SEED = 20260509
N_SAMPLES = 50


def main() -> None:
    chemicals = parse_chemicals(DATA_DIR / "good_chems.txt")
    reactions = parse_reactions(DATA_DIR / "good_reactions.txt", chemicals)
    native = parse_metabolite_list(DATA_DIR / "minimal_metabolites.txt", chemicals)
    universal = parse_metabolite_list(
        DATA_DIR / "ubiquitous_metabolites.txt", chemicals
    )
    print(f"Loaded {len(chemicals)} chemicals, {len(reactions)} reactions")
    hg = synthesize(reactions, native, universal, verbose=True)

    pool = sorted(
        (c for c, s in hg.chemical_to_shell.items() if s > 0),
        key=lambda c: c.id,
    )
    print(f"Non-shell-0 reachables: {len(pool)}")

    rng = random.Random(SEED)
    k = min(N_SAMPLES, len(pool))
    chosen = rng.sample(pool, k)

    with OUT_PATH.open("w") as f:
        f.write("id\tshell\tname\tinchi_stripped\tinchi\tsmiles\n")
        for c in chosen:
            shell = hg.chemical_to_shell[c]
            inchi_stripped = _strip_stereo(c.inchi) if c.inchi else ""
            f.write(
                f"{c.id}\t{shell}\t{c.name}\t{inchi_stripped}\t{c.inchi}\t{c.smiles}\n"
            )
    print(f"Wrote {k} samples to {OUT_PATH}")


if __name__ == "__main__":
    main()
