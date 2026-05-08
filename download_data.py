"""Download upstream data for the EnzymeMap-shell0 pipeline.

Pulls the curated metabolite seed lists from Google Drive and the EnzymeMap
processed-reactions CSV from GitHub. After downloading, builds the EnzymeMap
TSV corpus (chems / reactions / EC + organism sidecars) by invoking
``build_enzymemap_corpus.py``.

The MetaCyc good_chems / good_reactions files are no longer pulled — the
pipeline derives its chemical and reaction inventories from EnzymeMap.
"""

from __future__ import annotations

import os
import subprocess
import sys
import urllib.request
from pathlib import Path

import gdown

DATA_DIR = Path(__file__).parent / "data"

# Curated organism-agnostic seed lists (still hand-maintained, not in EnzymeMap).
GDRIVE_FILES = {
    "minimal_metabolites.txt":    "17rmvSCeBm0ZRdfFmnOGwhLluF7rIz6Dz",
    "ubiquitous_metabolites.txt": "1n2YZcBSJI9hemkwZNrOO4WiJC4IbXB_f",
}

ENZYMEMAP_URL = (
    "https://github.com/hesther/enzymemap/raw/main/data/processed_reactions.csv.gz"
)
ENZYMEMAP_PATH = DATA_DIR / "enzymemap_processed_reactions.csv.gz"

ENZYMEMAP_TSVS = [
    "enzymemap_chems.tsv",
    "enzymemap_reactions.tsv",
    "enzymemap_reaction_ecs.tsv",
    "enzymemap_reaction_organisms.tsv",
]


def main() -> None:
    DATA_DIR.mkdir(exist_ok=True)

    for filename, file_id in GDRIVE_FILES.items():
        out = DATA_DIR / filename
        if out.exists():
            print(f"Already present: {filename}")
            continue
        print(f"Downloading {filename} ...")
        gdown.download(id=file_id, output=str(out), quiet=False)

    if ENZYMEMAP_PATH.exists():
        print(f"Already present: {ENZYMEMAP_PATH.name}")
    else:
        print(f"Downloading {ENZYMEMAP_PATH.name} ...")
        urllib.request.urlretrieve(ENZYMEMAP_URL, ENZYMEMAP_PATH)

    if all((DATA_DIR / t).exists() for t in ENZYMEMAP_TSVS):
        print("EnzymeMap TSV corpus already built.")
    else:
        print("Building EnzymeMap corpus (chems + reactions + sidecars) ...")
        subprocess.run(
            [sys.executable, str(Path(__file__).parent / "build_enzymemap_corpus.py")],
            check=True,
        )

    print("\nAvailable data files:")
    for filename in sorted(os.listdir(DATA_DIR)):
        print(f"  {DATA_DIR / filename}")


if __name__ == "__main__":
    main()
