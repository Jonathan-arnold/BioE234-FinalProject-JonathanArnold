"""Build data/rhea_carrier_map.tsv from Rhea reaction SMILES.

For each reaction in good_reactions.txt that can be matched to a Rhea
reaction, identifies which substrates are carriers (cofactors / electron
donors / group donors) rather than skeleton substrates.

Carrier detection — a substrate is a carrier if ANY of:
  1. It has a product on the same reaction with Tanimoto similarity ≥ threshold
     (cofactor cycling: NAD+→NADH, FAD→FADH2, ATP→ADP, etc.)
  2. It has ≤ MAX_SMALL_ATOMS heavy atoms (H2O, H+, O2, CO2, Pi, PPi, …)

Matching strategy:
  - Our reactions have a single EC number (or a multi-EC list).
  - Rhea reactions have EC numbers in rhea2ec.tsv.
  - For each (our_rxn, EC) → candidate Rhea reactions with the same EC.
  - Best candidate = highest Jaccard overlap of substrate InChI sets.
  - Match accepted if Jaccard ≥ MIN_JACCARD.

Output format (rhea_carrier_map.tsv):
  rxnid<TAB>norm_inchi1,norm_inchi2,...

Usage:
    uv run python scripts/build_rhea_carrier_map.py
    uv run python scripts/build_rhea_carrier_map.py --data-dir /path/to/data
    uv run python scripts/build_rhea_carrier_map.py --tanimoto 0.8 --verbose
"""

from __future__ import annotations

import argparse
import io
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, DataStructs
    from rdkit.Chem.inchi import MolToInchi
    RDLogger.DisableLog("rdApp.*")
except ImportError:
    sys.exit("rdkit is required. Install with: pip install rdkit")

from synthesis_helper.parser import parse_chemicals, parse_reactions, _read_lines


DEFAULT_DATA_DIR = Path(__file__).resolve().parents[1] / "data"

RHEA_EC_URL    = "https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv"
RHEA_SMILES_URL = "https://ftp.expasy.org/databases/rhea/tsv/rhea-reaction-smiles.tsv"

MAX_SMALL_ATOMS = 4   # substrates with ≤ this many heavy atoms are always carriers
MIN_JACCARD     = 0.3  # minimum substrate InChI overlap to accept a Rhea match


# ---------------------------------------------------------------------------
# InChI helpers
# ---------------------------------------------------------------------------

_smiles_inchi_cache: dict[str, str | None] = {}


def smiles_to_inchi(smiles: str) -> str | None:
    if smiles in _smiles_inchi_cache:
        return _smiles_inchi_cache[smiles]
    mol = Chem.MolFromSmiles(smiles)
    result = MolToInchi(mol) if mol is not None else None
    _smiles_inchi_cache[smiles] = result
    return result


def normalize_inchi(inchi: str) -> str:
    s = re.sub(r"/p[+-]\d+", "", inchi.strip('"'))
    return re.sub(r"/q[+-]\d+", "", s)


# ---------------------------------------------------------------------------
# Rhea download helpers
# ---------------------------------------------------------------------------

def _fetch(url: str, label: str) -> str:
    print(f"  Fetching {label} ...", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=120) as resp:
        return resp.read().decode("utf-8")


def load_rhea_ec(url: str) -> dict[str, list[int]]:
    """EC → list of Rhea master IDs (UN direction only)."""
    text = _fetch(url, "rhea2ec.tsv")
    ec_to_rhea: dict[str, list[int]] = defaultdict(list)
    for line in text.splitlines():
        if not line or line.startswith("RHEA_ID"):
            continue
        parts = line.split("\t")
        if len(parts) < 4 or parts[1] != "UN":
            continue
        try:
            master_id = int(parts[2])
        except ValueError:
            continue
        ec = parts[3].strip()
        if ec:
            ec_to_rhea[ec].append(master_id)
    print(f"    {len(ec_to_rhea)} ECs, {sum(len(v) for v in ec_to_rhea.values())} Rhea reactions",
          file=sys.stderr)
    return dict(ec_to_rhea)


def load_rhea_smiles(url: str) -> dict[int, str]:
    """Rhea directional ID → reaction SMILES.

    rhea-reaction-smiles.tsv has two columns: RHEA_ID<TAB>REACTION_SMILES.
    Rhea's numbering: master_id = UN, master_id+1 = LR, master_id+2 = RL,
    master_id+3 = BI. We load all IDs; callers look up master_id+1 for the
    canonical left→right (substrates >> products) SMILES.
    """
    text = _fetch(url, "rhea-reaction-smiles.tsv")
    id_to_smiles: dict[int, str] = {}
    for line in text.splitlines():
        if not line or line.lower().startswith("rhea"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        try:
            rhea_id = int(parts[0])
        except ValueError:
            continue
        smiles = parts[1].strip()
        if smiles and ">>" in smiles:
            id_to_smiles[rhea_id] = smiles
    print(f"    {len(id_to_smiles)} reaction SMILES loaded", file=sys.stderr)
    return id_to_smiles


# ---------------------------------------------------------------------------
# EC number parsing
# ---------------------------------------------------------------------------

_EC_RE = re.compile(r"\d+\.\d+\.\d+\.(?:\d+|[a-zA-Z-]+)")


def extract_ecs(raw: str) -> list[str]:
    """Extract valid EC numbers from a possibly list-formatted string."""
    return _EC_RE.findall(raw)


def is_complete_ec(ec: str) -> bool:
    """Return True only for fully specified ECs (no wildcards like 1.7.1.-)."""
    return bool(re.fullmatch(r"\d+\.\d+\.\d+\.\d+", ec))


# ---------------------------------------------------------------------------
# Reaction SMILES parsing
# ---------------------------------------------------------------------------

def split_rxn_smiles(rxn_smiles: str) -> tuple[list[str], list[str]]:
    left, right = rxn_smiles.split(">>", 1)
    return [s for s in left.split(".") if s], [s for s in right.split(".") if s]


def smiles_list_to_inchis(smiles_list: list[str]) -> list[str | None]:
    return [smiles_to_inchi(s) for s in smiles_list]


# ---------------------------------------------------------------------------
# Carrier detection
# ---------------------------------------------------------------------------

def fp_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)


def detect_carrier_inchis(
    sub_smiles: list[str],
    prod_smiles: list[str],
    tanimoto_threshold: float,
) -> frozenset[str]:
    """Return normalized InChIs of substrates identified as carriers."""
    carriers: set[str] = set()

    sub_fps = [(smi, fp_from_smiles(smi)) for smi in sub_smiles]
    prod_fps = [(smi, fp_from_smiles(smi)) for smi in prod_smiles]

    for sub_smi, sub_fp in sub_fps:
        sub_mol = Chem.MolFromSmiles(sub_smi)
        if sub_mol is None:
            continue

        # Rule 1: small molecule → always a carrier
        if sub_mol.GetNumHeavyAtoms() <= MAX_SMALL_ATOMS:
            inchi = smiles_to_inchi(sub_smi)
            if inchi:
                carriers.add(normalize_inchi(inchi))
            continue

        # Rule 2: structural cofactor pairing with any product
        if sub_fp is not None:
            for _, prod_fp in prod_fps:
                if prod_fp is None:
                    continue
                sim = DataStructs.TanimotoSimilarity(sub_fp, prod_fp)
                if sim >= tanimoto_threshold:
                    inchi = smiles_to_inchi(sub_smi)
                    if inchi:
                        carriers.add(normalize_inchi(inchi))
                    break

    return frozenset(carriers)


# ---------------------------------------------------------------------------
# Reaction matching
# ---------------------------------------------------------------------------

def reaction_sub_inchis(our_rxn, chemicals_by_id: dict) -> frozenset[str]:
    result = set()
    for chem in our_rxn.substrates:
        raw = chemicals_by_id.get(chem.id, {}).get("inchi", "")
        if raw:
            result.add(normalize_inchi(raw.strip('"')))
    return frozenset(result)


def jaccard(a: frozenset, b: frozenset) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    if not union:
        return 0.0
    return len(a & b) / len(union)


def best_rhea_match(
    our_sub_inchis: frozenset[str],
    candidate_ids: list[int],
    rhea_smiles: dict[int, str],
) -> tuple[int | None, float]:
    """Return (best_rhea_master_id, jaccard_score) or (None, 0).

    Looks up master_id+1 (LR direction) in the SMILES table.
    """
    best_id, best_j = None, 0.0
    for rid in candidate_ids:
        rxn_smi = rhea_smiles.get(rid + 1)  # LR direction = master_id + 1
        if not rxn_smi:
            continue
        sub_smis, _ = split_rxn_smiles(rxn_smi)
        rhea_inchis = frozenset(
            normalize_inchi(i)
            for i in smiles_list_to_inchis(sub_smis)
            if i is not None
        )
        j = jaccard(our_sub_inchis, rhea_inchis)
        if j > best_j:
            best_j, best_id = j, rid
    return best_id, best_j


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=str(DEFAULT_DATA_DIR))
    ap.add_argument("--tanimoto", type=float, default=0.7,
                    help="Tanimoto threshold for cofactor-pair detection (default 0.7)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    if not (data_dir / "good_chems.txt").exists():
        sys.exit(f"Data not found at {data_dir}. Pass --data-dir or run download_data.py.")

    print("Loading MetaCyc data ...", file=sys.stderr)
    chemicals = parse_chemicals(data_dir / "good_chems.txt")
    reactions = parse_reactions(data_dir / "good_reactions.txt", chemicals)
    chem_by_id = {c.id: {"inchi": c.inchi, "smiles": c.smiles} for c in chemicals.values()}
    print(f"  {len(chemicals):,} chemicals, {len(reactions):,} reactions", file=sys.stderr)

    print("Fetching Rhea data ...", file=sys.stderr)
    ec_to_rhea = load_rhea_ec(RHEA_EC_URL)
    rhea_smiles = load_rhea_smiles(RHEA_SMILES_URL)

    print("Matching reactions and detecting carriers ...", file=sys.stderr)
    output: dict[int, frozenset[str]] = {}
    n_matched = 0
    n_skipped_no_ec = 0
    n_skipped_no_candidate = 0
    n_skipped_low_jaccard = 0

    for rxn in reactions:
        ecs = [ec for ec in extract_ecs(rxn.ecnum or "") if is_complete_ec(ec)]
        if not ecs:
            n_skipped_no_ec += 1
            continue

        our_sub_inchis = reaction_sub_inchis(rxn, chem_by_id)

        best_id, best_j = None, 0.0
        for ec in ecs:
            candidates = ec_to_rhea.get(ec, [])
            if not candidates:
                continue
            rid, j = best_rhea_match(our_sub_inchis, candidates, rhea_smiles)
            if j > best_j:
                best_j, best_id = j, rid

        if best_id is None:
            n_skipped_no_candidate += 1
            continue
        if best_j < MIN_JACCARD:
            n_skipped_low_jaccard += 1
            if args.verbose:
                print(f"    rxn {rxn.id} (EC {rxn.ecnum}): best Jaccard {best_j:.2f} < {MIN_JACCARD}, skipping",
                      file=sys.stderr)
            continue

        rxn_smi = rhea_smiles[best_id + 1]
        sub_smis, prod_smis = split_rxn_smiles(rxn_smi)
        carriers = detect_carrier_inchis(sub_smis, prod_smis, args.tanimoto)
        output[rxn.id] = carriers
        n_matched += 1

        if args.verbose and carriers:
            print(f"    rxn {rxn.id} (EC {rxn.ecnum}, Rhea {best_id}, J={best_j:.2f}): "
                  f"{len(carriers)} carrier(s)", file=sys.stderr)

    print(f"\n  Matched: {n_matched}", file=sys.stderr)
    print(f"  Skipped (no EC): {n_skipped_no_ec}", file=sys.stderr)
    print(f"  Skipped (no Rhea candidate): {n_skipped_no_candidate}", file=sys.stderr)
    print(f"  Skipped (low Jaccard): {n_skipped_low_jaccard}", file=sys.stderr)

    out_path = data_dir / "rhea_carrier_map.tsv"
    with open(out_path, "w") as fh:
        fh.write("rxnid\tcarrier_inchis\n")
        for rxn_id in sorted(output):
            inchis = ",".join(sorted(output[rxn_id]))
            fh.write(f"{rxn_id}\t{inchis}\n")
    print(f"\nWrote {len(output)} entries → {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
