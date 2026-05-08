"""Build the EnzymeMap-derived chemical and reaction corpus.

Reads ``data/enzymemap_processed_reactions.csv.gz``, keeps only rows with
``natural == True``, canonicalizes every molecule via RDKit InChI, assigns
stable integer ids (sorted by canonical InChI so ids are reproducible across
runs), and emits four parser-compatible TSVs:

  data/enzymemap_chems.tsv               id  name  inchi  smiles
  data/enzymemap_reactions.tsv           rxnid  ecnum  substrates  products
  data/enzymemap_reaction_ecs.tsv        rxnid  ecnums           (comma-sep)
  data/enzymemap_reaction_organisms.tsv  rxnid  organisms        (comma-sep)

Reactions are deduped on ``(frozenset(substrate_ids), frozenset(product_ids))``;
all EC numbers and organisms seen across collapsed rows are preserved in the
sidecars. ``enzymemap_reactions.tsv`` carries the lex-first EC so the existing
``parse_reactions`` consumer keeps working unchanged.
"""

from __future__ import annotations

import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.inchi import MolToInchi
except ImportError:
    sys.stderr.write("rdkit is required. Install with: pip install rdkit\n")
    raise

RDLogger.DisableLog("rdApp.*")

DATA_DIR = Path(__file__).parent / "data"
ENZYMEMAP_URL = (
    "https://github.com/hesther/enzymemap/raw/main/data/processed_reactions.csv.gz"
)
ENZYMEMAP_PATH = DATA_DIR / "enzymemap_processed_reactions.csv.gz"

CHEMS_OUT = DATA_DIR / "enzymemap_chems.tsv"
REACTIONS_OUT = DATA_DIR / "enzymemap_reactions.tsv"
RXN_ECS_OUT = DATA_DIR / "enzymemap_reaction_ecs.tsv"
RXN_ORGS_OUT = DATA_DIR / "enzymemap_reaction_organisms.tsv"


def ensure_csv() -> None:
    if ENZYMEMAP_PATH.exists():
        return
    print(f"Downloading EnzymeMap CSV to {ENZYMEMAP_PATH} ...")
    ENZYMEMAP_PATH.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(ENZYMEMAP_URL, ENZYMEMAP_PATH)


def split_reaction_smiles(rxn_smiles: str) -> tuple[list[str], list[str]]:
    """'A.B>>C.D' -> (['A','B'], ['C','D']). Agents (middle of A>B>C) are dropped."""
    if ">>" in rxn_smiles:
        left, right = rxn_smiles.split(">>", 1)
    else:
        parts = rxn_smiles.split(">")
        if len(parts) != 3:
            return [], []
        left, _, right = parts
    return [s for s in left.split(".") if s], [s for s in right.split(".") if s]


def split_orig_text(orig: str) -> tuple[list[str], list[str]] | None:
    """'a + b = c + d' -> (['a','b'], ['c','d']). Returns None if unparseable.

    Splits on ' = ' then ' + '. Names containing literal '+' (e.g. 'NAD+', 'H+')
    survive because the separator includes surrounding spaces.
    """
    if " = " not in orig:
        return None
    left, right = orig.split(" = ", 1)
    subs = [s.strip() for s in left.split(" + ") if s.strip()]
    prods = [s.strip() for s in right.split(" + ") if s.strip()]
    if not subs or not prods:
        return None
    return subs, prods


def smiles_to_canon(
    smiles: str,
    cache: dict[str, tuple[str, str] | None],
) -> tuple[str, str] | None:
    """Return (canonical_inchi, canonical_smiles) or None if RDKit can't parse."""
    if smiles in cache:
        return cache[smiles]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        cache[smiles] = None
        return None
    inchi = MolToInchi(mol)
    if not inchi.startswith("InChI="):
        cache[smiles] = None
        return None
    canon_smi = Chem.MolToSmiles(mol)
    cache[smiles] = (inchi, canon_smi)
    return cache[smiles]


def _safe(text: str) -> str:
    """Strip characters that would corrupt a TSV row."""
    return text.replace("\t", " ").replace("\n", " ").replace("\r", " ").strip()


def main() -> None:
    ensure_csv()

    print(f"Loading EnzymeMap from {ENZYMEMAP_PATH} ...")
    df = pd.read_csv(ENZYMEMAP_PATH, compression="gzip", low_memory=False)
    print(f"  {len(df)} rows total")

    df = df[df["natural"] == True].copy()  # noqa: E712 — pandas bool compare
    print(f"  {len(df)} rows with natural=True")

    inchi_to_smiles: dict[str, str] = {}      # canonical InChI -> canonical SMILES
    name_votes: dict[str, Counter] = defaultdict(Counter)  # inchi -> {name: count}
    smiles_cache: dict[str, tuple[str, str] | None] = {}

    # rxn key = (frozenset(substrate_inchi), frozenset(product_inchi))
    rxn_to_ecs: dict[tuple[frozenset[str], frozenset[str]], set[str]] = defaultdict(set)
    rxn_to_orgs: dict[tuple[frozenset[str], frozenset[str]], set[str]] = defaultdict(set)
    rxn_first_seen: list[tuple[frozenset[str], frozenset[str]]] = []
    seen: set[tuple[frozenset[str], frozenset[str]]] = set()

    bad_smiles = 0
    bad_text = 0
    name_count_mismatch = 0

    for _, row in df.iterrows():
        rxn_smiles = row["unmapped"]
        if not isinstance(rxn_smiles, str):
            bad_smiles += 1
            continue
        sub_smi, prod_smi = split_reaction_smiles(rxn_smiles)
        if not sub_smi or not prod_smi:
            bad_smiles += 1
            continue

        sub_inchis: list[str] = []
        prod_inchis: list[str] = []
        ok = True
        for smi_list, bucket in ((sub_smi, sub_inchis), (prod_smi, prod_inchis)):
            for smi in smi_list:
                canon = smiles_to_canon(smi, smiles_cache)
                if canon is None:
                    ok = False
                    break
                inchi, canon_smi = canon
                bucket.append(inchi)
                inchi_to_smiles.setdefault(inchi, canon_smi)
            if not ok:
                break
        if not ok or not sub_inchis or not prod_inchis:
            bad_smiles += 1
            continue

        # Harvest names from orig_rxn_text only when a side has exactly one
        # molecule — then the (smiles, name) alignment is unambiguous. Multi-
        # molecule sides are skipped because EnzymeMap's SMILES order doesn't
        # reliably match the textual order, so positional zip mis-labels
        # cofactors with each other (e.g. arsenate gets named "NAD+").
        orig = row.get("orig_rxn_text")
        if isinstance(orig, str):
            parsed = split_orig_text(orig)
            if parsed is None:
                bad_text += 1
            else:
                sub_names, prod_names = parsed
                if len(sub_inchis) == 1 and len(sub_names) == 1:
                    name_votes[sub_inchis[0]][_safe(sub_names[0])] += 1
                if len(prod_inchis) == 1 and len(prod_names) == 1:
                    name_votes[prod_inchis[0]][_safe(prod_names[0])] += 1
                if len(sub_names) != len(sub_smi) or len(prod_names) != len(prod_smi):
                    name_count_mismatch += 1

        key = (frozenset(sub_inchis), frozenset(prod_inchis))
        if key not in seen:
            seen.add(key)
            rxn_first_seen.append(key)

        ec = row.get("ec_num")
        if isinstance(ec, str) and ec:
            rxn_to_ecs[key].add(ec.strip())
        org = row.get("organism")
        if isinstance(org, str) and org:
            rxn_to_orgs[key].add(_safe(org))

    print(f"  Skipped {bad_smiles} rows with bad/unresolvable SMILES")
    print(f"  Skipped name extraction on {bad_text} unparseable + {name_count_mismatch} count-mismatched orig texts")
    print(f"  Unique chemicals: {len(inchi_to_smiles)}")
    print(f"  Unique reactions: {len(rxn_first_seen)}")

    # Stable IDs: sort chemicals by canonical InChI.
    sorted_inchis = sorted(inchi_to_smiles)
    inchi_to_id: dict[str, int] = {inchi: i + 1 for i, inchi in enumerate(sorted_inchis)}

    # Pick most-voted name per chemical; fallback to canonical SMILES.
    def chem_name(inchi: str) -> str:
        votes = name_votes.get(inchi)
        if votes:
            name, _ = votes.most_common(1)[0]
            if name:
                return name
        return inchi_to_smiles[inchi]

    print(f"Writing chemicals → {CHEMS_OUT}")
    with open(CHEMS_OUT, "w") as f:
        f.write("id\tname\tinchi\tsmiles\n")
        for inchi in sorted_inchis:
            cid = inchi_to_id[inchi]
            name = _safe(chem_name(inchi))
            smi = _safe(inchi_to_smiles[inchi])
            f.write(f"{cid}\t{name}\t{inchi}\t{smi}\n")

    # Reactions get ids in first-seen order — stable as long as the input rows are.
    rxn_ids: dict[tuple[frozenset[str], frozenset[str]], int] = {
        key: i + 1 for i, key in enumerate(rxn_first_seen)
    }

    print(f"Writing reactions → {REACTIONS_OUT}")
    with open(REACTIONS_OUT, "w") as f:
        f.write("rxnid\tecnum\tsubstrates\tproducts\n")
        for key in rxn_first_seen:
            rid = rxn_ids[key]
            ecs = sorted(rxn_to_ecs.get(key, set()))
            primary_ec = ecs[0] if ecs else ""
            sub_ids = " ".join(str(inchi_to_id[i]) for i in sorted(key[0], key=inchi_to_id.get))
            prod_ids = " ".join(str(inchi_to_id[i]) for i in sorted(key[1], key=inchi_to_id.get))
            f.write(f"{rid}\t{primary_ec}\t{sub_ids}\t{prod_ids}\n")

    print(f"Writing EC sidecar → {RXN_ECS_OUT}")
    with open(RXN_ECS_OUT, "w") as f:
        f.write("rxnid\tecnums\n")
        for key in rxn_first_seen:
            ecs = sorted(rxn_to_ecs.get(key, set()))
            if not ecs:
                continue
            f.write(f"{rxn_ids[key]}\t{','.join(ecs)}\n")

    print(f"Writing organism sidecar → {RXN_ORGS_OUT}")
    with open(RXN_ORGS_OUT, "w") as f:
        f.write("rxnid\torganisms\n")
        for key in rxn_first_seen:
            orgs = sorted(rxn_to_orgs.get(key, set()))
            if not orgs:
                continue
            # Organisms can contain commas — escape as ';' to keep the format simple.
            cleaned = [o.replace(",", ";") for o in orgs]
            f.write(f"{rxn_ids[key]}\t{','.join(cleaned)}\n")

    print("Done.")


if __name__ == "__main__":
    main()
