"""Hand-drawn 'small cell' fixture for algorithm verification.

A tiny synthetic metabolism where every shell number, cascade member, and
pathway is pre-computable on paper. See validation/small_cell.md for the
diagram and ground-truth tables.

Topology: 15 chemicals (14 active + 1 unreachable substrate), 12 reactions.
Natives (shell 0):  N1, N2
Universal (shell 0): U1
Fed (only when explicitly added to natives): F1
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from synthesis_helper.models import Chemical, Reaction


def _chem(
    id_: int,
    name: str,
    smiles: str = "",
    inchi: str | None = None,
) -> Chemical:
    """Build a test Chemical. By default uses a synthetic InChI that RDKit
    won't parse (``InChI=1S/SMALL_CELL/<name>``); pass a real ``inchi`` /
    ``smiles`` when a test needs RDKit to actually fingerprint the chem.
    """
    return Chemical(
        id=id_,
        name=name,
        inchi=inchi if inchi is not None else f"InChI=1S/SMALL_CELL/{name}",
        smiles=smiles,
    )


@dataclass
class SmallCell:
    chems: dict[str, Chemical] = field(default_factory=dict)
    rxns: dict[str, Reaction] = field(default_factory=dict)
    natives: set[Chemical] = field(default_factory=set)
    universals: set[Chemical] = field(default_factory=set)
    fed: set[Chemical] = field(default_factory=set)

    @property
    def all_reactions(self) -> list[Reaction]:
        return list(self.rxns.values())


def build_small_cell() -> SmallCell:
    # A few chems carry REAL SMILES / InChI so RDKit-based tools
    # (find_by_structure) can fingerprint them. Ethanol / acetaldehyde /
    # acetate form a structurally related triplet so Tanimoto scoring is
    # meaningful. The others keep the synthetic ``InChI=1S/SMALL_CELL/...``
    # which RDKit rejects — intentional, exercises the skip-on-parse-fail
    # path without needing a pristine 15-molecule dataset.
    N1 = _chem(101, "N1", smiles="CCO",
               inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")  # ethanol
    N2 = _chem(102, "N2")
    U1 = _chem(103, "U1")
    F1 = _chem(104, "F1")
    M1 = _chem(201, "M1", smiles="CC=O")     # acetaldehyde
    M2 = _chem(202, "M2", smiles="CC(=O)O")  # acetic acid
    M_deep = _chem(203, "M_deep")
    M4 = _chem(204, "M4")
    M5 = _chem(205, "M5")
    M3a = _chem(206, "M3a")
    M3b = _chem(207, "M3b")
    M6 = _chem(208, "M6")
    M7 = _chem(209, "M7")
    T = _chem(300, "T")
    X_missing = _chem(999, "X_missing")

    chems = {
        c.name: c
        for c in [N1, N2, U1, F1, M1, M2, M_deep, M4, M5, M3a, M3b, M6, M7, T, X_missing]
    }

    # Real ECs on three reactions so compare_pathways / pathway_to_composition
    # have signal to test against:
    #   R6 = 1.14.13.1  (cytochrome P450 monooxygenase — triggers is_p450)
    #   R7 = 2.7.1.1    (hexokinase — has entry in ec_names → NOT orphan)
    #   R9 = 1.1.1.1    (alcohol dehydrogenase — has entry in ec_names)
    # R9 is off the main T pathway (feeds M3b), so the T pathway containing
    # R7 is the "no-orphan" path and the T pathway containing R6+R5 is
    # "orphan" — ideal for sort-order tests.
    rxns = {
        "R1": Reaction(id=1, substrates=frozenset([N1]), products=frozenset([M1])),
        "R2": Reaction(id=2, substrates=frozenset([M1]), products=frozenset([M2])),
        "R2b": Reaction(id=3, substrates=frozenset([M2]), products=frozenset([M_deep])),
        "R3": Reaction(id=4, substrates=frozenset([N1, N2]), products=frozenset([M4])),
        "R4": Reaction(id=5, substrates=frozenset([M4, U1]), products=frozenset([M5])),
        "R5": Reaction(id=6, substrates=frozenset([N2]), products=frozenset([M3a])),
        "R6": Reaction(id=7, substrates=frozenset([M5, M3a]), products=frozenset([T]),
                       ecnum="1.14.13.1"),
        "R7": Reaction(id=8, substrates=frozenset([N2]), products=frozenset([T]),
                       ecnum="2.7.1.1"),
        "R8": Reaction(id=9, substrates=frozenset([N1]), products=frozenset([M3b])),
        "R9": Reaction(id=10, substrates=frozenset([M1]), products=frozenset([M3b]),
                       ecnum="1.1.1.1"),
        "R10": Reaction(id=11, substrates=frozenset([F1]), products=frozenset([M6])),
        "R11": Reaction(id=12, substrates=frozenset([X_missing]), products=frozenset([M7])),
    }

    return SmallCell(
        chems=chems,
        rxns=rxns,
        natives={N1, N2},
        universals={U1},
        fed={F1},
    )


def write_small_cell_tsvs(sc: SmallCell, out_dir: Path) -> None:
    """Dump the fixture as TSV files matching parser.py formats.

    Used by the parser round-trip test and the end-to-end integration test.
    Keeping this in the fixture module ensures the TSV is always derived from
    the canonical Python fixture — no hand-maintained sync burden.
    """
    chem_lines = ["id\tname\tinchi\tsmiles"]
    for c in sc.chems.values():
        chem_lines.append(f"{c.id}\t{c.name}\t{c.inchi}\t{c.smiles or 'null'}")
    (out_dir / "enzymemap_chems.tsv").write_text("\n".join(chem_lines) + "\n")

    rxn_lines = ["rxnid\tecnum\tsubstrates\tproducts"]
    for r in sc.rxns.values():
        subs = " ".join(str(s.id) for s in r.substrates)
        prods = " ".join(str(p.id) for p in r.products)
        rxn_lines.append(f"{r.id}\t{r.ecnum or ''}\t{subs}\t{prods}")
    (out_dir / "enzymemap_reactions.tsv").write_text("\n".join(rxn_lines) + "\n")

    # Organism sidecar: leave every reaction untagged so the two-pass BFS
    # in MCP state.py reduces to single-pass over the full reaction set —
    # pass 1 enables nothing (no E. coli reactions), so its reachables are
    # just the natives ∪ universals, and pass 2 reproduces the original
    # MetaCyc-era shell numbering the small-cell tests assert against.
    (out_dir / "enzymemap_reaction_organisms.tsv").write_text(
        "rxnid\torganisms\n"
    )

    native_lines = ["name\tinchi\tdescriptor"]
    for c in sc.natives:
        native_lines.append(f"{c.name}\t{c.inchi}\tfeedstock")
    (out_dir / "minimal_metabolites.txt").write_text("\n".join(native_lines) + "\n")

    univ_lines = ["name\tinchi\tdescriptor"]
    for c in sc.universals:
        univ_lines.append(f"{c.name}\t{c.inchi}\tcofactor")
    (out_dir / "ubiquitous_metabolites.txt").write_text("\n".join(univ_lines) + "\n")

    # ec_names.tsv — map only a subset of the fixture's ECs so tests can
    # cover both named-enzyme and orphan-EC branches. R6's 1.14.13.1 is
    # intentionally OMITTED → becomes an orphan P450 step.
    ec_lines = [
        "ecnum\tname",
        "2.7.1.1\thexokinase",
        "1.1.1.1\talcohol dehydrogenase",
    ]
    (out_dir / "ec_names.tsv").write_text("\n".join(ec_lines) + "\n")
