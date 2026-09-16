#!/usr/bin/env python3
"""Regenerate the bundled example fixtures.

The fixtures are synthetic on purpose. They are built to exercise every label in
the closed vocabulary exactly once, they are small enough to read, and they carry
no licence encumbrance — which is what lets `actstate run` and the golden tests
work with no network and no downloaded database.

Run from the `activity-state` directory:

    python3 examples/build_fixtures.py
"""
from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
STRUCTURES = HERE / "structures"


def pdb_atom(serial, name, res_name, chain, seq_id, x, y, z, b=50.0, element="C"):
    return (
        f"ATOM  {serial:>5} {name:<4}{res_name:>4} {chain}{seq_id:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00{b:>6.2f}          {element:>2}\n"
    )


def hetatm(serial, name, res_name, chain, seq_id, x, y, z, element="ZN"):
    return (
        f"HETATM{serial:>5} {name:<4}{res_name:>4} {chain}{seq_id:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 30.00          {element:>2}\n"
    )


def build_structure(path, header, residues, heteroatoms=(), b_factor=50.0):
    """residues: list of (seq_id, three_letter, (x, y, z))."""
    lines = [header]
    serial = 1
    for seq_id, res_name, (x, y, z) in residues:
        # A CA and a CB, so the reader's CB-preferred path is exercised.
        lines.append(pdb_atom(serial, "CA", res_name, "A", seq_id, x, y, z, b_factor, "C"))
        serial += 1
        if res_name != "GLY":
            lines.append(
                pdb_atom(serial, "CB", res_name, "A", seq_id, x + 0.5, y, z, b_factor, "C")
            )
            serial += 1
    for index, (res_name, (x, y, z), element) in enumerate(heteroatoms, start=1):
        lines.append(hetatm(serial, element, res_name, "A", 900 + index, x, y, z, element))
        serial += 1
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def spread(seq_ids, spacing):
    """Place residues along a line at a fixed spacing."""
    return [(seq_id, i * spacing) for i, seq_id in enumerate(seq_ids)]


def main() -> None:
    STRUCTURES.mkdir(parents=True, exist_ok=True)

    # Backbone filler so each structure has more than just the catalytic residues.
    def scaffold(exclude):
        return [
            (i, "ALA", (float(i) * 3.0, 25.0, 0.0))
            for i in range(1, 21)
            if i not in exclude
        ]

    # 1. P_ACTIVE — catalytic triad clustered, zinc present, experimental header.
    triad = {5: "HIS", 9: "ASP", 14: "SER"}
    build_structure(
        STRUCTURES / "P_ACTIVE.pdb",
        "HEADER    HYDROLASE                               01-JAN-26   0AAA\n",
        scaffold(triad)
        + [
            (5, "HIS", (0.0, 0.0, 0.0)),
            (9, "ASP", (4.0, 0.0, 0.0)),
            (14, "SER", (2.0, 3.5, 0.0)),
        ],
        heteroatoms=[("ZN", (2.0, 1.5, 1.0), "ZN")],
    )

    # 2. P_PREDICTED — same intact geometry, but an AlphaFold model.
    build_structure(
        STRUCTURES / "AF-P_PREDICTED-F1.pdb",
        "HEADER    PREDICTED MODEL                         01-JAN-26   AFDB\n"
        "REMARK   1 ALPHAFOLD DB PREDICTION\n",
        scaffold(triad)
        + [
            (5, "HIS", (0.0, 0.0, 0.0)),
            (9, "ASP", (4.0, 0.0, 0.0)),
            (14, "SER", (2.0, 3.5, 0.0)),
        ],
        b_factor=92.4,  # pLDDT
    )

    # 3. P_APO — intact site, declares zinc, no heteroatom but ordered water.
    build_structure(
        STRUCTURES / "P_APO.pdb",
        "HEADER    HYDROLASE                               01-JAN-26   0CCC\n",
        scaffold(triad)
        + [
            (5, "HIS", (0.0, 0.0, 0.0)),
            (9, "ASP", (4.0, 0.0, 0.0)),
            (14, "SER", (2.0, 3.5, 0.0)),
        ],
        # Water and glycerol only: neither is cofactor occupancy, and treating
        # them as such is the standard way to call an apo structure holo.
        heteroatoms=[("HOH", (8.0, 8.0, 8.0), "O"), ("GOL", (9.0, 9.0, 9.0), "C")],
    )

    # 4. P_DISPERSED — the same three residues, far apart.
    build_structure(
        STRUCTURES / "P_DISPERSED.pdb",
        "HEADER    HYDROLASE                               01-JAN-26   0DDD\n",
        scaffold(triad)
        + [
            (5, "HIS", (0.0, 0.0, 0.0)),
            (9, "ASP", (30.0, 0.0, 0.0)),
            (14, "SER", (0.0, 30.0, 0.0)),
        ],
    )

    # 5. P_DISRUPTED — annotated catalytic positions, but the sequence has Ala there.
    build_structure(
        STRUCTURES / "P_DISRUPTED.pdb",
        "HEADER    PSEUDOENZYME                            01-JAN-26   0EEE\n",
        scaffold(triad)
        + [
            (5, "ALA", (0.0, 0.0, 0.0)),
            (9, "ALA", (4.0, 0.0, 0.0)),
            (14, "SER", (2.0, 3.5, 0.0)),
        ],
    )

    # 6. P_NOSITE has no structure at all, and no annotated site.

    # --- sequences -------------------------------------------------------
    # 20 residues each. Positions 5, 9 and 14 carry the catalytic residues.
    def sequence(p5, p9, p14):
        seq = list("A" * 20)
        seq[4], seq[8], seq[13] = p5, p9, p14
        return "".join(seq)

    sequences = {
        "P_ACTIVE": sequence("H", "D", "S"),
        "P_PREDICTED": sequence("H", "D", "S"),
        "P_APO": sequence("H", "D", "S"),
        "P_DISPERSED": sequence("H", "D", "S"),
        "P_DISRUPTED": sequence("A", "A", "S"),   # the degraded site
        "P_NOSITE": "A" * 20,
    }
    (HERE / "sequences.fasta").write_text(
        "".join(f">sp|{acc}|EXAMPLE\n{seq}\n" for acc, seq in sequences.items()),
        encoding="utf-8",
    )

    # --- annotation table ------------------------------------------------
    act = (
        'ACT_SITE 5; /note="Charge relay system"; /evidence="ECO:0000269|PubMed:1"; '
        'ACT_SITE 9; /note="Charge relay system"; /evidence="ECO:0000269|PubMed:1"; '
        'ACT_SITE 14; /note="Nucleophile"; /evidence="ECO:0000269|PubMed:1"'
    )
    cofactor_zn = 'COFACTOR: Name=Zn(2+); Xref=ChEBI:CHEBI:29105; Evidence=ECO:0000250;'
    binding = 'BINDING 5..9; /ligand="Zn(2+)"; /ligand_id="ChEBI:CHEBI:29105"'

    rows = [
        ("P_ACTIVE", act, binding, cofactor_zn, "3.4.21.1", "Example active hydrolase"),
        ("P_PREDICTED", act, "", "", "3.4.21.1", "Example predicted hydrolase"),
        ("P_APO", act, binding, cofactor_zn, "3.4.21.1", "Example apo hydrolase"),
        ("P_DISPERSED", act, "", "", "3.4.21.1", "Example dispersed-site hydrolase"),
        ("P_DISRUPTED", act, "", "", "", "Example pseudoenzyme"),
        ("P_NOSITE", "", "", "", "", "Example protein with no annotated site"),
    ]
    header = "Entry\tft_act_site\tft_binding\tcc_cofactor\tec\tProtein names\n"
    (HERE / "annotations.tsv").write_text(
        header + "".join("\t".join(r) + "\n" for r in rows), encoding="utf-8"
    )

    # --- optional sidecars ----------------------------------------------
    # Supplied only for P_ACTIVE, so the fixture also demonstrates that the other
    # entries correctly report these signals as unavailable rather than neutral.
    (HERE / "fold_state.json").write_text(
        '{\n  "P_ACTIVE": {"state": "active_assembly", "source": "example fixture"}\n}\n',
        encoding="utf-8",
    )
    (HERE / "reference_comparison.json").write_text(
        '{\n  "P_ACTIVE": {"reference": "0REF", "state": "active", "score": 0.94}\n}\n',
        encoding="utf-8",
    )
    # Supplied only for P_DISRUPTED, because it is the only fixture making the
    # disruption claim. Without it that entry reports its Ala positions as
    # contradicting evidence and stops at `indeterminate` -- which is the point
    # of the fixture as much as the label is.
    (HERE / "expected_residues.json").write_text(
        '{\n  "P_DISRUPTED": {\n    "5": "H",\n    "9": "D",\n    "14": "S"\n  }\n}\n',
        encoding="utf-8",
    )

    print(f"wrote fixtures under {HERE}")


if __name__ == "__main__":
    main()
