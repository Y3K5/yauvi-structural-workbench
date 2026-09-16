"""The structure reader."""
from __future__ import annotations

import pytest

from actstate.structure import StructureError, read_structure


def test_reads_residues_and_prefers_cb(write_pdb, clustered_triad):
    structure = read_structure(write_pdb("x.pdb", clustered_triad))
    assert len(structure.residues) == 3
    assert {r.atom_used for r in structure.residues} == {"CB"}
    assert structure.by_seq_id()[5].name == "HIS"


def test_glycine_falls_back_to_ca(write_pdb):
    structure = read_structure(write_pdb("g.pdb", [(1, "GLY", (0.0, 0.0, 0.0))]))
    assert structure.residues[0].atom_used == "CA"


def test_distance_between_residues(write_pdb, clustered_triad):
    residues = read_structure(write_pdb("x.pdb", clustered_triad)).by_seq_id()
    assert residues[5].distance_to(residues[9]) == pytest.approx(4.0, abs=0.01)


def test_one_letter_translation(write_pdb):
    structure = read_structure(write_pdb("x.pdb", [(1, "HIS", (0.0, 0.0, 0.0))]))
    assert structure.residues[0].one_letter == "H"


def test_water_and_buffer_are_not_candidate_cofactors(write_pdb, clustered_triad):
    structure = read_structure(
        write_pdb(
            "x.pdb",
            clustered_triad,
            heteroatoms=[("HOH", (9.0, 9.0, 9.0), "O"), ("GOL", (8.0, 8.0, 8.0), "C")],
        )
    )
    assert len(structure.heteroatoms) == 2
    # Counting a cryoprotectant as occupancy is how an apo structure gets
    # miscalled holo.
    assert structure.candidate_cofactors() == []


def test_a_real_ligand_is_a_candidate_cofactor(write_pdb, clustered_triad):
    structure = read_structure(
        write_pdb("x.pdb", clustered_triad, heteroatoms=[("ZN", (2.0, 1.5, 1.0), "ZN")])
    )
    assert [h.name for h in structure.candidate_cofactors()] == ["ZN"]


def test_alphafold_header_marks_a_prediction(write_pdb, clustered_triad):
    structure = read_structure(
        write_pdb(
            "m.pdb",
            clustered_triad,
            header="HEADER    PREDICTED MODEL\nREMARK   1 ALPHAFOLD DB PREDICTION\n",
        )
    )
    assert structure.is_predicted
    assert "AlphaFold" in structure.source_note


def test_alphafold_filename_marks_a_prediction(write_pdb, clustered_triad):
    structure = read_structure(write_pdb("AF-P12345-F1.pdb", clustered_triad))
    assert structure.is_predicted


def test_an_experimental_header_is_not_a_prediction(write_pdb, clustered_triad):
    structure = read_structure(
        write_pdb("e.pdb", clustered_triad, header="HEADER    HYDROLASE    01-JAN-26   0AAA\n")
    )
    assert not structure.is_predicted


def test_only_the_first_model_is_read(tmp_path):
    path = tmp_path / "nmr.pdb"
    path.write_text(
        "MODEL        1\n"
        "ATOM      1  CA  HIS A   5       0.000   0.000   0.000  1.00 50.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "ATOM      2  CA  HIS A   5      99.000  99.000  99.000  1.00 50.00           C\n"
        "ENDMDL\n",
        encoding="utf-8",
    )
    structure = read_structure(path)
    assert structure.by_seq_id()[5].x == pytest.approx(0.0)


def test_missing_file_is_an_error(tmp_path):
    with pytest.raises(StructureError, match="not found"):
        read_structure(tmp_path / "absent.pdb")


def test_a_file_with_no_residues_is_an_error(tmp_path):
    path = tmp_path / "empty.pdb"
    path.write_text("HEADER    NOTHING\nEND\n", encoding="utf-8")
    with pytest.raises(StructureError, match="no amino-acid residues"):
        read_structure(path)


def test_malformed_coordinate_lines_are_skipped(tmp_path):
    path = tmp_path / "mixed.pdb"
    path.write_text(
        "ATOM      1  CA  HIS A   5     BAD_X   0.000   0.000  1.00 50.00           C\n"
        "ATOM      2  CA  ASP A   9       4.000   0.000   0.000  1.00 50.00           C\n",
        encoding="utf-8",
    )
    structure = read_structure(path)
    assert list(structure.by_seq_id()) == [9]


# --- mmCIF ---------------------------------------------------------------

MMCIF = """data_TEST
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
ATOM 1 CB HIS A 5 0.000 0.000 0.000 50.00
ATOM 2 CB ASP A 9 4.000 0.000 0.000 50.00
ATOM 3 CB SER A 14 2.000 3.500 0.000 50.00
HETATM 4 ZN ZN A 901 2.000 1.500 1.000 30.00
#
"""


def test_mmcif_is_read_by_column_name(tmp_path):
    path = tmp_path / "x.cif"
    path.write_text(MMCIF, encoding="utf-8")
    structure = read_structure(path)
    assert sorted(structure.by_seq_id()) == [5, 9, 14]
    assert [h.name for h in structure.candidate_cofactors()] == ["ZN"]


def test_mmcif_without_atom_site_is_an_error(tmp_path):
    path = tmp_path / "bad.cif"
    path.write_text("data_TEST\n#\n_entry.id TEST\n", encoding="utf-8")
    with pytest.raises(StructureError, match="no atom_site"):
        read_structure(path)
