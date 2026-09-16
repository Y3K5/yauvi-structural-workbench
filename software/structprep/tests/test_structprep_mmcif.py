"""mmCIF support tests.

These are deliberately self-contained. The PDB tests in `test_structprep.py` skip
unless STRUCTPREP_TEST_STRUCTURES points at unshipped RCSB coordinates, so they do
not run in a clean checkout; the fixtures here are synthetic and always run.

One test encodes a mistake made while adding this support: `label_asym_id` was
rewritten along with `auth_asym_id`, which silently reassigned an entity label on
every ligand whose two identifiers differ. It changed 52 atom records in a real
cryo-EM entry while the run record claimed only a chain relabel.
"""
import json
import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from structprep.core import InputError, prepare, write_outputs

TAGS = """loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""

# The ZN sits 0.5 A from CA, so it is "contacting" and must be refused, not removed.
# Its label_asym_id is G while its auth_asym_id is A — the case that caused the bug.
ROWS_ONE_MODEL = """ATOM   1 N  N   . ALA A 1 1 ? 0.000 0.000 0.000 1.00 10.00 1   ALA A N  1
ATOM   2 C  CA  . ALA A 1 1 ? 1.500 0.000 0.000 1.00 10.00 1   ALA A CA 1
HETATM 3 O  O   . HOH B 2 . ? 50.000 50.000 50.000 1.00 20.00 101 HOH A O  1
HETATM 4 ZN ZN  . ZN  G 3 . ? 2.000 0.000 0.000 1.00 15.00 201 ZN  A ZN 1
HETATM 5 C  C1  . GOL C 4 . ? 1.600 1.000 0.000 1.00 25.00 301 GOL A C1 1
"""

ROWS_TWO_MODELS = ROWS_ONE_MODEL + """ATOM   6 N  N   . ALA A 1 1 ? 9.000 0.000 0.000 1.00 10.00 1   ALA A N  2
ATOM   7 C  CA  . ALA A 1 1 ? 10.500 0.000 0.000 1.00 10.00 1  ALA A CA 2
HETATM 8 ZN ZN  . ZN  G 3 . ? 11.000 0.000 0.000 1.00 15.00 201 ZN A ZN 2
"""


def _cif(rows: str) -> str:
    return ("data_TEST\n#\n_entry.id TEST\n#\n" + TAGS + rows
            + "#\n_struct.title 'content after the loop'\n#\n")


def write(tmp_path, text, name="case.cif"):
    path = tmp_path / name
    path.write_text(text)
    return path


# --- the contract the PDB path gives, on mmCIF ---------------------------------

def test_kept_atom_rows_are_reemitted_byte_for_byte(tmp_path):
    src = _cif(ROWS_ONE_MODEL)
    result = prepare(write(tmp_path, src))
    write_outputs(result, tmp_path / "out")
    prepared = (tmp_path / "out" / "PREPARED.cif").read_text()

    source_rows = {l for l in src.splitlines() if l.startswith(("ATOM", "HETATM"))}
    kept_rows = [l for l in prepared.splitlines() if l.startswith(("ATOM", "HETATM"))]
    assert kept_rows, "no atom rows survived"
    assert all(row in source_rows for row in kept_rows)


def test_everything_outside_the_loop_survives(tmp_path):
    src = _cif(ROWS_ONE_MODEL)
    result = prepare(write(tmp_path, src))
    write_outputs(result, tmp_path / "out")
    prepared = (tmp_path / "out" / "PREPARED.cif").read_text()

    assert "_entry.id TEST" in prepared                      # before the loop
    assert "_struct.title 'content after the loop'" in prepared   # after it
    before = [l for l in src.splitlines() if not l.startswith(("ATOM", "HETATM"))]
    after = [l for l in prepared.splitlines() if not l.startswith(("ATOM", "HETATM"))]
    assert before == after


def test_output_is_named_and_recorded_as_mmcif(tmp_path):
    result = prepare(write(tmp_path, _cif(ROWS_ONE_MODEL)))
    write_outputs(result, tmp_path / "out")
    assert (tmp_path / "out" / "PREPARED.cif").is_file()
    assert not (tmp_path / "out" / "PREPARED.pdb").exists()
    record = json.loads((tmp_path / "out" / "PREPARATION.json").read_text())
    assert record["coordinate_format"] == "mmcif"
    manifest = json.loads((tmp_path / "out" / "RUN_MANIFEST.json").read_text())
    assert "PREPARED.cif" in manifest["outputs"]


def test_solvent_goes_contacting_additive_is_refused_ion_is_kept(tmp_path):
    """The same decisions the PDB path makes, reached from mmCIF columns.

    HOH is solvent and far from the polymer, so it goes. GOL is a cryoprotectant —
    in DEFAULT_DROP — but sits 1.2 A from CA, so it is refused rather than removed.
    ZN is an ion, which is not dropped by default, so it is simply kept.
    """
    result = prepare(write(tmp_path, _cif(ROWS_ONE_MODEL)))
    removed = {a["resname"] for a in result["removed"]}
    kept = {a["resname"] for a in result["kept"]}
    assert "HOH" in removed
    assert {"ZN", "GOL"} <= kept
    assert [r["component"] for r in result["refusals"]] == ["GOL"]


def test_chain_selection_reads_auth_asym_id(tmp_path):
    result = prepare(write(tmp_path, _cif(ROWS_ONE_MODEL)), {"chains": ["A"]})
    assert result["counts"]["atoms_kept"] == 4      # ZN and GOL are auth chain A
    result_b = prepare(write(tmp_path, _cif(ROWS_ONE_MODEL)), {"chains": ["G"]})
    assert result_b["counts"]["atoms_kept"] == 0


# --- the regression -------------------------------------------------------------

def test_flatten_relabels_auth_asym_id_and_leaves_label_asym_id_alone(tmp_path):
    """label_asym_id is an entity label, not a chain name, and the two differ.

    Rewriting both made the ZN of auth chain A claim entity A as well, changing the
    deposited entity assignment under a record that said only chains were relabelled.
    """
    result = prepare(write(tmp_path, _cif(ROWS_TWO_MODELS)))
    write_outputs(result, tmp_path / "out")
    prepared = (tmp_path / "out" / "PREPARED.cif").read_text()
    zinc = [l.split() for l in prepared.splitlines() if " ZN " in f" {l} "]
    assert len(zinc) == 2, "both models' zinc atoms should survive"

    label_ids = {row[6] for row in zinc}
    auth_ids = {row[17] for row in zinc}
    assert label_ids == {"G"}, f"label_asym_id was rewritten: {label_ids}"
    assert len(auth_ids) == 2, f"auth_asym_id was not flattened to unique ids: {auth_ids}"
    assert result["counts"]["models_in"] == 2


def test_untouched_rows_are_not_rewritten_when_nothing_is_flattened(tmp_path):
    """A single-model file needs no relabelling, so every row is the source row."""
    src = _cif(ROWS_ONE_MODEL)
    result = prepare(write(tmp_path, src))
    write_outputs(result, tmp_path / "out")
    prepared = (tmp_path / "out" / "PREPARED.cif").read_text()
    for row in src.splitlines():
        if row.startswith("ATOM") or (row.startswith("HETATM") and " ZN " in f" {row} "):
            assert row in prepared


# --- refusals rather than guesses -----------------------------------------------

def test_multiline_value_in_the_loop_is_refused(tmp_path):
    rows = ROWS_ONE_MODEL + ";a multi-line value\n;\n"
    with pytest.raises(InputError, match="multi-line"):
        prepare(write(tmp_path, _cif(rows)))


def test_two_atom_site_loops_are_refused(tmp_path):
    text = _cif(ROWS_ONE_MODEL) + TAGS + ROWS_ONE_MODEL + "#\n"
    with pytest.raises(InputError, match="more than one"):
        prepare(write(tmp_path, text))


def test_missing_atom_site_loop_is_refused(tmp_path):
    with pytest.raises(InputError, match="no _atom_site loop"):
        prepare(write(tmp_path, "data_TEST\n#\n_entry.id TEST\n#\n"))


def test_ragged_row_count_is_refused(tmp_path):
    rows = ROWS_ONE_MODEL + "ATOM 5 N N . ALA A 1 1 ? 0.0 0.0\n"
    with pytest.raises(InputError, match="not whole"):
        prepare(write(tmp_path, _cif(rows)))


def test_missing_coordinate_column_is_refused(tmp_path):
    tags = TAGS.replace("_atom_site.Cartn_z\n", "")
    rows = "\n".join(" ".join(l.split()[:12] + l.split()[13:])
                     for l in ROWS_ONE_MODEL.strip().splitlines()) + "\n"
    with pytest.raises(InputError, match="Cartn_z"):
        prepare(write(tmp_path, "data_T\n#\n" + tags + rows + "#\n"))


def test_quoted_value_is_parsed_and_reemitted(tmp_path):
    """A quoted token holds a space, so the row's column count only survives if the
    tokeniser respects the quotes. auth_atom_id is read before label_atom_id."""
    rows = ROWS_ONE_MODEL.replace("201 ZN  A ZN 1", "201 ZN  A 'ZN A' 1")
    path = write(tmp_path, _cif(rows))
    result = prepare(path)
    names = {a["name"] for a in result["kept"] + result["removed"]}
    assert "ZN A" in names
    write_outputs(result, tmp_path / "out")
    assert "'ZN A'" in (tmp_path / "out" / "PREPARED.cif").read_text()
