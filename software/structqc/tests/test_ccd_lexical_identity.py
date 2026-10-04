"""Synthetic CIF identity and reference-relative geometry regressions."""
import copy
import hashlib
import json

import numpy as np
import pytest

from structqc import chirality, coordinate_scope


def fixture(folder, cid="ALA", style="bare", reflect=False, primed=False, ambiguous=False):
    quote = lambda s: str(s) if style == "bare" else "'" + str(s) + "'" if style == "single" else '"' + str(s) + '"'
    names = ["C1'", "N", "C2'", "CB"] if primed else ["CA", "N", "C", "CB"]
    xyz = [[0, 0, 0], [-1, 0, 0], [0, 1, 0], [0, 0, 1]]
    elements = ["C", "N", "C", "C"]
    pdb = folder / "synthetic.pdb"
    pdb.write_text(''.join(
        f'HETATM{i+1:5d} {name:^4} {cid:>3} A{1:4d}    {x:8.3f}{y:8.3f}{(-z if reflect else z):8.3f}{1:6.2f}{20:6.2f}          {element:>2}\n'
        for i, (name, (x, y, z), element) in enumerate(zip(names, xyz, elements))) + "END\n")
    reference = folder / "component.cif"
    lines = ["data_fixture", "_chem_comp.id " + quote(cid), "_chem_comp.type 'NON-POLYMER'",
             "_chem_comp.pdbx_ambiguous_flag " + quote("Y" if ambiguous else "N"), "loop_",
             "_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol", "_chem_comp_atom.pdbx_stereo_config",
             "_chem_comp_atom.pdbx_model_Cartn_x_ideal", "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
             "_chem_comp_atom.pdbx_model_Cartn_z_ideal"]
    lines += [' '.join([quote(n), quote(e), quote("S" if i == 0 else "N"), *map(str, pos)])
              for i, (n, e, pos) in enumerate(zip(names, elements, xyz))]
    lines += [quote("H") + " H N 1 -1 -1", "loop_", "_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2"]
    lines += [quote(names[0]) + " " + quote(n) for n in [*names[1:], "H"]]
    reference.write_text('\n'.join(lines) + '\n')
    manifest = folder / "references.json"
    manifest.write_text(json.dumps({"schema_version": "1.0", "components": [{
        "component_id": cid, "path": reference.name, "sha256": hashlib.sha256(reference.read_bytes()).hexdigest()}]}))
    return coordinate_scope.load_scope(pdb), manifest


def relock(manifest, reference):
    record = json.loads(manifest.read_text())
    record["components"][0]["sha256"] = hashlib.sha256(reference.read_bytes()).hexdigest()
    manifest.write_text(json.dumps(record))


@pytest.mark.parametrize("cid", ["ALA", "DAL", "LIG", "SUG", "LIP"])
@pytest.mark.parametrize("style", ["bare", "single", "double"])
@pytest.mark.parametrize("reflect", [False, True])
def test_lexical_quoting_does_not_change_geometry(tmp_path, cid, style, reflect):
    scope, refs = fixture(tmp_path, cid, style, reflect)
    finding = chirality.check_scope(scope, refs)["findings"][0]
    assert finding["state"] == ("reference_mismatch" if reflect else "matches_reference")
    assert finding["observed_configuration"] is None
    assert finding["atom_id"] == "CA"
    assert finding["neighbor_atom_ids"] == ["C", "CB", "N"]


@pytest.mark.parametrize("style", ["bare", "single", "double"])
def test_apostrophes_remain_part_of_atom_identity(tmp_path, style):
    scope, refs = fixture(tmp_path, "LIG", style, primed=True)
    references, _ = chirality.read_references(refs)
    assert {"C1'", "C2'"} <= set(references["LIG"]["atoms"])
    assert chirality.check_scope(scope, refs)["findings"][0]["state"] == "matches_reference"


@pytest.mark.parametrize("style", ["single", "double"])
def test_quoted_ambiguity_is_not_a_clean_result(tmp_path, style):
    scope, refs = fixture(tmp_path, "LIG", style, ambiguous=True)
    finding = chirality.check_scope(scope, refs)["findings"][0]
    assert (finding["state"], finding["reason"]) == ("unevaluated", "ambiguous_reference_chemistry")


def test_missing_atoms_alternates_and_references_stay_unevaluated(tmp_path):
    scope, refs = fixture(tmp_path, "SUG", "double", primed=True)
    missing = copy.deepcopy(scope)
    component = missing["components"][0]
    component["atom_indices"] = [i for i in component["atom_indices"] if missing["atoms"][i]["atom"] != "C2'"]
    finding = chirality.check_scope(missing, refs)["findings"][0]
    assert (finding["state"], finding["reason"]) == ("unevaluated", "missing_required_atoms")
    alternate = copy.deepcopy(scope)
    alternate["components"][0]["alternate_locations"] = ["A", "B"]
    assert chirality.check_scope(alternate, refs)["findings"][0]["reason"] == "alternative_conformers"
    assert chirality.check_scope(scope)["state"] == "partial"


def test_duplicate_decoded_atom_id_rejected(tmp_path):
    _, refs = fixture(tmp_path, "LIG", "double")
    reference = tmp_path / "component.cif"
    reference.write_text(reference.read_text().replace('"H" H N 1 -1 -1', 'CA C N 1 -1 -1'))
    relock(refs, reference)
    with pytest.raises(ValueError, match="Duplicate CCD atom identity"):
        chirality.read_references(refs)


def test_checksum_change_rejected_before_decoding(tmp_path):
    scope, refs = fixture(tmp_path)
    record = json.loads(refs.read_text())
    record["components"][0]["sha256"] = "0" * 64
    refs.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="checksum mismatch"):
        chirality.check_scope(scope, refs)


def test_proper_transform_preserves_handedness_and_reflection_path_rejected(tmp_path):
    scope, refs = fixture(tmp_path, "LIG", "single", primed=True)
    rotation = coordinate_scope.proper_rotation([[0, -1, 0], [1, 0, 0], [0, 0, 1]], [20, -4, 8])
    moved = copy.deepcopy(scope)
    for atom in moved["atoms"]:
        xyz = rotation @ np.array([atom[x] for x in "xyz"]) + [20, -4, 8]
        atom.update(zip("xyz", xyz))
    before = chirality.check_scope(scope, refs)["findings"][0]
    after = chirality.check_scope(moved, refs)["findings"][0]
    assert before["state"] == after["state"]
    assert after["observed_volume_A3"] == pytest.approx(before["observed_volume_A3"], abs=1e-9)
    with pytest.raises(ValueError):
        coordinate_scope.proper_rotation(np.diag([-1, 1, 1]))


def test_null_tokens_and_unknown_stereochemistry_remain_unevaluated(tmp_path):
    scope, refs = fixture(tmp_path)
    reference = tmp_path / "component.cif"
    reference.write_text(reference.read_text().replace("CA C S 0 0 0", "CA C ? ? . 0"))
    relock(refs, reference)
    references, _ = chirality.read_references(refs)
    atom = references["ALA"]["atoms"]["CA"]
    assert atom["pdbx_model_Cartn_x_ideal"] == "?"
    assert atom["pdbx_model_Cartn_y_ideal"] == "."
    assert chirality.check_scope(scope, refs)["state"] == "partial"
