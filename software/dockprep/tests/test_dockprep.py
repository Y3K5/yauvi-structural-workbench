from __future__ import annotations

import json
from pathlib import Path

from dockprep.cli import main
from dockprep.core import ScientificBlock, prepare_job, sha256


def _inputs(tmp_path: Path, ligand_suffix: str = ".pdbqt") -> dict:
    receptor = tmp_path / "receptor.pdbqt"; receptor.write_text("ATOM\n")
    ligand = tmp_path / f"ligand{ligand_suffix}"; ligand.write_text("MOLECULE\n")
    approval = tmp_path / "approval.json"
    approval.write_text(json.dumps({"approval_id": "approval-1", "status": "approved_for_docking", "job_type": "small_molecule",
                                    "receptor_sha256": sha256(receptor), "ligand_sha256": sha256(ligand)}))
    return {"job_type": "small_molecule",
            "receptor": {"path": receptor.name, "assembly_id": "1", "conformational_state": "observed",
                         "selected_site_residues": [{"chain_id": "A", "auth_seq_id": 10, "mapping_state": "exact"}]},
            "ligand": {"path": ligand.name, "identifier": "CHEBI:1", "chemistry": {"identifier": "CHEBI:1", "connectivity": "declared",
                       "bond_orders": "declared", "charge": 0, "protonation": "pH_user_declared"}},
            "source_terms": {"derivatives_allowed": True, "docking_allowed": True}, "approval": {"path": approval.name},
            "search": {"center": [0, 0, 0], "size": [20, 20, 20], "seed": 7}}


def test_approved_job_is_checksum_bound_and_path_clean(tmp_path):
    spec = _inputs(tmp_path)
    out = tmp_path / "out"
    document = prepare_job(spec, base_dir=tmp_path, out_dir=out)
    assert document["approval"]["status"] == "approved_for_docking"
    assert document["receptor"]["path"].startswith("inputs/")
    assert str(tmp_path) not in json.dumps(document)


def test_ambiguous_site_blocks(tmp_path):
    spec = _inputs(tmp_path)
    spec["receptor"]["selected_site_residues"][0]["mapping_state"] = "ambiguous"
    try:
        prepare_job(spec, base_dir=tmp_path, out_dir=tmp_path / "out")
    except ScientificBlock as exc:
        assert "exact residue mapping" in str(exc)
    else:
        raise AssertionError("ambiguous mapping should block docking")


def test_coordinate_only_ligand_blocks(tmp_path):
    spec = _inputs(tmp_path, ".pdb")
    try:
        prepare_job(spec, base_dir=tmp_path, out_dir=tmp_path / "out")
    except ScientificBlock as exc:
        assert "coordinate-only" in str(exc)
    else:
        raise AssertionError("coordinate-only ligand should block docking")
