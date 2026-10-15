from __future__ import annotations

import json

from pose_evidence.cli import main
from pose_evidence.core import analyze


def atom(serial, name, resname, chain, resseq, x, y, z, element="C"):
    return f"ATOM  {serial:5d} {name:^4s} {resname:>3s} {chain:1s}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}\n"


def test_atom_order_does_not_change_pose_evidence(tmp_path):
    receptor = tmp_path / "receptor.pdb"; receptor.write_text(atom(1, "CA", "ALA", "A", 1, 0, 0, 0) + atom(2, "N", "ALA", "A", 1, 0, 1, 0, "N"))
    lines = [atom(1, "C1", "LIG", "L", 1, 2.5, 0, 0), atom(2, "O1", "LIG", "L", 1, 2.5, 1, 0, "O")]
    first = tmp_path / "first.pdbqt"; first.write_text("".join(lines))
    second = tmp_path / "second.pdbqt"; second.write_text("".join(reversed(lines)))
    base = {"subject": {"id": "p"}, "receptor": {"path": receptor.name}, "poses": [
        {"pose_id": "p1", "path": first.name, "engine": "autodock_vina", "job_type": "small_molecule", "ligand_chains": ["L"], "engine_score": -7, "engine_score_unit": "kcal/mol"},
        {"pose_id": "p2", "path": second.name, "engine": "autodock_vina", "job_type": "small_molecule", "ligand_chains": ["L"], "engine_score": -6.9, "engine_score_unit": "kcal/mol"}]}
    document, _ = analyze(base, base_dir=tmp_path)
    assert len({row["cluster_id"] for row in document["clusters"]}) == 1
    assert document["poses"][0]["contact_count"] == document["poses"][1]["contact_count"]


def test_cli_exit_one_preserves_missing_sasa(tmp_path):
    receptor = tmp_path / "r.pdb"; receptor.write_text(atom(1, "CA", "ALA", "A", 1, 0, 0, 0))
    pose = tmp_path / "p.pdb"; pose.write_text(atom(1, "C1", "LIG", "L", 1, 3, 0, 0))
    config = tmp_path / "config.json"; config.write_text(json.dumps({"subject": {"id": "p"}, "receptor": {"path": "r.pdb"},
        "poses": [{"pose_id": "p1", "path": "p.pdb", "engine": "cluspro", "job_type": "protein_peptide", "ligand_chains": ["L"], "engine_score_unit": "ClusPro cluster rank"}]}))
    out = tmp_path / "out"
    assert main(["run", "--config", str(config), "--out", str(out)]) == 1
    assert (out / "POSE_EVIDENCE.json").is_file()
    assert "affinity" not in json.loads((out / "POSE_EVIDENCE.json").read_text())["limitations"][0].lower()
