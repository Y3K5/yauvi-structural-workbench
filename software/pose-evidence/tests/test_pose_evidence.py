from __future__ import annotations

import csv
from importlib.resources import files
import json
import os
from pathlib import Path

import pytest

from pose_evidence.cli import describe, main
from pose_evidence.core import analyze, sha256, write_outputs


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


LIMITATIONS = [
    "Docking scores remain engine-specific and are not affinities.",
    "Hydrogen bonds, salt bridges, clashes, and metal coordination are geometry proxies.",
    "Imported ClusPro results are read-only and are never uploaded by this module.",
]
SASA_REASON = "FreeSASA or an explicitly validated SASA adapter was not invoked"


def _pdb(path: Path, rows: list[tuple]) -> None:
    path.write_text("".join(atom(*row) for row in rows), encoding="utf-8")


def _pose(pose_id: str, path: Path, *, engine: str = "autodock_vina", job_type: str = "small_molecule",
          chains: list[str] | None = None, score: float | None = -1.0, unit: str = "kcal/mol") -> dict:
    return {"pose_id": pose_id, "path": path.name, "engine": engine, "job_type": job_type,
            "ligand_chains": ["L"] if chains is None else chains, "engine_score": score, "engine_score_unit": unit}


def _config(tmp_path: Path, document: dict, name: str = "config.json") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(document), encoding="utf-8")
    return path


def _cli(capsys, argv: list[str]) -> tuple[int, str, str]:
    code = main(argv)
    captured = capsys.readouterr()
    return code, captured.out, captured.err


def _tree(root: Path) -> dict[str, bytes]:
    return {path.relative_to(root).as_posix(): path.read_bytes() for path in sorted(root.rglob("*")) if path.is_file()}


def _simple(tmp_path: Path) -> Path:
    _pdb(tmp_path / "receptor.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    _pdb(tmp_path / "pose.pdb", [(1, "C1", "LIG", "L", 1, 3, 0, 0, "C")])
    return _config(tmp_path, {"subject": {"id": "synthetic"}, "receptor": {"path": "receptor.pdb", "sha256": sha256(tmp_path / "receptor.pdb")},
                              "poses": [_pose("p1", tmp_path / "pose.pdb")]})


def test_describe_lists_geometry_limits(capsys):
    code, out, err = _cli(capsys, ["describe"])
    assert code == 0 and err == ""
    payload = json.loads(out)
    assert payload == describe()
    assert payload["commands"] == ["describe", "validate", "fetch", "run"]
    assert payload["deterministic_output"] is True
    assert "geometry proxies" in payload["limitations"][2]


def test_fetch_prints_the_packaged_source_plan(capsys):
    code, out, err = _cli(capsys, ["fetch", "--plan"])
    assert code == 0 and err == ""
    assert out == (files("pose_evidence") / "sources.yaml").read_text()
    assert "cluspro" in out


def test_validate_counts_poses_and_does_not_write_outputs(tmp_path, capsys):
    config = _simple(tmp_path)
    code, out, err = _cli(capsys, ["validate", "--config", str(config)])
    payload = json.loads(out)
    assert code == 1 and err == ""
    assert payload == {"clusters": 1, "poses": 1, "valid": True}
    assert not (tmp_path / "POSE_EVIDENCE.json").exists()


def test_run_keeps_the_claim_ceiling_and_leaves_sasa_unevaluated(tmp_path, capsys):
    config = _simple(tmp_path)
    out = tmp_path / "out"
    code, _stdout, err = _cli(capsys, ["run", "--config", str(config), "--out", str(out)])
    assert code == 1 and err == ""
    document = json.loads((out / "POSE_EVIDENCE.json").read_text(encoding="utf-8"))
    manifest = json.loads((out / "RUN_MANIFEST.json").read_text(encoding="utf-8"))
    assert document["limitations"] == LIMITATIONS
    assert document["buried_sasa"] == {"state": "unevaluated", "reason": SASA_REASON}
    assert manifest["missing_evidence"] == ["buried_sasa"]
    assert manifest["optional_runtimes"] == {"freesasa": "not_invoked"}
    for name in ("POSE_CLUSTERS.tsv", "POSE_CONTACTS.tsv", "POSE_LAYER.json"):
        assert (out / name).is_file()


def test_hand_built_pose_records_contacts_clashes_and_proxies(tmp_path):
    _pdb(tmp_path / "receptor.pdb", [
        (1, "N", "ALA", "A", 1, 0, 0, 0, "N"),
        (2, "OD1", "ASP", "A", 2, 10, 0, 0, "O"),
        (3, "ZN", "ZN", "A", 3, 20, 0, 0, "ZN"),
    ])
    _pdb(tmp_path / "pose.pdb", [
        (1, "O1", "LIG", "L", 1, 3, 0, 0, "O"),
        (2, "NZ", "LYS", "L", 2, 14, 0, 0, "N"),
        (3, "C1", "LIG", "L", 3, 21.2, 0, 0, "C"),
    ])
    config = {"subject": {"id": "synthetic"}, "receptor": {"path": "receptor.pdb"}, "contact_cutoff_A": 5.0,
              "poses": [_pose("p1", tmp_path / "pose.pdb", score=-4.5)]}
    document, digests = analyze(config, base_dir=tmp_path)
    pose = document["poses"][0]
    by_ligand = {row["ligand_auth_seq_id"]: row for row in document["contacts"]}
    assert set(by_ligand) == {1, 2, 3}
    assert by_ligand[1]["contact_state"] == "contact"
    assert by_ligand[1]["hydrogen_bond_proxy"] is True
    assert by_ligand[1]["minimum_distance_A"] == 3.0
    assert by_ligand[2]["salt_bridge_proxy"] is True
    assert by_ligand[2]["hydrogen_bond_proxy"] is False
    assert by_ligand[2]["minimum_distance_A"] == 4.0
    assert by_ligand[3]["contact_state"] == "clash"
    assert by_ligand[3]["metal_coordination_proxy"] is True
    assert by_ligand[3]["minimum_distance_A"] == 1.2
    assert pose["contact_count"] == 3
    assert pose["clash_count"] == 1
    assert pose["hydrogen_bond_proxy_count"] == 1
    assert pose["salt_bridge_proxy_count"] == 1
    assert pose["metal_coordination_proxy_count"] == 1
    assert document["buried_sasa"]["state"] == "unevaluated"
    assert digests["receptor"] == sha256(tmp_path / "receptor.pdb")
    out = tmp_path / "out"
    write_outputs(out, document, input_sha256=digests, config_sha256="a" * 64)
    with (out / "POSE_CONTACTS.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["contact_state"] for row in rows] == ["contact", "contact", "clash"]
    layer = json.loads((out / "POSE_LAYER.json").read_text(encoding="utf-8"))
    assert {row["metric"] for row in layer["records"]} == {"docking_contact_count"}
    assert all(row["evidence_class"] == "docking_geometry" for row in layer["records"])


def test_small_molecule_poses_cluster_by_ligand_rmsd(tmp_path):
    _pdb(tmp_path / "receptor.pdb", [(1, "CA", "ALA", "A", 1, 0, 100, 0, "C")])
    for name, x in (("near.pdb", 0.0), ("boundary.pdb", 2.0), ("far.pdb", 5.0)):
        _pdb(tmp_path / name, [(1, "C1", "LIG", "L", 1, x, 0, 0, "C")])
    config = {"subject": {"id": "synthetic"}, "receptor": {"path": "receptor.pdb"}, "ligand_rmsd_cutoff_A": 2.0,
              "poses": [_pose("p1", tmp_path / "near.pdb"), _pose("p2", tmp_path / "boundary.pdb", score=-0.5),
                        _pose("p3", tmp_path / "far.pdb", score=-0.1)]}
    document, _digests = analyze(config, base_dir=tmp_path)
    clusters = {row["pose_id"]: row for row in document["clusters"]}
    assert clusters["p1"]["cluster_id"] == clusters["p2"]["cluster_id"] == "cluster_001"
    assert clusters["p3"]["cluster_id"] == "cluster_002"
    assert {row["method"] for row in document["clusters"]} == {"ligand_heavy_atom_rmsd_le_2.0_A"}
    assert document["repeated_seed_consistency"] == [
        {"engine": "autodock_vina", "receptor_state": "unresolved", "distinct_clusters": 2},
    ]


def test_macromolecular_poses_cluster_by_common_contacts(tmp_path):
    _pdb(tmp_path / "receptor.pdb", [
        (1, "CA", "ALA", "A", 1, 0, 0, 0, "C"),
        (2, "CA", "ALA", "A", 2, 10, 0, 0, "C"),
        (3, "CA", "ALA", "A", 3, 20, 0, 0, "C"),
    ])
    _pdb(tmp_path / "left.pdb", [
        (1, "CA", "ALA", "B", 1, 4, 0, 0, "C"),
        (2, "CA", "ALA", "B", 2, 14, 0, 0, "C"),
    ])
    _pdb(tmp_path / "overlap.pdb", [
        (1, "CA", "ALA", "B", 1, 4, 0, 0, "C"),
        (2, "CA", "ALA", "B", 2, 24, 0, 0, "C"),
    ])
    _pdb(tmp_path / "third.pdb", [(1, "CA", "ALA", "B", 2, 24, 0, 0, "C")])
    poses = [
        _pose("p1", tmp_path / "left.pdb", engine="haddock3", job_type="protein_protein", chains=["B"], score=None, unit="haddock"),
        _pose("p2", tmp_path / "overlap.pdb", engine="haddock3", job_type="protein_protein", chains=["B"], score=None, unit="haddock"),
        _pose("p3", tmp_path / "third.pdb", engine="haddock3", job_type="protein_protein", chains=["B"], score=None, unit="haddock"),
    ]
    config = {"subject": {"id": "synthetic"}, "receptor": {"path": "receptor.pdb"}, "poses": poses}
    document, _digests = analyze(config, base_dir=tmp_path)
    clusters = {row["pose_id"]: row["cluster_id"] for row in document["clusters"]}
    assert clusters == {"p1": "cluster_001", "p2": "cluster_001", "p3": "cluster_002"}
    assert {row["method"] for row in document["clusters"]} == {"fraction_common_contacts_ge_0.5"}
    stricter = analyze({**config, "fraction_common_contacts_cutoff": 0.51}, base_dir=tmp_path)[0]
    stricter_ids = {row["pose_id"]: row["cluster_id"] for row in stricter["clusters"]}
    assert stricter_ids == {"p1": "cluster_001", "p2": "cluster_002", "p3": "cluster_002"}


def test_sasa_stays_unevaluated_when_freesasa_is_on_path(tmp_path, capsys, monkeypatch):
    marker = tmp_path / "freesasa-was-called"
    binary = tmp_path / "bin" / "freesasa"
    binary.parent.mkdir()
    binary.write_text("#!/bin/sh\ntouch " + str(marker) + "\n", encoding="utf-8")
    binary.chmod(0o755)
    monkeypatch.setenv("PATH", str(binary.parent) + os.pathsep + os.environ.get("PATH", ""))
    out = tmp_path / "out"
    code, _stdout, err = _cli(capsys, ["run", "--config", str(_simple(tmp_path)), "--out", str(out)])
    assert code == 1 and err == ""
    assert not marker.exists()
    document = json.loads((out / "POSE_EVIDENCE.json").read_text(encoding="utf-8"))
    manifest = json.loads((out / "RUN_MANIFEST.json").read_text(encoding="utf-8"))
    assert document["buried_sasa"]["state"] == "unevaluated"
    assert manifest["optional_runtimes"]["freesasa"] == "not_invoked"
    assert manifest["missing_evidence"] == ["buried_sasa"]


def test_repeated_runs_are_byte_identical(tmp_path):
    _pdb(tmp_path / "receptor.pdb", [(1, "N", "ALA", "A", 1, 0, 0, 0, "N")])
    _pdb(tmp_path / "first.pdb", [(1, "C1", "LIG", "L", 1, 3, 0, 0, "C")])
    _pdb(tmp_path / "second.pdb", [(1, "C1", "LIG", "L", 1, 3.4, 0, 0, "C")])
    config = _config(tmp_path, {"subject": {"id": "synthetic"}, "receptor": {"path": "receptor.pdb"},
                                 "poses": [_pose("p2", tmp_path / "second.pdb"), _pose("p1", tmp_path / "first.pdb")]})
    first, second = tmp_path / "out-a", tmp_path / "out-b"
    assert main(["run", "--config", str(config), "--out", str(first)]) == 1
    assert main(["run", "--config", str(config), "--out", str(second)]) == 1
    assert _tree(first) == _tree(second)
    assert len(_tree(first)) == 5


def _invalid_pose_case(tmp_path: Path, kind: str) -> tuple[Path, str]:
    _pdb(tmp_path / "receptor.pdb", [(1, "CA", "ALA", "A", 1, 0, 0, 0, "C")])
    _pdb(tmp_path / "pose.pdb", [(1, "C1", "LIG", "L", 1, 3, 0, 0, "C")])
    _pdb(tmp_path / "other.pdb", [(1, "C2", "LIG", "L", 1, 3, 0, 0, "C")])
    _pdb(tmp_path / "duplicate.pdb", [
        (1, "C1", "LIG", "L", 1, 3, 0, 0, "C"),
        (2, "C1", "LIG", "L", 1, 4, 0, 0, "C"),
    ])
    _pdb(tmp_path / "hydrogens.pdb", [(1, "H1", "LIG", "L", 1, 3, 0, 0, "H")])
    document = {"subject": {"id": "synthetic"}, "receptor": {"path": "receptor.pdb", "sha256": sha256(tmp_path / "receptor.pdb")},
                "poses": [_pose("p1", tmp_path / "pose.pdb"), _pose("p2", tmp_path / "pose.pdb", score=-0.2)]}
    fragment = "pose_id values must be non-empty and unique"
    if kind == "checksum":
        document["receptor"]["sha256"] = "0" * 64
        fragment = "receptor checksum mismatch"
    elif kind == "engine":
        document["poses"][0]["engine"] = "glide"
        fragment = "unsupported pose engine"
    elif kind == "duplicate":
        document["poses"][1]["pose_id"] = "p1"
    elif kind == "empty_id":
        document["poses"][0]["pose_id"] = ""
    elif kind == "no_poses":
        document["poses"] = []
        fragment = "at least one pose is required"
    elif kind == "cutoff":
        document["contact_cutoff_A"] = 0
        fragment = "contact cutoff must be positive"
    elif kind == "negative_cutoff":
        document["contact_cutoff_A"] = -1
        fragment = "contact cutoff must be positive"
    elif kind == "score":
        document["poses"][0]["engine_score"] = "high"
        fragment = "engine_score must be numeric"
    elif kind == "job_type":
        document["poses"][0]["job_type"] = "nucleic_acid"
        fragment = "unsupported pose job_type"
    elif kind == "chain":
        document["poses"][0]["ligand_chains"] = ["Z"]
        fragment = "selected no ligand atoms"
    elif kind == "pose_shape":
        document["poses"] = ["p1"]
        fragment = "pose records must be objects"
    elif kind == "no_receptor":
        document = {"poses": document["poses"]}
        fragment = "receptor object is required"
    elif kind == "names":
        document["poses"][1] = _pose("p2", tmp_path / "other.pdb")
        fragment = "identical named heavy atoms"
    elif kind == "unique_atoms":
        document["poses"] = [_pose("p1", tmp_path / "duplicate.pdb"), _pose("p2", tmp_path / "duplicate.pdb")]
        fragment = "ligand atom identifiers are not unique"
    elif kind == "hydrogens":
        document["poses"] = [_pose("p1", tmp_path / "hydrogens.pdb")]
        fragment = "no heavy atoms"
    elif kind == "missing":
        document["receptor"]["path"] = "missing-receptor.pdb"
        fragment = "file not found"
    path = _config(tmp_path, document)
    if kind == "malformed":
        path.write_text("{", encoding="utf-8")
        fragment = "cannot read JSON"
    return path, fragment


@pytest.mark.parametrize("kind", [
    "checksum", "engine", "duplicate", "empty_id", "no_poses", "cutoff", "negative_cutoff", "score",
    "job_type", "chain", "pose_shape", "no_receptor", "names", "unique_atoms", "hydrogens", "missing", "malformed",
])
def test_invalid_pose_inputs_exit_two(tmp_path, capsys, kind):
    config, fragment = _invalid_pose_case(tmp_path, kind)
    out = tmp_path / "out"
    code, _stdout, err = _cli(capsys, ["run", "--config", str(config), "--out", str(out)])
    assert code == 2
    assert fragment in err
    assert not (out / "POSE_EVIDENCE.json").exists()


@pytest.mark.parametrize("payload", ["[]", "null", "\"config\""])
@pytest.mark.parametrize("command", ["validate", "run"])
def test_non_object_pose_config_is_invalid_input(tmp_path, capsys, command, payload):
    config = tmp_path / "config.json"
    config.write_text(payload + "\n", encoding="utf-8")
    out = tmp_path / "out"
    argv = ["validate", "--config", str(config)] if command == "validate" else ["run", "--config", str(config), "--out", str(out)]
    code, _stdout, err = _cli(capsys, argv)
    assert code == 2
    assert err.strip() == "pose-evidence: pose config must be a JSON object"
    assert "Traceback" not in err
    assert not (out / "POSE_EVIDENCE.json").exists()
