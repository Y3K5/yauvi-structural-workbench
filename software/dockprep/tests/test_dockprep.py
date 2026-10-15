from __future__ import annotations

from importlib.resources import files
import json
from pathlib import Path

import pytest

from dockprep.cli import describe, main
from dockprep.core import HADDOCK_CONTAINER_DIGEST, HADDOCK_CONTAINER_IMAGE, ScientificBlock, prepare_job, sha256


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


EXECUTION_CLAIM = "A docking score is not affinity, engagement, activity, efficacy, or mechanism."
PREPARATION_CLAIM = "Docking preparation and completion do not establish affinity or biological activity."


def _spec(tmp_path: Path, document: dict) -> Path:
    path = tmp_path / "job.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return path


def _cli(capsys, argv: list[str]) -> tuple[int, str, str]:
    code = main(argv)
    captured = capsys.readouterr()
    return code, captured.out, captured.err


def _hide_runtimes(monkeypatch) -> None:
    monkeypatch.setattr("dockprep.core.shutil.which", lambda _name: None)
    monkeypatch.delenv("YAUVI_HADDOCK_MATRIX_RUNNER", raising=False)


def _script(path: Path, body: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/usr/bin/env python3\n" + body, encoding="utf-8")
    path.chmod(0o755)
    return str(path)


def _tree(root: Path) -> dict[str, bytes]:
    return {path.relative_to(root).as_posix(): path.read_bytes() for path in sorted(root.rglob("*")) if path.is_file()}


def _macromolecular(tmp_path: Path) -> dict:
    receptor = tmp_path / "receptor.pdb"
    ligand = tmp_path / "partner.pdb"
    receptor.write_text("ATOM\n", encoding="utf-8")
    ligand.write_text("ATOM\n", encoding="utf-8")
    manifests = {}
    for key, name in (
        ("workspace_manifest", "workspace.json"),
        ("receptor_manifest", "receptor_manifest.json"),
        ("ligand_manifest", "ligand_manifest.json"),
    ):
        (tmp_path / name).write_text("{}\n", encoding="utf-8")
        manifests[key] = name
    approval = tmp_path / "approval.json"
    approval.write_text(json.dumps({
        "approval_id": "approval-macro", "status": "approved_for_docking", "job_type": "protein_protein",
        "receptor_sha256": sha256(receptor), "ligand_sha256": sha256(ligand),
    }), encoding="utf-8")
    return {
        "job_type": "protein_protein",
        "receptor": {"path": receptor.name, "assembly_id": "1", "conformational_state": "observed",
                     "selected_site_residues": [{"chain_id": "A", "auth_seq_id": 1, "mapping_state": "exact"}]},
        "ligand": {"path": ligand.name, "identifier": "synthetic-partner"},
        "source_terms": {"derivatives_allowed": True, "docking_allowed": True},
        "approval": {"path": approval.name},
        "haddock": {"container_image": HADDOCK_CONTAINER_IMAGE, "container_image_digest": HADDOCK_CONTAINER_DIGEST, **manifests},
    }


def _non_pdbqt_small_molecule(tmp_path: Path) -> dict:
    receptor = tmp_path / "receptor.pdb"
    ligand = tmp_path / "ligand.sdf"
    receptor.write_text("ATOM\n", encoding="utf-8")
    ligand.write_text("synthetic\n", encoding="utf-8")
    approval = tmp_path / "approval.json"
    approval.write_text(json.dumps({
        "approval_id": "approval-sdf", "status": "approved_for_docking", "job_type": "small_molecule",
        "receptor_sha256": sha256(receptor), "ligand_sha256": sha256(ligand),
    }), encoding="utf-8")
    return {
        "job_type": "small_molecule",
        "receptor": {"path": receptor.name, "assembly_id": "1", "conformational_state": "observed",
                     "selected_site_residues": [{"chain_id": "A", "auth_seq_id": 10, "mapping_state": "exact"}]},
        "ligand": {"path": ligand.name, "identifier": "synthetic-ligand", "chemistry": {
            "identifier": "synthetic-ligand", "connectivity": "declared", "bond_orders": "declared",
            "charge": 0, "protonation": "declared"}},
        "source_terms": {"derivatives_allowed": True, "docking_allowed": True},
        "approval": {"path": approval.name},
        "search": {"center": [0, 0, 0], "size": [20, 20, 20], "seed": 7},
    }


def test_describe_lists_the_execution_command_and_claim_limits(capsys):
    code, out, err = _cli(capsys, ["describe"])
    assert code == 0 and err == ""
    payload = json.loads(out)
    assert payload == describe()
    assert payload["commands"] == ["describe", "validate", "fetch", "run", "execute"]
    assert payload["deterministic_output"] is True
    assert "not affinity or mechanism" in payload["limitations"][1]


def test_fetch_prints_the_packaged_source_plan(capsys):
    code, out, err = _cli(capsys, ["fetch", "--plan"])
    assert code == 0 and err == ""
    assert out == (files("dockprep") / "sources.yaml").read_text()
    assert "autodock_vina" in out and "haddock3" in out


def test_validate_reports_missing_vina_and_writes_no_job(tmp_path, capsys, monkeypatch):
    _hide_runtimes(monkeypatch)
    config = _spec(tmp_path, _inputs(tmp_path))
    code, out, err = _cli(capsys, ["validate", "--config", str(config)])
    payload = json.loads(out)
    assert code == 1 and err == ""
    assert payload["valid"] is True
    assert payload["missing_evidence"] == ["autodock_vina_runtime"]
    assert payload["job_id"].startswith("dock_")
    assert not (tmp_path / "DOCKING_JOB.json").exists()


def test_validate_is_clean_when_vina_is_reported_present(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr("dockprep.core.shutil.which", lambda name: "/mock/vina" if name == "vina" else None)
    config = _spec(tmp_path, _inputs(tmp_path))
    code, out, err = _cli(capsys, ["validate", "--config", str(config)])
    payload = json.loads(out)
    assert code == 0 and err == ""
    assert payload["missing_evidence"] == []


def test_run_records_the_preparation_ceiling_when_vina_is_absent(tmp_path, capsys, monkeypatch):
    _hide_runtimes(monkeypatch)
    config = _spec(tmp_path, _inputs(tmp_path))
    out = tmp_path / "out"
    code, _stdout, err = _cli(capsys, ["run", "--config", str(config), "--out", str(out)])
    assert code == 1 and err == ""
    job = json.loads((out / "DOCKING_JOB.json").read_text(encoding="utf-8"))
    manifest = json.loads((out / "RUN_MANIFEST.json").read_text(encoding="utf-8"))
    assert job["missing_evidence"] == ["autodock_vina_runtime"]
    assert job["execution"]["status"] == "ready"
    assert job["execution"]["runtime_available_at_preparation"] is False
    assert PREPARATION_CLAIM in job["limitations"]
    assert manifest["missing_evidence"] == ["autodock_vina_runtime"]
    assert manifest["optional_runtimes"] == {"autodock_vina": "not_available"}
    assert job["receptor"]["path"].startswith("inputs/")


@pytest.mark.parametrize("kind,fragment", [
    ("job_type", "job_type must be one of"),
    ("size", "search-space sizes must be positive"),
    ("center", "three-value center and size"),
    ("search", "search must be an object"),
    ("identifier", "ligand.identifier is required"),
    ("missing_approval", "does not exist"),
    ("missing_receptor", "does not exist"),
    ("malformed", "cannot read JSON"),
    ("empty_object", "job_type must be one of"),
])
def test_invalid_job_specs_exit_two(tmp_path, capsys, kind, fragment):
    if kind == "malformed":
        config = tmp_path / "job.json"
        config.write_text("{", encoding="utf-8")
    elif kind == "empty_object":
        config = _spec(tmp_path, {})
    else:
        spec = _inputs(tmp_path)
        if kind == "job_type":
            spec["job_type"] = "not-a-job"
        elif kind == "size":
            spec["search"]["size"] = [0, 20, 20]
        elif kind == "center":
            spec["search"]["center"] = [0, 0]
        elif kind == "search":
            spec["search"] = "box"
        elif kind == "identifier":
            spec["ligand"]["identifier"] = ""
        elif kind == "missing_approval":
            spec["approval"]["path"] = "no-such-approval.json"
        else:
            (tmp_path / "receptor.pdbqt").unlink()
        config = _spec(tmp_path, spec)
    out = tmp_path / "out"
    code, _stdout, err = _cli(capsys, ["run", "--config", str(config), "--out", str(out)])
    assert code == 2
    assert fragment in err
    assert not (out / "DOCKING_JOB.json").exists()


@pytest.mark.parametrize("kind,fragment", [
    ("checksum", "checksum-bound docking approval"),
    ("source_terms", "source terms"),
    ("mapping", "exact residue mapping"),
    ("digest", "pinned canonical runtime"),
    ("no_haddock", "canonical HADDOCK3"),
    ("approval_job_type", "approval job_type"),
])
def test_scientific_blocks_exit_one_without_writing_a_job(tmp_path, capsys, kind, fragment):
    if kind in {"digest", "no_haddock"}:
        spec = _macromolecular(tmp_path)
        if kind == "digest":
            spec["haddock"]["container_image_digest"] = "sha256:" + ("ab" * 32)
        else:
            del spec["haddock"]
    else:
        spec = _inputs(tmp_path)
        if kind in {"checksum", "approval_job_type"}:
            approval = json.loads((tmp_path / "approval.json").read_text(encoding="utf-8"))
            if kind == "checksum":
                approval["ligand_sha256"] = "f" * 64
            else:
                approval["job_type"] = "protein_protein"
            (tmp_path / "approval.json").write_text(json.dumps(approval), encoding="utf-8")
        elif kind == "source_terms":
            spec["source_terms"] = {"derivatives_allowed": False, "docking_allowed": False}
        else:
            spec["receptor"]["selected_site_residues"][0]["mapping_state"] = "ambiguous"
    config = _spec(tmp_path, spec)
    out = tmp_path / "out"
    code, _stdout, err = _cli(capsys, ["run", "--config", str(config), "--out", str(out)])
    assert code == 1
    assert "scientifically blocked" in err
    assert fragment in err
    assert not (out / "DOCKING_JOB.json").exists()
    assert not (out / "inputs").exists()
    assert not any(out.rglob("*")) if out.exists() else True


def test_execute_refuses_before_the_job_bundle_exists(tmp_path, capsys):
    missing = tmp_path / "DOCKING_JOB.json"
    out = tmp_path / "execution"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(missing), "--out", str(out)])
    assert code == 2
    assert "cannot read JSON" in err
    assert not out.exists()


def test_execute_refuses_a_bundle_that_is_not_ready(tmp_path, capsys):
    job = tmp_path / "DOCKING_JOB.json"
    job.write_text(json.dumps({"module_id": "docking_preparation", "execution": {"status": "blocked"}}), encoding="utf-8")
    out = tmp_path / "execution"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(job), "--out", str(out)])
    assert code == 2
    assert "not a ready dockprep job" in err
    assert not (out / "DOCKING_EXECUTION.json").exists()


@pytest.mark.parametrize("payload", ["[]", "null", "\"spec\""])
@pytest.mark.parametrize("command", ["validate", "run"])
def test_non_object_job_spec_is_invalid_input(tmp_path, capsys, command, payload):
    config = tmp_path / "job.json"
    config.write_text(payload + "\n", encoding="utf-8")
    out = tmp_path / "out"
    argv = ["validate", "--config", str(config)] if command == "validate" else ["run", "--config", str(config), "--out", str(out)]
    code, _stdout, err = _cli(capsys, argv)
    assert code == 2
    assert err.strip() == "dockprep: job spec must be a JSON object"
    assert "Traceback" not in err
    assert not (out / "DOCKING_JOB.json").exists()


@pytest.mark.parametrize("payload", ["[]", "null"])
def test_execute_rejects_a_non_object_bundle(tmp_path, capsys, payload):
    job = tmp_path / "DOCKING_JOB.json"
    job.write_text(payload + "\n", encoding="utf-8")
    out = tmp_path / "execution"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(job), "--out", str(out)])
    assert code == 2
    assert err.strip() == "dockprep: job bundle must be a JSON object"
    assert "Traceback" not in err
    assert not (out / "DOCKING_EXECUTION.json").exists()


def test_execute_rejects_a_ready_job_missing_receptor(tmp_path, capsys):
    job = tmp_path / "DOCKING_JOB.json"
    job.write_text(json.dumps({"module_id": "docking_preparation", "execution": {"status": "ready", "engine": "autodock_vina"}}), encoding="utf-8")
    out = tmp_path / "execution"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(job), "--out", str(out)])
    assert code == 2
    assert err.strip() == "dockprep: job bundle is missing receptor"
    assert "Traceback" not in err
    assert not (out / "DOCKING_EXECUTION.json").exists()


def test_execute_rejects_checksum_mismatch_and_path_escape(tmp_path, capsys, monkeypatch):
    _hide_runtimes(monkeypatch)
    config = _spec(tmp_path, _inputs(tmp_path))
    prepared = tmp_path / "prepared"
    assert main(["run", "--config", str(config), "--out", str(prepared)]) == 1
    capsys.readouterr()
    job_path = prepared / "DOCKING_JOB.json"
    (prepared / "inputs" / "receptor.pdbqt").write_text("CHANGED\n", encoding="utf-8")
    mismatched = tmp_path / "mismatched"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(job_path), "--out", str(mismatched)])
    assert code == 2
    assert "checksum mismatch" in err
    assert not (mismatched / "DOCKING_EXECUTION.json").exists()

    job = json.loads(job_path.read_text(encoding="utf-8"))
    job["receptor"]["path"] = "../outside.pdbqt"
    escaped_job = tmp_path / "escaped.json"
    escaped_job.write_text(json.dumps(job), encoding="utf-8")
    # The tampered receptor still has to be a file inside this bundle's parent
    # only after the path check. Escape is rejected before the bytes are read.
    (prepared / "inputs" / "receptor.pdbqt").write_bytes((tmp_path / "receptor.pdbqt").read_bytes())
    escaped_out = tmp_path / "escaped"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(escaped_job), "--out", str(escaped_out)])
    assert code == 2
    assert "escapes the job bundle" in err
    assert not (escaped_out / "DOCKING_EXECUTION.json").exists()


def test_absent_vina_blocks_execution_without_crashing(tmp_path, capsys, monkeypatch):
    _hide_runtimes(monkeypatch)
    config = _spec(tmp_path, _inputs(tmp_path))
    prepared = tmp_path / "prepared"
    assert main(["run", "--config", str(config), "--out", str(prepared)]) == 1
    capsys.readouterr()
    out = tmp_path / "execution"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(prepared / "DOCKING_JOB.json"), "--out", str(out)])
    assert code == 1 and err == ""
    execution = json.loads((out / "DOCKING_EXECUTION.json").read_text(encoding="utf-8"))
    assert execution["status"] == "blocked"
    assert execution["missing_evidence"] == ["autodock_vina_runtime"]
    assert execution["limitations"] == [EXECUTION_CLAIM]
    assert execution["outputs"] == []


def test_absent_meeko_blocks_non_pdbqt_inputs_without_crashing(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr("dockprep.core.shutil.which", lambda name: "/mock/vina" if name == "vina" else None)
    config = _spec(tmp_path, _non_pdbqt_small_molecule(tmp_path))
    prepared = tmp_path / "prepared"
    assert main(["run", "--config", str(config), "--out", str(prepared)]) == 0
    capsys.readouterr()
    out = tmp_path / "execution"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(prepared / "DOCKING_JOB.json"), "--out", str(out)])
    assert code == 1 and err == ""
    execution = json.loads((out / "DOCKING_EXECUTION.json").read_text(encoding="utf-8"))
    assert execution["status"] == "blocked"
    assert execution["missing_evidence"] == [
        "meeko_ligand_preparation_runtime", "meeko_receptor_preparation_runtime",
    ]
    assert execution["limitations"] == [EXECUTION_CLAIM]
    assert not (out / "VINA_POSES.pdbqt").exists()


def test_absent_haddock_runner_is_named_missing_evidence(tmp_path, capsys, monkeypatch):
    _hide_runtimes(monkeypatch)
    config = _spec(tmp_path, _macromolecular(tmp_path))
    prepared = tmp_path / "prepared"
    code, _stdout, err = _cli(capsys, ["run", "--config", str(config), "--out", str(prepared)])
    assert code == 1 and err == ""
    job = json.loads((prepared / "DOCKING_JOB.json").read_text(encoding="utf-8"))
    assert job["missing_evidence"] == ["haddock3_runtime"]
    assert job["execution"]["engine"] == "haddock3"
    assert job["execution"]["haddock"]["container_image_digest"] == HADDOCK_CONTAINER_DIGEST
    out = tmp_path / "execution"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(prepared / "DOCKING_JOB.json"), "--out", str(out)])
    assert code == 1 and err == ""
    execution = json.loads((out / "DOCKING_EXECUTION.json").read_text(encoding="utf-8"))
    assert execution["status"] == "blocked"
    assert execution["missing_evidence"] == ["canonical_haddock3_matrix_runner"]
    assert execution["limitations"] == [EXECUTION_CLAIM]


def test_mocked_haddock_runner_completes_without_a_container(tmp_path, capsys, monkeypatch):
    _hide_runtimes(monkeypatch)
    runner = tmp_path / "run_haddock_matrix.py"
    runner.write_text("raise SystemExit(0)\n", encoding="utf-8")
    monkeypatch.setenv("YAUVI_HADDOCK_MATRIX_RUNNER", str(runner))
    config = _spec(tmp_path, _macromolecular(tmp_path))
    prepared = tmp_path / "prepared"
    assert main(["run", "--config", str(config), "--out", str(prepared)]) == 1
    capsys.readouterr()
    out = tmp_path / "execution"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(prepared / "DOCKING_JOB.json"), "--out", str(out)])
    assert code == 0 and err == ""
    execution = json.loads((out / "DOCKING_EXECUTION.json").read_text(encoding="utf-8"))
    assert execution["status"] == "complete"
    assert execution["missing_evidence"] == []
    assert "haddock_runner.log" in execution["outputs"]
    assert execution["limitations"] == [EXECUTION_CLAIM]


def test_failing_mocked_vina_is_reported_rather_than_raised(tmp_path, capsys, monkeypatch):
    vina = _script(tmp_path / "bin" / "vina", "raise SystemExit(3)\n")
    monkeypatch.setattr("dockprep.core.shutil.which", lambda name: vina if name == "vina" else None)
    config = _spec(tmp_path, _inputs(tmp_path))
    prepared = tmp_path / "prepared"
    assert main(["run", "--config", str(config), "--out", str(prepared)]) == 0
    capsys.readouterr()
    out = tmp_path / "execution"
    code, _stdout, err = _cli(capsys, ["execute", "--job", str(prepared / "DOCKING_JOB.json"), "--out", str(out)])
    assert code == 1 and err == ""
    execution = json.loads((out / "DOCKING_EXECUTION.json").read_text(encoding="utf-8"))
    assert execution["status"] == "failed"
    assert execution["missing_evidence"] == ["successful_vina_execution"]
    assert execution["limitations"] == [EXECUTION_CLAIM]


def test_blocked_runs_are_byte_identical(tmp_path, monkeypatch):
    _hide_runtimes(monkeypatch)
    config = _spec(tmp_path, _inputs(tmp_path))
    first, second = tmp_path / "first", tmp_path / "second"
    assert main(["run", "--config", str(config), "--out", str(first)]) == 1
    assert main(["run", "--config", str(config), "--out", str(second)]) == 1
    assert _tree(first) == _tree(second)
    executed = [tmp_path / "exec-first", tmp_path / "exec-second"]
    assert main(["execute", "--job", str(first / "DOCKING_JOB.json"), "--out", str(executed[0])]) == 1
    assert main(["execute", "--job", str(first / "DOCKING_JOB.json"), "--out", str(executed[1])]) == 1
    assert _tree(executed[0]) == _tree(executed[1])
    assert b"autodock_vina_runtime" in (executed[0] / "DOCKING_EXECUTION.json").read_bytes()


def test_mocked_vina_and_meeko_complete_and_repeat(tmp_path, capsys, monkeypatch):
    vina = _script(tmp_path / "bin" / "vina", """
import sys
from pathlib import Path
args = sys.argv[1:]
if "version" in args or "--version" in args:
    print("vina mock-1")
    raise SystemExit(0)
Path(args[args.index("--out") + 1]).write_text("REMARK synthetic pose\\n", encoding="utf-8")
raise SystemExit(0)
""")
    meeko = _script(tmp_path / "bin" / "mk_prepare_receptor.py", """
import sys
from pathlib import Path
args = sys.argv[1:]
Path(args[args.index("-o") + 1]).write_text("ATOM\\n", encoding="utf-8")
raise SystemExit(0)
""")
    ligand_prep = _script(tmp_path / "bin" / "mk_prepare_ligand.py", """
import sys
from pathlib import Path
args = sys.argv[1:]
Path(args[args.index("-o") + 1]).write_text("ATOM\\n", encoding="utf-8")
raise SystemExit(0)
""")

    def which(name):
        return {"vina": vina, "mk_prepare_receptor.py": meeko, "mk_prepare_ligand.py": ligand_prep}.get(name)

    monkeypatch.setattr("dockprep.core.shutil.which", which)
    config = _spec(tmp_path, _non_pdbqt_small_molecule(tmp_path))
    prepared = tmp_path / "prepared"
    assert main(["run", "--config", str(config), "--out", str(prepared)]) == 0
    capsys.readouterr()
    outs = [tmp_path / "exec-a", tmp_path / "exec-b"]
    assert main(["execute", "--job", str(prepared / "DOCKING_JOB.json"), "--out", str(outs[0])]) == 0
    assert main(["execute", "--job", str(prepared / "DOCKING_JOB.json"), "--out", str(outs[1])]) == 0
    assert _tree(outs[0]) == _tree(outs[1])
    execution = json.loads((outs[0] / "DOCKING_EXECUTION.json").read_text(encoding="utf-8"))
    assert execution["status"] == "complete"
    assert execution["missing_evidence"] == []
    assert execution["runtime"]["version"] == "vina mock-1"
    assert "VINA_POSES.pdbqt" in execution["outputs"]
    assert execution["limitations"] == [EXECUTION_CLAIM]
    assert (outs[0] / "VINA_POSES.pdbqt").read_text(encoding="utf-8") == "REMARK synthetic pose\n"
