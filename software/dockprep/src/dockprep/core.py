from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
from typing import Any, Mapping

SCHEMA_VERSION = "1.0"
JOB_TYPES = {"small_molecule", "protein_protein", "protein_peptide"}
HADDOCK_CONTAINER_IMAGE = "protein-platform-haddock3:2026.7.0"
HADDOCK_CONTAINER_DIGEST = "sha256:c3a49243bb3c2cdf5be489096b75501ccb83dcc8cc2442532fa82afdb785c6b7"


class InputError(RuntimeError):
    pass


class ScientificBlock(RuntimeError):
    """Valid request that cannot yet produce defensible scientific evidence."""


def sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json_sha256(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read_json(path: str | Path) -> Any:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(f"cannot read JSON {path}: {exc}") from exc


def _resolve(base: Path, value: Any, label: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise InputError(f"{label}.path is required")
    path = Path(value)
    resolved = path.resolve() if path.is_absolute() else (base / path).resolve()
    if not resolved.is_file():
        raise InputError(f"{label} file does not exist: {value}")
    return resolved


def _copy_input(source: Path, out: Path, role: str) -> dict[str, Any]:
    destination_dir = out / "inputs"
    destination_dir.mkdir(parents=True, exist_ok=True)
    suffix = "".join(source.suffixes[-2:]) if source.name.endswith(".pdbqt.gz") else source.suffix
    destination = destination_dir / f"{role}{suffix.lower()}"
    shutil.copyfile(source, destination)
    return {"path": destination.relative_to(out).as_posix(), "sha256": sha256(destination), "format": source.suffix.lower().lstrip(".")}


def _approval(spec: Mapping[str, Any], receptor_sha: str, ligand_sha: str, base: Path) -> tuple[Path, dict[str, Any]]:
    approval_path = _resolve(base, spec.get("path"), "approval")
    record = read_json(approval_path)
    if not isinstance(record, Mapping):
        raise InputError("approval must be a JSON object")
    expected = {
        "status": "approved_for_docking",
        "receptor_sha256": receptor_sha,
        "ligand_sha256": ligand_sha,
    }
    mismatches = [field for field, value in expected.items() if str(record.get(field, "")) != value]
    if mismatches:
        raise ScientificBlock("checksum-bound docking approval does not match: " + ", ".join(mismatches))
    return approval_path, {"approval_id": str(record.get("approval_id", "")), "status": expected["status"],
                           "approved_job_type": str(record.get("job_type", ""))}


def _exact_sites(spec: Mapping[str, Any]) -> list[dict[str, Any]]:
    residues = spec.get("selected_site_residues", [])
    if not isinstance(residues, list):
        raise InputError("receptor.selected_site_residues must be a list")
    normalized = []
    for residue in residues:
        if not isinstance(residue, Mapping):
            raise InputError("selected-site residues must be objects")
        state = str(residue.get("mapping_state", "unresolved"))
        if state != "exact":
            raise ScientificBlock("selected-site residues must have exact residue mapping")
        try:
            position = int(residue["auth_seq_id"])
        except (KeyError, TypeError, ValueError) as exc:
            raise InputError("selected-site residue requires integer auth_seq_id") from exc
        normalized.append({"chain_id": str(residue.get("chain_id", "")), "auth_seq_id": position,
                           "insertion_code": str(residue.get("insertion_code", "")), "mapping_state": state})
    return sorted(normalized, key=lambda row: (row["chain_id"], row["auth_seq_id"], row["insertion_code"]))


def _chemistry(spec: Mapping[str, Any], ligand_path: Path) -> dict[str, Any]:
    if ligand_path.suffix.lower() in {".pdb", ".cif", ".mmcif", ".xyz"}:
        raise ScientificBlock("a coordinate-only small-molecule ligand is insufficient; declared connectivity and bond orders are required")
    chemistry = spec.get("chemistry")
    if not isinstance(chemistry, Mapping):
        raise ScientificBlock("small-molecule chemistry provenance is unresolved")
    required = ["identifier", "connectivity", "bond_orders", "charge", "protonation"]
    missing = [field for field in required if chemistry.get(field) in (None, "", [])]
    if missing:
        raise ScientificBlock("small-molecule chemistry is missing: " + ", ".join(missing))
    return {field: chemistry[field] for field in required} | {"provenance": str(chemistry.get("provenance", "user_declared"))}


def _search(spec: Mapping[str, Any], job_type: str) -> dict[str, Any]:
    search = spec.get("search", {})
    if not isinstance(search, Mapping):
        raise InputError("search must be an object")
    seed = int(search.get("seed", 0))
    if job_type == "small_molecule":
        center = search.get("center")
        size = search.get("size")
        if not (isinstance(center, list) and isinstance(size, list) and len(center) == len(size) == 3):
            raise InputError("small-molecule search requires three-value center and size arrays")
        if any(float(value) <= 0 for value in size):
            raise InputError("search-space sizes must be positive")
        return {"center": [round(float(v), 6) for v in center], "size": [round(float(v), 6) for v in size],
                "seed": seed, "exhaustiveness": int(search.get("exhaustiveness", 8)), "num_modes": int(search.get("num_modes", 9))}
    mode = str(search.get("mode", "ab_initio"))
    if mode not in {"experimental_restraints", "annotated_site_restraints", "ab_initio"}:
        raise InputError("unsupported macromolecular docking search mode")
    return {"mode": mode, "seed": seed, "restraints": list(search.get("restraints", []))}


def prepare_job(spec: Mapping[str, Any], *, base_dir: str | Path, out_dir: str | Path) -> dict[str, Any]:
    if not isinstance(spec, Mapping):
        raise InputError("job spec must be a JSON object")
    base = Path(base_dir).resolve()
    out = Path(out_dir)
    job_type = str(spec.get("job_type", ""))
    if job_type not in JOB_TYPES:
        raise InputError(f"job_type must be one of {sorted(JOB_TYPES)}")
    receptor_spec = spec.get("receptor")
    ligand_spec = spec.get("ligand")
    if not isinstance(receptor_spec, Mapping) or not isinstance(ligand_spec, Mapping):
        raise InputError("receptor and ligand objects are required")
    receptor_path = _resolve(base, receptor_spec.get("path"), "receptor")
    ligand_path = _resolve(base, ligand_spec.get("path"), "ligand")
    # Every gate runs before the first copy, so a block cannot leave inputs/.
    sites = _exact_sites(receptor_spec)
    assembly_id = str(receptor_spec.get("assembly_id", ""))
    conformational_state = str(receptor_spec.get("conformational_state", ""))
    if not assembly_id or not conformational_state:
        raise ScientificBlock("an exact assembly and explicit conformational state are required")
    identifier = str(ligand_spec.get("identifier", ""))
    if not identifier:
        raise InputError("ligand.identifier is required")
    if not bool(spec.get("source_terms", {}).get("derivatives_allowed", False)) or not bool(spec.get("source_terms", {}).get("docking_allowed", False)):
        raise ScientificBlock("source terms do not explicitly permit derivative preparation and docking")
    chemistry = _chemistry(ligand_spec, ligand_path) if job_type == "small_molecule" else None
    haddock_paths = None
    if job_type != "small_molecule":
        haddock = spec.get("haddock")
        if not isinstance(haddock, Mapping):
            raise ScientificBlock("macromolecular docking requires canonical HADDOCK3 workspace, receptor, and ligand manifests")
        if str(haddock.get("container_image", "")) != HADDOCK_CONTAINER_IMAGE or str(haddock.get("container_image_digest", "")) != HADDOCK_CONTAINER_DIGEST:
            raise ScientificBlock("HADDOCK3 container identity does not match the pinned canonical runtime")
        haddock_paths = {key: _resolve(base, haddock.get(key), f"haddock.{key}") for key in ("workspace_manifest", "receptor_manifest", "ligand_manifest")}
    approval_path, approval = _approval(spec.get("approval", {}), sha256(receptor_path), sha256(ligand_path), base)
    if approval["approved_job_type"] and approval["approved_job_type"] != job_type:
        raise ScientificBlock("approval job_type does not match the requested job")
    preparation = dict(spec.get("preparation", {}))
    for key in ("protonation_assumptions", "retained_waters", "metals", "cofactors", "repaired_atoms", "flexible_residues"):
        preparation.setdefault(key, [] if key != "protonation_assumptions" else "unresolved")
    search = _search(spec, job_type)
    restraints = list(spec.get("restraints", []))
    receptor = _copy_input(receptor_path, out, "receptor")
    ligand = _copy_input(ligand_path, out, "ligand")
    receptor.update({"assembly_id": assembly_id, "conformational_state": conformational_state, "selected_site_residues": sites})
    ligand["identifier"] = identifier
    if chemistry is None:
        chemistry = {"kind": job_type, "sequence_or_structure_checksum": ligand["sha256"]}
    canonical_haddock = None
    if haddock_paths is not None:
        copied = {key: _copy_input(path, out, f"haddock_{key}") for key, path in haddock_paths.items()}
        canonical_haddock = {"adapter": "canonical_haddock3_matrix_runner", "container_image": HADDOCK_CONTAINER_IMAGE,
                             "container_image_digest": HADDOCK_CONTAINER_DIGEST, "manifests": copied}
    approval = {**approval, "record": _copy_input(approval_path, out, "approval")}
    engine = "autodock_vina" if job_type == "small_molecule" else "haddock3"
    runtime_available = bool(shutil.which("vina")) if engine == "autodock_vina" else bool(shutil.which("haddock3") or shutil.which("docker"))
    content = {
        "job_type": job_type, "receptor_sha256": receptor["sha256"], "ligand_sha256": ligand["sha256"],
        "assembly_id": receptor["assembly_id"], "state": receptor["conformational_state"], "search": search,
        "approval_sha256": approval["record"]["sha256"],
    }
    return {
        "schema_version": SCHEMA_VERSION, "module_id": "docking_preparation", "job_id": "dock_" + _json_sha256(content)[:16],
        "job_type": job_type, "receptor": receptor, "ligand": ligand | {"chemistry": chemistry},
        "preparation": preparation, "search": search, "restraints": restraints,
        "approval": approval, "source_terms": {"derivatives_allowed": True, "docking_allowed": True},
        "execution": {"engine": engine, "status": "ready", "runtime_available_at_preparation": runtime_available,
                      "canonical_adapter": "refinement_workflow/run_haddock_matrix.py" if engine == "haddock3" else "meeko_vina",
                      **({"haddock": canonical_haddock} if canonical_haddock else {})},
        "missing_evidence": [] if runtime_available else [f"{engine}_runtime"],
        "limitations": ["Docking preparation and completion do not establish affinity or biological activity.",
                        "Engine scores and units remain method-specific.",
                        "No restraint is inferred from a favorable prior pose."],
    }


def write_job(out_dir: str | Path, document: Mapping[str, Any], *, config_sha256: str) -> list[Path]:
    out = Path(out_dir)
    job = out / "DOCKING_JOB.json"
    run = out / "RUN_MANIFEST.json"
    job.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    run.write_text(json.dumps({
        "schema_version": SCHEMA_VERSION, "module_id": "docking_preparation", "version": "0.1.0",
        "input_sha256": {"job_spec": config_sha256, "receptor": document["receptor"]["sha256"], "ligand": document["ligand"]["sha256"],
                          "approval": document["approval"]["record"]["sha256"]},
        "parameters": {"job_type": document["job_type"], "search": document["search"], "preparation": document["preparation"]},
        "runtime_versions": {"python": platform.python_version()},
        "optional_runtimes": {document["execution"]["engine"]: "available_not_invoked" if document["execution"]["runtime_available_at_preparation"] else "not_available"},
        "outputs": [job.name], "missing_evidence": document["missing_evidence"],
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return [job, run]


def _runtime_version(binary: str) -> str:
    result = subprocess.run([binary, "--version"], capture_output=True, text=True, check=False, timeout=30)
    return ((result.stdout or result.stderr).strip().splitlines() or ["version_unreported"])[0]


def execute_job(job_path: str | Path, out_dir: str | Path) -> tuple[dict[str, Any], int]:
    job_file = Path(job_path).resolve()
    job = read_json(job_file)
    if not isinstance(job, Mapping):
        raise InputError("job bundle must be a JSON object")
    execution = job.get("execution", {})
    if not isinstance(execution, Mapping):
        raise InputError("job bundle is not a ready dockprep job")
    if job.get("module_id") != "docking_preparation" or execution.get("status") != "ready":
        raise InputError("job bundle is not a ready dockprep job")
    receptor_record = job.get("receptor")
    if not isinstance(receptor_record, Mapping) or not isinstance(receptor_record.get("path"), str) or not isinstance(receptor_record.get("sha256"), str):
        raise InputError("job bundle is missing receptor")
    base = job_file.parent
    receptor = (base / receptor_record["path"]).resolve()
    ligand = (base / job["ligand"]["path"]).resolve()
    if not receptor.is_relative_to(base) or not ligand.is_relative_to(base):
        raise InputError("job input path escapes the job bundle")
    if sha256(receptor) != job["receptor"]["sha256"] or sha256(ligand) != job["ligand"]["sha256"]:
        raise InputError("job input checksum mismatch")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    engine = job["execution"]["engine"]
    logs: list[str] = []
    output_files: list[str] = []
    status = "blocked"
    runtime: dict[str, Any] = {"engine": engine, "version": "not_available"}
    missing: list[str] = []
    returncode = 1
    if engine == "autodock_vina":
        vina = shutil.which("vina")
        if vina is None:
            missing.append("autodock_vina_runtime")
        else:
            receptor_pdbqt = receptor
            ligand_pdbqt = ligand
            if receptor.suffix.lower() != ".pdbqt":
                prepare_receptor = shutil.which("mk_prepare_receptor.py") or shutil.which("mk_prepare_receptor")
                if prepare_receptor is None:
                    missing.append("meeko_receptor_preparation_runtime")
                else:
                    receptor_pdbqt = out / "receptor.pdbqt"
                    result = subprocess.run([prepare_receptor, "-i", str(receptor), "-o", str(receptor_pdbqt)], capture_output=True, text=True, check=False)
                    (out / "meeko_receptor.log").write_text(result.stdout + result.stderr, encoding="utf-8")
                    logs.append("meeko_receptor.log")
                    if result.returncode:
                        missing.append("successful_receptor_preparation")
            if ligand.suffix.lower() != ".pdbqt":
                prepare_ligand = shutil.which("mk_prepare_ligand.py") or shutil.which("mk_prepare_ligand")
                if prepare_ligand is None:
                    missing.append("meeko_ligand_preparation_runtime")
                else:
                    ligand_pdbqt = out / "ligand.pdbqt"
                    result = subprocess.run([prepare_ligand, "-i", str(ligand), "-o", str(ligand_pdbqt)], capture_output=True, text=True, check=False)
                    (out / "meeko_ligand.log").write_text(result.stdout + result.stderr, encoding="utf-8")
                    logs.append("meeko_ligand.log")
                    if result.returncode:
                        missing.append("successful_ligand_preparation")
            if not missing:
                search = job["search"]
                pose_path = out / "VINA_POSES.pdbqt"
                log_path = out / "vina.log"
                command = [vina, "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_pdbqt),
                           "--center_x", str(search["center"][0]), "--center_y", str(search["center"][1]), "--center_z", str(search["center"][2]),
                           "--size_x", str(search["size"][0]), "--size_y", str(search["size"][1]), "--size_z", str(search["size"][2]),
                           "--seed", str(search["seed"]), "--exhaustiveness", str(search["exhaustiveness"]), "--num_modes", str(search["num_modes"]),
                           "--out", str(pose_path), "--log", str(log_path)]
                result = subprocess.run(command, cwd=out, capture_output=True, text=True, check=False)
                (out / "vina_process.log").write_text(result.stdout + result.stderr, encoding="utf-8")
                logs.append("vina_process.log")
                if log_path.exists(): logs.append(log_path.name)
                if result.returncode == 0 and pose_path.is_file():
                    status, returncode = "complete", 0
                    output_files.append(pose_path.name)
                else:
                    status, returncode = "failed", 1
                    missing.append("successful_vina_execution")
                runtime = {"engine": engine, "version": _runtime_version(vina)}
    else:
        runner = os.environ.get("YAUVI_HADDOCK_MATRIX_RUNNER", "")
        runner_path = Path(runner).resolve() if runner else None
        haddock_record = job.get("execution", {}).get("haddock", {})
        workspace_record = (haddock_record.get("manifests", {}) or {}).get("workspace_manifest", {})
        workspace_manifest = (base / workspace_record.get("path", "")).resolve() if workspace_record.get("path") else None
        if runner_path is None or not runner_path.is_file():
            missing.append("canonical_haddock3_matrix_runner")
        elif workspace_manifest is None or not workspace_manifest.is_file():
            missing.append("canonical_haddock3_workspace_manifest")
        else:
            python = shutil.which("python3")
            command = [python or "python3", str(runner_path), "--workspace-manifest", str(workspace_manifest),
                       "--container-image", HADDOCK_CONTAINER_IMAGE]
            result = subprocess.run(command, cwd=base, capture_output=True, text=True, check=False)
            (out / "haddock_runner.log").write_text(result.stdout + result.stderr, encoding="utf-8")
            logs.append("haddock_runner.log")
            status, returncode = ("complete", 0) if result.returncode == 0 else ("failed", 1)
            if result.returncode:
                missing.append("successful_haddock3_execution")
            runtime = {"engine": engine, "version": "canonical_runner", "runner_sha256": sha256(runner_path)}
    manifest = {
        "schema_version": SCHEMA_VERSION, "module_id": "docking_preparation", "job_id": job["job_id"], "engine": engine,
        "status": status, "input_sha256": {"job": sha256(job_file), "receptor": job["receptor"]["sha256"], "ligand": job["ligand"]["sha256"]},
        "runtime": runtime, "outputs": sorted(output_files + logs), "missing_evidence": sorted(set(missing)),
        "limitations": ["A docking score is not affinity, engagement, activity, efficacy, or mechanism."],
    }
    (out / "DOCKING_EXECUTION.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest, returncode
