from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import platform
from typing import Any, Iterable, Mapping

import numpy as np

SCHEMA_VERSION = "1.0"
METALS = {"ZN", "MG", "MN", "FE", "CO", "NI", "CU", "CA"}
POSITIVE = {"ARG", "LYS", "HIS"}
NEGATIVE = {"ASP", "GLU"}


class InputError(RuntimeError):
    pass


def sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path: str | Path) -> Any:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(f"cannot read JSON {path}: {exc}") from exc


def _resolve(base: Path, value: Any, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise InputError(f"{label}.path is required")
    path = Path(value)
    path = path.resolve() if path.is_absolute() else (base / path).resolve()
    if not path.is_file():
        raise InputError(f"{label} file not found: {value}")
    return path


def parse_atoms(path: str | Path) -> list[dict[str, Any]]:
    atoms = []
    for index, line in enumerate(Path(path).read_text(errors="replace").splitlines(), start=1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            residue_number = int(line[22:26])
        except ValueError as exc:
            raise InputError(f"invalid coordinate record at {Path(path).name}:{index}") from exc
        element = line[76:78].strip().upper() or line[12:16].strip()[0].upper()
        if element == "H":
            continue
        atoms.append({"record": line[:6].strip(), "atom_name": line[12:16].strip(), "resname": line[17:20].strip(),
                      "chain": line[21:22].strip() or "_", "resseq": residue_number, "icode": line[26:27].strip(),
                      "xyz": xyz, "element": element})
    if not atoms:
        raise InputError(f"no heavy atoms found in {path}")
    return atoms


def _distance(left: Mapping[str, Any], right: Mapping[str, Any]) -> float:
    return float(np.linalg.norm(np.asarray(left["xyz"], dtype=float) - np.asarray(right["xyz"], dtype=float)))


def _residue_key(atom: Mapping[str, Any]) -> tuple[str, int, str]:
    return str(atom["chain"]), int(atom["resseq"]), str(atom["icode"])


def _pose_contacts(receptor_atoms: list[dict[str, Any]], ligand_atoms: list[dict[str, Any]], pose_id: str, cutoff: float) -> list[dict[str, Any]]:
    residue_pairs: dict[tuple[tuple[str, int, str], tuple[str, int, str]], dict[str, Any]] = {}
    for receptor_atom in receptor_atoms:
        for ligand_atom in ligand_atoms:
            distance = _distance(receptor_atom, ligand_atom)
            if distance > cutoff:
                continue
            key = (_residue_key(receptor_atom), _residue_key(ligand_atom))
            prior = residue_pairs.get(key)
            if prior is None or distance < prior["minimum_distance_A"]:
                salt = ((receptor_atom["resname"] in POSITIVE and ligand_atom["resname"] in NEGATIVE) or
                        (receptor_atom["resname"] in NEGATIVE and ligand_atom["resname"] in POSITIVE)) and distance <= 4.0
                hbond = receptor_atom["element"] in {"N", "O", "S"} and ligand_atom["element"] in {"N", "O", "S"} and distance <= 3.5
                metal = (receptor_atom["element"] in METALS or ligand_atom["element"] in METALS) and distance <= 3.0
                residue_pairs[key] = {
                    "pose_id": pose_id, "receptor_chain": key[0][0], "receptor_auth_seq_id": key[0][1],
                    "receptor_insertion_code": key[0][2], "receptor_resname": receptor_atom["resname"],
                    "ligand_chain": key[1][0], "ligand_auth_seq_id": key[1][1], "ligand_insertion_code": key[1][2],
                    "ligand_resname": ligand_atom["resname"], "minimum_distance_A": round(distance, 6),
                    "contact_state": "clash" if distance < 2.0 else "contact", "hydrogen_bond_proxy": hbond,
                    "salt_bridge_proxy": salt, "metal_coordination_proxy": metal,
                }
    return sorted(residue_pairs.values(), key=lambda row: (row["pose_id"], row["receptor_chain"], row["receptor_auth_seq_id"], row["ligand_chain"], row["ligand_auth_seq_id"]))


def _ligand_rmsd(left: Mapping[str, Any], right: Mapping[str, Any]) -> float:
    def keyed(pose: Mapping[str, Any]) -> dict[tuple[str, int, str, str, str], list[float]]:
        result = {}
        for atom in pose["ligand_atoms"]:
            key = (atom["chain"], atom["resseq"], atom["icode"], atom["resname"], atom["atom_name"])
            if key in result:
                raise InputError("ligand atom identifiers are not unique, so RMSD is unresolved")
            result[key] = atom["xyz"]
        return result
    first, second = keyed(left), keyed(right)
    common = sorted(set(first) & set(second))
    if not common or set(first) != set(second):
        raise InputError("small-molecule pose RMSD requires identical named heavy atoms")
    delta = np.asarray([first[key] for key in common]) - np.asarray([second[key] for key in common])
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _contact_set(pose: Mapping[str, Any]) -> set[tuple[str, int, str, int]]:
    return {(row["receptor_chain"], row["receptor_auth_seq_id"], row["ligand_chain"], row["ligand_auth_seq_id"]) for row in pose["contacts"]}


def _cluster(poses: list[dict[str, Any]], rmsd_cutoff: float, fcc_cutoff: float) -> list[dict[str, Any]]:
    assigned: dict[str, str] = {}
    clusters: list[list[dict[str, Any]]] = []
    for pose in sorted(poses, key=lambda value: value["pose_id"]):
        placed = False
        for index, members in enumerate(clusters, start=1):
            representative = members[0]
            if pose["job_type"] == "small_molecule" and representative["job_type"] == "small_molecule":
                similar = _ligand_rmsd(pose, representative) <= rmsd_cutoff
                method = f"ligand_heavy_atom_rmsd_le_{rmsd_cutoff}_A"
            else:
                left, right = _contact_set(pose), _contact_set(representative)
                fraction = len(left & right) / max(1, min(len(left), len(right)))
                similar = fraction >= fcc_cutoff
                method = f"fraction_common_contacts_ge_{fcc_cutoff}"
            if similar:
                members.append(pose); assigned[pose["pose_id"]] = f"cluster_{index:03d}"; pose["cluster_method"] = method; placed = True; break
        if not placed:
            clusters.append([pose]); assigned[pose["pose_id"]] = f"cluster_{len(clusters):03d}"
            pose["cluster_method"] = (f"ligand_heavy_atom_rmsd_le_{rmsd_cutoff}_A" if pose["job_type"] == "small_molecule" else f"fraction_common_contacts_ge_{fcc_cutoff}")
    return [{"pose_id": pose["pose_id"], "cluster_id": assigned[pose["pose_id"]], "engine": pose["engine"],
             "method": pose["cluster_method"], "engine_score": pose["engine_score"], "engine_score_unit": pose["engine_score_unit"],
             "seed": pose["seed"], "receptor_state": pose["receptor_state"]} for pose in sorted(poses, key=lambda row: row["pose_id"])]


def analyze(config: Mapping[str, Any], *, base_dir: str | Path) -> tuple[dict[str, Any], dict[str, str]]:
    base = Path(base_dir).resolve()
    receptor_spec = config.get("receptor")
    if not isinstance(receptor_spec, Mapping):
        raise InputError("receptor object is required")
    receptor_path = _resolve(base, receptor_spec.get("path"), "receptor")
    receptor_atoms = parse_atoms(receptor_path)
    receptor_sha = sha256(receptor_path)
    expected = str(receptor_spec.get("sha256", receptor_sha))
    if expected != receptor_sha:
        raise InputError("receptor checksum mismatch")
    cutoff = float(config.get("contact_cutoff_A", 5.0))
    if cutoff <= 0:
        raise InputError("contact cutoff must be positive")
    raw_poses = config.get("poses")
    if not isinstance(raw_poses, list) or not raw_poses:
        raise InputError("at least one pose is required")
    poses: list[dict[str, Any]] = []
    input_digests = {"receptor": receptor_sha}
    pose_ids = set()
    for raw in raw_poses:
        if not isinstance(raw, Mapping):
            raise InputError("pose records must be objects")
        pose_id = str(raw.get("pose_id", ""))
        if not pose_id or pose_id in pose_ids:
            raise InputError("pose_id values must be non-empty and unique")
        pose_ids.add(pose_id)
        engine = str(raw.get("engine", ""))
        if engine not in {"autodock_vina", "haddock3", "cluspro"}:
            raise InputError(f"unsupported pose engine {engine!r}")
        job_type = str(raw.get("job_type", ""))
        if job_type not in {"small_molecule", "protein_protein", "protein_peptide"}:
            raise InputError("unsupported pose job_type")
        path = _resolve(base, raw.get("path"), f"pose {pose_id}")
        all_atoms = parse_atoms(path)
        ligand_chains = {str(value) for value in raw.get("ligand_chains", [])}
        ligand_atoms = [atom for atom in all_atoms if atom["chain"] in ligand_chains] if ligand_chains else all_atoms
        if not ligand_atoms:
            raise InputError(f"pose {pose_id} selected no ligand atoms")
        contacts = _pose_contacts(receptor_atoms, ligand_atoms, pose_id, cutoff)
        score = raw.get("engine_score")
        try:
            score = None if score in (None, "") else float(score)
        except (TypeError, ValueError) as exc:
            raise InputError("engine_score must be numeric") from exc
        restraint_pairs = {tuple(item) for item in raw.get("required_contact_pairs", [])}
        observed_pairs = {(row["receptor_chain"], row["receptor_auth_seq_id"]) for row in contacts}
        native_site = {tuple(item) for item in config.get("native_site_residues", [])}
        poses.append({"pose_id": pose_id, "engine": engine, "job_type": job_type, "source_sha256": sha256(path),
                      "engine_score": score, "engine_score_unit": str(raw.get("engine_score_unit", "unreported")),
                      "seed": int(raw.get("seed", 0)), "receptor_state": str(raw.get("receptor_state", "unresolved")),
                      "ligand_atoms": ligand_atoms, "contacts": contacts, "contact_count": len(contacts),
                      "clash_count": sum(row["contact_state"] == "clash" for row in contacts),
                      "hydrogen_bond_proxy_count": sum(row["hydrogen_bond_proxy"] for row in contacts),
                      "salt_bridge_proxy_count": sum(row["salt_bridge_proxy"] for row in contacts),
                      "metal_coordination_proxy_count": sum(row["metal_coordination_proxy"] for row in contacts),
                      "restraint_satisfaction": None if not restraint_pairs else round(len(restraint_pairs & observed_pairs) / len(restraint_pairs), 6),
                      "native_site_overlap": None if not native_site else round(len(native_site & observed_pairs) / len(native_site), 6)})
        input_digests[f"pose_{pose_id}"] = sha256(path)
    clusters = _cluster(poses, float(config.get("ligand_rmsd_cutoff_A", 2.0)), float(config.get("fraction_common_contacts_cutoff", 0.5)))
    contacts = [contact for pose in poses for contact in pose["contacts"]]
    pose_records = [{key: value for key, value in pose.items() if key not in {"ligand_atoms", "contacts", "cluster_method"}} for pose in sorted(poses, key=lambda row: row["pose_id"])]
    seed_groups: dict[tuple[str, str], set[str]] = {}
    cluster_by_pose = {row["pose_id"]: row["cluster_id"] for row in clusters}
    for pose in pose_records:
        seed_groups.setdefault((pose["engine"], pose["receptor_state"]), set()).add(cluster_by_pose[pose["pose_id"]])
    consistency = [{"engine": key[0], "receptor_state": key[1], "distinct_clusters": len(value)} for key, value in sorted(seed_groups.items())]
    return ({
        "schema_version": SCHEMA_VERSION, "module_id": "pose_evidence", "subject": dict(config.get("subject", {"id": "unknown"})),
        "receptor_sha256": receptor_sha, "contact_cutoff_A": cutoff, "poses": pose_records, "clusters": clusters,
        "contacts": contacts, "repeated_seed_consistency": consistency,
        "state_sensitivity": {state: [row["pose_id"] for row in pose_records if row["receptor_state"] == state]
                              for state in sorted({row["receptor_state"] for row in pose_records})},
        "buried_sasa": {"state": "unevaluated", "reason": "FreeSASA or an explicitly validated SASA adapter was not invoked"},
        "limitations": ["Docking scores remain engine-specific and are not affinities.",
                        "Hydrogen bonds, salt bridges, clashes, and metal coordination are geometry proxies.",
                        "Imported ClusPro results are read-only and are never uploaded by this module."],
    }, input_digests)


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_tsv(path: Path, rows: Iterable[Mapping[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows: writer.writerow(row)


def write_outputs(out_dir: str | Path, document: Mapping[str, Any], *, input_sha256: Mapping[str, str], config_sha256: str) -> list[Path]:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    summary, cluster_table, contact_table, layer, run = (out / name for name in (
        "POSE_EVIDENCE.json", "POSE_CLUSTERS.tsv", "POSE_CONTACTS.tsv", "POSE_LAYER.json", "RUN_MANIFEST.json"))
    _write_json(summary, document)
    _write_tsv(cluster_table, document["clusters"], ["pose_id", "cluster_id", "engine", "method", "engine_score", "engine_score_unit", "seed", "receptor_state"])
    _write_tsv(contact_table, document["contacts"], ["pose_id", "receptor_chain", "receptor_auth_seq_id", "receptor_insertion_code", "receptor_resname",
                                                             "ligand_chain", "ligand_auth_seq_id", "ligand_insertion_code", "ligand_resname", "minimum_distance_A",
                                                             "contact_state", "hydrogen_bond_proxy", "salt_bridge_proxy", "metal_coordination_proxy"])
    residues = {}
    for row in document["contacts"]:
        key = (row["receptor_chain"], row["receptor_auth_seq_id"], row["receptor_insertion_code"])
        residues.setdefault(key, []).append(row["pose_id"])
    _write_json(layer, {"schema_version": SCHEMA_VERSION, "contract_id": "structure_layer_bundle", "subject": document["subject"],
                        "coordinate_sha256": document["receptor_sha256"], "layer_id": "pose_contacts",
                        "records": [{"chain_id": key[0], "auth_seq_id": key[1], "insertion_code": key[2], "metric": "docking_contact_count",
                                     "state": "modeled_pose_contact", "value": len(set(values)), "detail": ";".join(sorted(set(values))),
                                     "source_digest": document["receptor_sha256"], "evidence_class": "docking_geometry"}
                                    for key, values in sorted(residues.items())]})
    _write_json(run, {"schema_version": SCHEMA_VERSION, "module_id": "pose_evidence", "version": "0.1.0",
                      "input_sha256": dict(sorted({**input_sha256, "config": config_sha256}.items())),
                      "parameters": {"contact_cutoff_A": document["contact_cutoff_A"]},
                      "runtime_versions": {"python": platform.python_version(), "numpy": np.__version__},
                      "optional_runtimes": {"freesasa": "not_invoked"},
                      "outputs": [summary.name, cluster_table.name, contact_table.name, layer.name],
                      "missing_evidence": ["buried_sasa"] if document["buried_sasa"]["state"] == "unevaluated" else []})
    return [summary, cluster_table, contact_table, layer, run]
