from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from pathlib import Path
import shutil
from typing import Any, Iterable, Mapping

import Bio
import numpy as np
import scipy
from Bio.PDB import MMCIFParser, PDBParser

SCHEMA_VERSION = "1.0"
SITE_TYPES = {"active_site", "binding_site", "metal_ligand", "disulfide", "other_site"}
ROLES = {"nucleophile", "acid_base", "charge_relay", "metal_ligand", "unspecified"}
WATER = {"HOH", "WAT", "DOD"}


class InputError(RuntimeError):
    pass


def sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def read_json(path: str | Path | None) -> Any:
    if path is None:
        return None
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(f"cannot read JSON {path}: {exc}") from exc


def read_annotations(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    if source.suffix.lower() == ".json":
        document = read_json(source)
        if not isinstance(document, dict):
            raise InputError("annotation JSON must be an object")
        return document
    delimiter = "," if source.suffix.lower() == ".csv" else "\t"
    with source.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    sites = []
    for row in rows:
        try:
            position = int(row.get("position", ""))
        except ValueError as exc:
            raise InputError("every annotation row requires an integer position") from exc
        sites.append({
            "position": position, "type": row.get("type", "other_site"),
            "role": row.get("role", "unspecified"),
            "expected_residues": [x.strip().upper() for x in row.get("expected_residues", "").replace(",", ";").split(";") if x.strip()],
            "ligand_id": row.get("ligand_id", ""), "detail": row.get("detail", ""),
        })
    return {"sites": sites, "declared_cofactors": []}


def _load_structure(path: Path):
    try:
        parser = MMCIFParser(QUIET=True) if path.name.lower().endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
        return next(parser.get_structure("query", str(path)).get_models())
    except Exception as exc:
        raise InputError(f"cannot parse structure: {type(exc).__name__}: {exc}") from exc


def _mapping(manifest: Mapping[str, Any]) -> dict[int, list[dict[str, Any]]]:
    found: dict[int, list[dict[str, Any]]] = {}
    for residue in manifest.get("residues", []):
        position = residue.get("sequence_index")
        if position is not None:
            found.setdefault(int(position), []).append(dict(residue))
    return found


def _structure_residues(model) -> dict[tuple[str, int, str], Any]:
    return {
        (str(chain.id), int(residue.id[1]), str(residue.id[2]).strip()): residue
        for chain in model for residue in chain if residue.id[0] == " "
    }


def _heteroatoms(model) -> list[dict[str, Any]]:
    out = []
    for chain in model:
        for residue in chain:
            if residue.id[0] == " " or residue.resname.upper() in WATER:
                continue
            atoms = [a for a in residue.get_atoms() if str(a.element).upper() not in {"H", "D"}]
            if atoms:
                out.append({
                    "component_id": residue.resname.upper(), "chain_id": str(chain.id),
                    "auth_seq_id": int(residue.id[1]), "insertion_code": str(residue.id[2]).strip(),
                    "atoms": atoms,
                })
    return out


def _distance(residue, ligand: Mapping[str, Any]) -> float | None:
    left = [a for a in residue.get_atoms() if str(a.element).upper() not in {"H", "D"}]
    right = ligand["atoms"]
    if not left or not right:
        return None
    return min(float(np.linalg.norm(a.coord - b.coord)) for a in left for b in right)


def _component_map(raw: Mapping[str, Any] | None) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for key, value in (raw or {}).items():
        values = value if isinstance(value, list) else [value]
        out[str(key).upper()] = {str(v).upper() for v in values}
    return out


def _cofactor_calls(declared: Iterable[Mapping[str, Any]], observed: list[dict[str, Any]], mapping: Mapping[str, set[str]]) -> list[dict[str, Any]]:
    present = {str(item["component_id"]).upper() for item in observed}
    calls = []
    for item in declared:
        component = str(item.get("component_id", "")).upper()
        chebi = str(item.get("chebi", "")).upper()
        expected: set[str] = {component} if component else set()
        if chebi:
            expected |= set(mapping.get(chebi, ()))
        if not expected:
            state, detail = "unresolved", "no exact component id or mapped ChEBI identity was supplied"
        elif present & expected:
            state, detail = "observed_match", "observed component(s): " + ",".join(sorted(present & expected))
        else:
            state, detail = "not_observed", "expected component(s) not present: " + ",".join(sorted(expected))
        calls.append({"name": str(item.get("name", "")), "component_id": component, "chebi": chebi,
                      "state": state, "detail": detail})
    return calls


def _pockets(raw_documents: Iterable[Any]) -> list[dict[str, Any]]:
    pockets: list[dict[str, Any]] = []
    for document in raw_documents:
        if document is None:
            continue
        rows = document.get("pockets", []) if isinstance(document, dict) else document
        if not isinstance(rows, list):
            raise InputError("pocket result must contain a pockets list")
        for row in rows:
            method = str(row.get("method", document.get("method", "") if isinstance(document, dict) else "")).strip()
            if not method:
                raise InputError("every pocket requires a method")
            residues = row.get("residues", [])
            pockets.append({
                "method": method, "pocket_id": str(row.get("pocket_id", row.get("id", ""))),
                "score": row.get("score"), "residues": sorted(residues, key=lambda r: (str(r.get("chain_id", "")), int(r.get("auth_seq_id", 0)))),
            })
    return sorted(pockets, key=lambda p: (p["method"], p["pocket_id"]))


def analyze(
    manifest: Mapping[str, Any], structure_path: str | Path, annotations: Mapping[str, Any], *,
    component_map: Mapping[str, Any] | None = None, pocket_results: Iterable[Any] = (),
) -> dict[str, Any]:
    pocket_documents = list(pocket_results)
    structure_path = Path(structure_path)
    if sha256(structure_path) != str(manifest.get("coordinate", {}).get("sha256", "")):
        raise InputError("StructQC manifest checksum does not match structure")
    model = _load_structure(structure_path)
    residue_lookup, sequence_map = _structure_residues(model), _mapping(manifest)
    observed = _heteroatoms(model)
    sites: list[dict[str, Any]] = []
    representative: list[tuple[int, np.ndarray]] = []
    for raw in annotations.get("sites", []):
        try:
            position = int(raw["position"])
        except (KeyError, ValueError, TypeError) as exc:
            raise InputError("site position must be an integer reference-sequence position") from exc
        site_type = str(raw.get("type", "other_site"))
        role = str(raw.get("role", "unspecified"))
        if site_type not in SITE_TYPES or role not in ROLES:
            raise InputError(f"unknown site type/role at position {position}: {site_type}/{role}")
        mapped = sequence_map.get(position, [])
        if len(mapped) != 1:
            sites.append({"position": position, "type": site_type, "role": role, "state": "unresolved_mapping",
                          "detail": f"expected one coordinate mapping, found {len(mapped)}", "expected_residues": raw.get("expected_residues", [])})
            continue
        residue_record = mapped[0]
        key = (str(residue_record["chain_id"]), int(residue_record["auth_seq_id"]), str(residue_record.get("insertion_code", "")))
        residue = residue_lookup.get(key)
        if residue is None:
            sites.append({"position": position, "type": site_type, "role": role, "state": "missing_coordinates",
                          "detail": "mapped residue is absent from parsed coordinates", "expected_residues": raw.get("expected_residues", [])})
            continue
        observed_aa = str(residue_record.get("one_letter", "X"))
        expected = sorted(str(x).upper() for x in raw.get("expected_residues", []))
        if expected:
            state = "role_compatible" if observed_aa in expected else "role_mismatch"
            detail = f"observed {observed_aa}; expected one of {','.join(expected)} for declared role {role}"
        else:
            state, detail = "role_unresolved", f"role {role} has no declared residue expectation"
        nearest = []
        for ligand in observed:
            distance = _distance(residue, ligand)
            if distance is not None:
                nearest.append({"component_id": ligand["component_id"], "chain_id": ligand["chain_id"],
                                "auth_seq_id": ligand["auth_seq_id"], "distance_A": round(distance, 6)})
        nearest.sort(key=lambda item: (item["distance_A"], item["component_id"]))
        ca = residue["CA"].coord if "CA" in residue else None
        if ca is not None:
            representative.append((position, np.asarray(ca, dtype=float)))
        sites.append({
            "position": position, "type": site_type, "role": role, "state": state, "detail": detail,
            "expected_residues": expected, "observed_residue": observed_aa,
            "chain_id": key[0], "auth_seq_id": key[1], "insertion_code": key[2],
            "nearest_heteroatoms": nearest[:5], "annotation_detail": str(raw.get("detail", "")),
        })
    distances = [float(np.linalg.norm(a[1] - b[1])) for i, a in enumerate(representative) for b in representative[i + 1:]]
    pockets = _pockets(pocket_documents)
    return {
        "schema_version": SCHEMA_VERSION, "module_id": "site_context", "subject": manifest["subject"],
        "coordinate_sha256": manifest["coordinate"]["sha256"],
        "input_sha256": {
            "structure_evidence_manifest": _json_sha256(manifest),
            "structure": sha256(structure_path),
            "annotations": _json_sha256(annotations),
            **({"component_map": _json_sha256(component_map)} if component_map is not None else {}),
            **{f"pocket_result_{index + 1}": _json_sha256(value)
               for index, value in enumerate(pocket_documents)},
        },
        "config": {
            "annotation_sites": len(annotations.get("sites", [])),
            "component_map_supplied": component_map is not None,
            "pocket_methods": sorted({item["method"] for item in pockets}),
        },
        "sites": sorted(sites, key=lambda x: (x["position"], x["type"])),
        "site_geometry": {"representative": "CA", "pair_count": len(distances),
                          "maximum_separation_A": round(max(distances), 6) if distances else None,
                          "interpretation": "descriptive geometry; no universal catalytic threshold applied"},
        "observed_heteroatoms": [
            {k: v for k, v in item.items() if k != "atoms"} for item in observed
        ],
        "cofactors": _cofactor_calls(annotations.get("declared_cofactors", []), observed, _component_map(component_map)),
        "pockets": pockets,
        "missing_evidence": ["pocket_results"] if not pockets else [],
        "limitations": [
            "Annotated roles, observed ligands, geometry, and predicted pockets are separate evidence legs.",
            "Pocket scores are method-specific and are not a druggability score.",
            "Coordinate proximity does not establish catalysis, occupancy, affinity, or function.",
        ],
    }


def _json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_outputs(out_dir: str | Path, document: Mapping[str, Any]) -> list[Path]:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    summary, sites_path, pockets_path, layer, run = (out / n for n in (
        "SITE_CONTEXT.json", "SITE_RESIDUES.tsv", "POCKETS.tsv", "SITE_LAYER.json", "RUN_MANIFEST.json"))
    _json(summary, document)
    site_fields = ["position", "type", "role", "state", "observed_residue", "chain_id", "auth_seq_id", "insertion_code", "detail"]
    with sites_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=site_fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader(); writer.writerows(document["sites"])
    with pockets_path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["method", "pocket_id", "score", "residue_count"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for pocket in document["pockets"]:
            writer.writerow({"method": pocket["method"], "pocket_id": pocket["pocket_id"], "score": pocket["score"],
                             "residue_count": len(pocket["residues"])})
    layer_records = [
        {"chain_id": s.get("chain_id"), "auth_seq_id": s.get("auth_seq_id"), "insertion_code": s.get("insertion_code", ""),
         "metric": s["type"], "state": s["state"], "value": None, "detail": s["detail"],
         "source_digest": document["input_sha256"]["annotations"], "evidence_class": "annotation_mapped"}
        for s in document["sites"] if s.get("chain_id") is not None
    ]
    for pocket in document["pockets"]:
        for residue in pocket["residues"]:
            layer_records.append({"chain_id": residue.get("chain_id"), "auth_seq_id": residue.get("auth_seq_id"),
                                  "insertion_code": residue.get("insertion_code", ""), "metric": "predicted_pocket",
                                  "state": "predicted", "value": pocket.get("score"),
                                  "detail": f"{pocket['method']} pocket {pocket['pocket_id']}",
                                  "source_digest": _json_sha256(pocket),
                                  "evidence_class": pocket["method"]})
    _json(layer, {"schema_version": SCHEMA_VERSION, "contract_id": "structure_layer_bundle",
                  "subject": document["subject"], "coordinate_sha256": document["coordinate_sha256"],
                  "layer_id": "site_context", "records": sorted(layer_records, key=lambda r: (str(r["chain_id"]), int(r["auth_seq_id"]), r["metric"]))})
    _json(run, {"schema_version": SCHEMA_VERSION, "module_id": "site_context", "version": "0.1.0",
                "input_sha256": document["input_sha256"],
                "parameters": document["config"],
                "runtime_versions": {"python": platform.python_version(), "biopython": Bio.__version__,
                                     "numpy": np.__version__, "scipy": scipy.__version__},
                "optional_runtimes": {
                    "fpocket": ("precomputed_result_supplied" if "fpocket" in document["config"]["pocket_methods"]
                                else "available_not_invoked" if shutil.which("fpocket") else "not_available"),
                    "p2rank": ("precomputed_result_supplied" if "p2rank" in document["config"]["pocket_methods"]
                               else "available_not_invoked" if shutil.which("prank") else "not_available"),
                },
                "outputs": [summary.name, sites_path.name, pockets_path.name, layer.name],
                "missing_evidence": document["missing_evidence"]})
    return [summary, sites_path, pockets_path, layer, run]
