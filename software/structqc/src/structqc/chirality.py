"""Reference-relative tetrahedral handedness; no repair or invented R/S labels."""
from __future__ import annotations
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from .coordinate_scope import load_scope


def _rows(block, prefix, fields):
    columns = [list(block.find_values(prefix + field)) for field in fields]
    if not columns[0]:
        return []
    if any(c and len(c) != len(columns[0]) for c in columns):
        raise ValueError("Incomplete CCD category")
    columns = [c or ["?"] * len(columns[0]) for c in columns]
    return [dict(zip(fields, values)) for values in zip(*columns)]


def read_references(manifest_path):
    """A local JSON manifest locks every CCD file; paths stay local and relative."""
    import gemmi
    path = Path(manifest_path)
    manifest = json.loads(path.read_text())
    if manifest.get("schema_version") != "1.0":
        raise ValueError("Unsupported chemical reference manifest")
    references, hashes = {}, {}
    for item in manifest.get("components", []):
        relative = Path(item["path"])
        file = (path.parent / relative).resolve()
        if relative.is_absolute() or ".." in relative.parts or not file.is_relative_to(path.parent.resolve()):
            raise ValueError("Chemical reference escapes manifest directory")
        raw = file.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != item["sha256"]:
            raise ValueError("Chemical reference checksum mismatch")
        cid = item["component_id"]
        if cid in references:
            raise ValueError("Duplicate chemical reference")
        document = gemmi.cif.read_string(raw.decode())
        blocks = [b for b in document if b.find_value("_chem_comp.id") == cid]
        if len(blocks) != 1:
            raise ValueError("CCD component identity mismatch")
        block = blocks[0]
        atoms = _rows(block, "_chem_comp_atom.", ["atom_id", "type_symbol", "pdbx_stereo_config",
                     "pdbx_model_Cartn_x_ideal", "pdbx_model_Cartn_y_ideal", "pdbx_model_Cartn_z_ideal"])
        bonds = _rows(block, "_chem_comp_bond.", ["atom_id_1", "atom_id_2"])
        atom_map = {a["atom_id"]: a for a in atoms}
        if len(atom_map) != len(atoms):
            raise ValueError("Duplicate CCD atom identity")
        references[cid] = {"atoms": atom_map, "bonds": bonds, "sha256": digest, "component_type": block.find_value("_chem_comp.type") or "unknown",
                           "ambiguous": block.find_value("_chem_comp.pdbx_ambiguous_flag") == "Y"}
        hashes[cid] = digest
    return references, {"manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "component_sha256": hashes}


def _volume(center, neighbors):
    vectors = np.array(neighbors) - np.array(center)
    return float(np.linalg.det(vectors))


def check_scope(scope, reference_manifest=None):
    references, identity = read_references(reference_manifest) if reference_manifest else ({}, {})
    findings, coverage = [], []
    for component in scope["components"]:
        cid = component["chemical_component_id"]
        reference = references.get(cid)
        row = {**{k: component[k] for k in ("component_id", "chemical_component_id", "chain_id", "auth_seq_id", "insertion_code", "selected_altloc")},
               "state": "reference_missing", "centers_declared": None, "centers_evaluated": 0}
        if reference:
            centers = [a for a in reference["atoms"].values() if a["pdbx_stereo_config"] in {"R", "S"}]
            row.update(state="no_declared_tetrahedral_centers" if not centers else "evaluated", centers_declared=len(centers))
            observed = {scope["atoms"][i]["atom"]: scope["atoms"][i] for i in component["atom_indices"]}
            for center in centers:
                finding = {**{k: row[k] for k in ("component_id", "chemical_component_id", "chain_id", "auth_seq_id", "insertion_code", "selected_altloc")},
                           "atom_id": center["atom_id"], "expected_reference_configuration": center["pdbx_stereo_config"],
                           "reference_sha256": reference["sha256"], "state": "unevaluated",
                           "reason": "unsupported_or_missing_reference_geometry", "observed_configuration": None}
                neighbors = sorted({b["atom_id_2"] if b["atom_id_1"] == center["atom_id"] else b["atom_id_1"]
                                    for b in reference["bonds"] if center["atom_id"] in (b["atom_id_1"], b["atom_id_2"])})
                heavy = [n for n in neighbors if n in reference["atoms"] and reference["atoms"][n]["type_symbol"] not in {"H", "D"}]
                if reference["ambiguous"]:
                    finding["reason"] = "ambiguous_reference_chemistry"
                elif component["alternate_locations"]:
                    finding["reason"] = "alternative_conformers"
                elif any(n not in observed for n in [center["atom_id"], *heavy]):
                    finding["reason"] = "missing_required_atoms"
                elif any(observed[n]["elem"] != reference["atoms"][n]["type_symbol"] for n in [center["atom_id"], *heavy]):
                    finding["reason"] = "element_identity_mismatch"
                elif any(not np.isfinite(observed[n]["occupancy"]) or observed[n]["occupancy"] <= 0 for n in [center["atom_id"], *heavy]):
                    finding["reason"] = "unrepresented_atom_occupancy"
                elif len(neighbors) == 4 and len(heavy) >= 3:
                    try:
                        ideal = lambda n: [float(reference["atoms"][n]["pdbx_model_Cartn_" + ax + "_ideal"]) for ax in "xyz"]
                        actual = lambda n: [observed[n][ax] for ax in "xyz"]
                        triples = [(abs(_volume(ideal(center["atom_id"]), [ideal(n) for n in triple])), triple)
                                   for triple in itertools.combinations(heavy, 3)]
                        _, triple = max(triples, key=lambda v: (v[0], v[1]))
                        target = _volume(ideal(center["atom_id"]), [ideal(n) for n in triple])
                        value = _volume(actual(center["atom_id"]), [actual(n) for n in triple])
                        if np.isfinite([target, value]).all() and abs(target) > 1e-6:
                            finding.update(reference_volume_A3=target, observed_volume_A3=value, neighbor_atom_ids=list(triple))
                            # A nearly planar center cannot substantiate handedness.
                            if abs(value) < 0.1 * abs(target):
                                finding["reason"] = "near_planar_geometry"
                            else:
                                finding.update(state="matches_reference" if value * target > 0 else "reference_mismatch", reason="signed_volume_comparison")
                                row["centers_evaluated"] += 1
                    except (ValueError, KeyError):
                        pass
                findings.append(finding)
            if row["centers_evaluated"] < len(centers):
                row["state"] = "partial_or_unevaluated"
            if (not centers and (component["alternate_locations"] or not reference["atoms"]
                    or any(a["pdbx_stereo_config"] not in {"N", "R", "S"} for a in reference["atoms"].values())
                    or any(n not in observed for n, a in reference["atoms"].items() if a["type_symbol"] not in {"H", "D"}))):
                row["state"] = "partial_or_unevaluated"
            row["reference_component_type"] = reference["component_type"]
        coverage.append(row)
    return {"state": "review_required" if any(f["state"] == "reference_mismatch" for f in findings) else "partial" if any(c["state"] in {"reference_missing", "partial_or_unevaluated"} for c in coverage) else "reference_checks_completed",
            "binding": {k: scope[k] for k in ("coordinate_sha256", "model_id", "model_index", "assembly_id", "frame", "conformer_policy")},
            "method": "ccd_atom_mapped_signed_volume_v1", "reference_identity": identity,
            "near_planar_relative_volume_threshold": 0.1, "findings": findings, "component_coverage": coverage,
            "limitations": ["Reference-relative tetrahedral handedness only; no absolute configuration is inferred.",
                            "Missing references or atoms and alternative conformers are not clean results.",
                            "No coordinates are repaired; mismatches are review flags, not automatic downstream holds.",
                            "Non-tetrahedral, linkage-dependent and unspecified stereochemistry are not qualified."]}


def analyze_chirality(path, *, model_index=0, chain=None, assembly_id="asu", reference_manifest=None):
    scope = load_scope(path, model_index, assembly_id, chain)
    result = check_scope(scope, reference_manifest)
    result["component_layer"] = {"binding": {k: scope[k] for k in ("coordinate_sha256", "model_id", "assembly_id", "frame", "conformer_policy")},
                                 "layer_kind": "component_chemistry", "atoms": scope["atoms"], "components": scope["components"]}
    return result
