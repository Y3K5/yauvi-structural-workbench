from __future__ import annotations
import csv
import hashlib
import json
import re
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from structqc.coordinate_scope import load_scope
from structqc.chirality import check_scope
from memorient.sasa import atom_sasa, _atom_radius
from memorient.geometry import ordered_intrinsic_rotation
from assembly_context.core import CONTACT_CUTOFF_A

VERSION = "0.1.0.dev0"
LIMITS = ["Experimental structural interpretation; not scientifically qualified.",
          "Geometry and solvent exposure do not establish biological reachability, affinity, function, or immunogenicity.",
          "Static coordinates do not establish dynamic probabilities or native assembly completeness.",
          "Unrepresented neighbors and uncomputed chemistry, curvature and probe reachability remain unknown."]


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def binding(scope):
    return {k: scope[k] for k in ("coordinate_sha256", "model_id", "assembly_id", "frame", "conformer_policy")}


def _identifier(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:-]{0,95}", value):
        raise ValueError("Invalid side or object identifier")
    return value


def validate_declaration(declaration, scope):
    if declaration.get("schema_version") != "1.0" or declaration.get("binding") != binding(scope):
        raise ValueError("Patch declaration must bind to the exact coordinate/model/assembly/frame/conformer scope")
    sides = declaration.get("sides", [])
    if not isinstance(sides, list) or not sides:
        raise ValueError("Declare at least one side")
    seen = set()
    copies = {c["chain_id"] for c in scope["copies"]}
    for side in sides:
        sid = _identifier(side.get("side_id"))
        if sid in seen:
            raise ValueError("Duplicate side identity")
        seen.add(sid)
        if not side.get("source", {}).get("id") or not side.get("source", {}).get("citation"):
            raise ValueError("Every patch needs its declaration or annotation source")
        if side.get("basis", "curated") not in {"curated", "experimental", "predicted", "derived"}:
            raise ValueError("Invalid patch evidence basis")
        roles = side.get("biological_roles", [])
        if not isinstance(roles, list) or any(not isinstance(role, str) or not role.strip() for role in roles):
            raise ValueError("Biological roles must be explicitly declared text labels")
        if side.get("equivalence_group") is not None and (not isinstance(side["equivalence_group"], str) or not side["equivalence_group"].strip()):
            raise ValueError("An equivalence group must be a declared text label")
        keys = set()
        for row in side.get("residue_set", []):
            if row.get("chain_id") not in copies or type(row.get("auth_seq_id")) is not int:
                raise ValueError("Patch requires an exact chain copy and author residue identity")
            key = (row["chain_id"], row["auth_seq_id"], row.get("insertion_code", ""))
            if key in keys:
                raise ValueError("Repeated patch residue identity")
            keys.add(key)
        if not keys:
            raise ValueError("Patch residue set is empty")
    return sides


def _evidence(path, expected):
    doc = json.loads(Path(path).read_text())
    if doc.get("binding") != expected:
        raise ValueError("Supporting evidence has a different coordinate/model/assembly/frame/conformer scope")
    if doc.get("engine") not in {"memorient", "assembly_context", "site_context", "state_atlas", "actstate", "patch_measurements"}:
        raise ValueError("Unsupported Bio-Orient evidence adapter")
    method = doc.get("method", {})
    if not method.get("id") or not re.fullmatch(r"[0-9a-f]{64}", str(method.get("source_sha256", ""))):
        raise ValueError("Supporting evidence requires method identity")
    if doc.get("output_sha256") != hashlib.sha256(canonical(doc.get("output"))).hexdigest():
        raise ValueError("Supporting evidence output checksum mismatch")
    output = doc["output"]
    if not isinstance(output, dict):
        raise ValueError("Supporting evidence output must be an object")
    internal = output.get("input_sha256", {}).get("structure")
    if internal and internal != expected["coordinate_sha256"]:
        raise ValueError("Engine coordinate identity conflicts with evidence binding")
    return doc


def _sampling_frame(coords):
    # Sampling in an intrinsic proper frame avoids a world-fixed sphere making
    # patch directions and burial depend on how the input file was rotated.
    try:
        rotation, origin, _ = ordered_intrinsic_rotation(coords)
    except ValueError:
        rotation, origin = np.eye(3), coords.mean(axis=0)
    return (coords - origin) @ rotation.T, rotation


def _surface(coords, radii):
    sampled_coords, rotation = _sampling_frame(coords)
    area, vectors = atom_sasa(sampled_coords, radii, return_directions=True)
    return area, vectors @ rotation


def _membrane(scope, declaration, topology, supporting):
    from memorient.contexts import get_context
    from memorient.geometry import load_structure
    from memorient.orientor import orient_structure
    found = [d for d in supporting if d["engine"] == "memorient"]
    if len(found) > 1:
        raise ValueError("Conflicting placements must be separate alternatives")
    membrane = None
    if found:
        membrane = found[0]["output"].get("input_coordinate_membrane")
        if membrane is None:
            membrane = found[0]["output"].get("membrane")
    elif declaration.get("membrane_context"):
        chain = declaration.get("membrane_chain")
        if not chain:
            raise ValueError("Membrane placement requires an explicitly selected source chain")
        context = get_context(declaration["membrane_context"])
        if not context.has_bilayer:
            return None
        structure = load_structure(declaration["_structure_path"], chain=chain, model=scope["model_index"])
        result = orient_structure(structure, context, topology_evidence=topology, validate=False)
        membrane = result.to_dict().get("input_coordinate_membrane")
    if not membrane:
        return None
    normal = np.asarray(membrane["normal"], dtype=float)
    center = np.asarray(membrane["center"], dtype=float)
    half = float(membrane["half_thickness"])
    if (normal.shape != (3,) or center.shape != (3,) or not np.isfinite([*normal, *center, half]).all()
            or half <= 0 or not np.isclose(np.linalg.norm(normal), 1, atol=1e-6)):
        raise ValueError("Invalid input-frame membrane")
    if membrane.get("frame") not in {"input_coordinates", "deposited_coordinates"}:
        raise ValueError("Membrane is in a different coordinate frame")
    from memorient.sidedness import verify_placement_sidedness
    authors = {c["chain_id"]: c["auth_asym_id"] for c in scope["copies"]}
    marker_atoms = [{"atom": a["atom"], "chain_id": authors[a["chain"]], "chain_copy_id": a["chain"],
                     "auth_seq_id": a["resi"], "insertion_code": a["icode"], "xyz": [a[k] for k in "xyz"]} for a in scope["atoms"]]
    assignment = verify_placement_sidedness(marker_atoms, membrane)
    return {**membrane, "side_assignment": assignment, "sidedness": assignment["state"], "frame": "deposited_coordinates", "binding": binding(scope)}


def analyze(structure_path, declaration_path, *, model_index=0, assembly_id="asu",
            chemical_reference=None, topology_path=None, evidence_paths=()):
    scope = load_scope(structure_path, model_index, assembly_id)
    declaration = json.loads(Path(declaration_path).read_text())
    patches = validate_declaration(declaration, scope)
    expected = binding(scope)
    supporting = [_evidence(p, expected) for p in evidence_paths]
    patch_by_id = {p["side_id"]: p for p in patches}
    for record in supporting:
        for item in record["output"].get("side_records", []):
            patch = patch_by_id.get(item.get("side_id"))
            if patch is None or sorted(item.get("residue_set", []), key=canonical) != sorted(patch["residue_set"], key=canonical):
                raise ValueError("Side evidence must name the exact declared residue set and chain copies")
    topology = json.loads(Path(topology_path).read_text()) if topology_path else None
    if topology and topology.get("coordinate_sha256") != scope["coordinate_sha256"]:
        raise ValueError("Topology coordinate checksum mismatch")
    if topology and topology.get("binding") != expected:
        raise ValueError("Bio-Orient topology requires an exact model/assembly/frame/conformer binding")
    if topology and not declaration.get("membrane_context"):
        raise ValueError("Topology needs a declared membrane context")
    declaration["_structure_path"] = str(structure_path)
    membrane = _membrane(scope, declaration, topology, supporting)
    chemistry = check_scope(scope, chemical_reference)
    heavy = [a for a in scope["atoms"] if a["elem"] not in {"H", "D"}]
    if not heavy:
        raise ValueError("No heavy atoms in the selected scope")
    coords = np.array([[a[k] for k in "xyz"] for a in heavy])
    radii = np.array([_atom_radius(a["elem"], a["atom"]) for a in heavy])
    sampled_coords, sampling_rotation = _sampling_frame(coords)
    areas, vectors = atom_sasa(sampled_coords, radii, return_directions=True)
    vectors = vectors @ sampling_rotation
    indices = {a["index"]: i for i, a in enumerate(heavy)}
    components = {c["component_id"]: c for c in scope["components"]}
    groups = {}
    for i, atom in enumerate(heavy):
        groups.setdefault(atom["chain"], []).append(i)
    isolated = {}
    for chain, members in groups.items():
        retained = set(members)
        excluded = [i for i in range(len(heavy)) if i not in retained]
        sampled = atom_sasa(sampled_coords, radii, target_indices=members, excluded_indices=excluded)
        isolated[chain] = {i: sampled[i] for i in members}
    objects = []
    for chain, members in sorted(groups.items()):
        polymer = [i for i in members if components[heavy[i]["properties"]["component_id"]]["is_polymer"]]
        if polymer:
            kinds = {components[heavy[i]["properties"]["component_id"]]["component_kind"] for i in polymer}
            objects.append({"object_id": "chain:" + chain, "kind": next(iter(kinds)) if len(kinds) == 1 else "polymer", "chain_id": chain, "indices": polymer})
    for cid, component in sorted(components.items()):
        if not component["is_polymer"]:
            members = [indices[i] for i in component["atom_indices"] if i in indices]
            if members:
                objects.append({"object_id": "component:" + cid, "kind": component["component_kind"],
                                "component_id": cid, "chemical_component_id": component["chemical_component_id"], "indices": members})
    for obj in objects:
        obj["centroid"] = coords[obj["indices"]].mean(axis=0).tolist()
        obj["atom_indices"] = [heavy[i]["index"] for i in obj["indices"]]
    for annotation in declaration.get("object_annotations", []):
        matches = [o for o in objects if o["object_id"] == annotation.get("object_id")]
        if len(matches) != 1 or not annotation.get("source", {}).get("id") or not annotation.get("source", {}).get("citation"):
            raise ValueError("Object annotation requires an exact object and source")
        if annotation.get("kind") not in {"protein", "nucleic_acid", "ligand", "glycan", "lipid", "solvent", "chemical_component"}:
            raise ValueError("Unsupported object category")
        matches[0].update(kind=annotation["kind"], annotation=annotation)
    context_objects = declaration.get("context_objects", [])
    for obj in context_objects:
        _identifier(obj.get("object_id"))
        if obj.get("kind") not in {"compartment", "cavity", "membrane", "solvent"} or not obj.get("source", {}).get("id") or not obj.get("source", {}).get("citation"):
            raise ValueError("Context objects require an explicit kind and cited source")
        if "centroid" in obj and (np.asarray(obj["centroid"]).shape != (3,) or not np.isfinite(obj["centroid"]).all()):
            raise ValueError("Invalid context centroid")
    if membrane:
        context_objects = [*context_objects, {"object_id": "modeled_membrane", "kind": "membrane", "centroid": membrane["center"], "normal": membrane["normal"]}]
    if len({o["object_id"] for o in [*objects, *context_objects]}) != len(objects) + len(context_objects):
        raise ValueError("Duplicate neighbor/context object identity")
    sides, edges = [], []
    for patch in patches:
        wanted = {(r["chain_id"], r["auth_seq_id"], r.get("insertion_code", "")) for r in patch["residue_set"]}
        matches = [c for c in scope["components"] if (c["chain_id"], c["auth_seq_id"], c["insertion_code"]) in wanted]
        # A single author residue identity cannot ambiguously name two components.
        matched_keys = [(c["chain_id"], c["auth_seq_id"], c["insertion_code"]) for c in matches]
        if len(set(matched_keys)) != len(matched_keys):
            raise ValueError("Patch residue identity maps to multiple components")
        members = sorted({indices[i] for c in matches for i in c["atom_indices"] if i in indices})
        missing = sorted(wanted - set(matched_keys))
        area = float(areas[members].sum())
        directional = vectors[members].sum(axis=0)
        coherence = float(np.linalg.norm(directional) / area) if area else 0.0
        normal = (directional / np.linalg.norm(directional)).tolist() if coherence >= 0.25 else None
        side = {"side_id": patch["side_id"], "residue_set": patch["residue_set"],
                "atom_indices": [heavy[i]["index"] for i in members],
                "centroid": coords[members].mean(axis=0).tolist() if members else None,
                "surface_normal": normal, "normal_state": "sampled" if normal else "unresolved",
                "directional_coherence": coherence, "coverage": {"declared_residues": len(wanted), "observed_residues": len(matches), "missing_residue_keys": [list(k) for k in missing]},
                "surface_chemistry": {"element_counts": {el: sum(heavy[i]["elem"] == el for i in members) for el in sorted({heavy[i]["elem"] for i in members})}, "component_ids": sorted({c["chemical_component_id"] for c in matches})},
                "solvent_exposure": {"sasa_A2": area, "isolated_chain_sasa_A2": sum(isolated[heavy[i]["chain"]][i] for i in members)},
                "curvature": {"state": "unevaluated"}, "mobility": {"state": "unknown"},
                "compartment_orientation": {"state": "unknown"}, "functional_annotation": patch.get("functional_annotation", []),
                "dynamic_state": {"state": "unknown"}, "functional_availability": {"state": "unknown"}, "equivalence_group": patch.get("equivalence_group"),
                "biological_roles": patch.get("biological_roles", []), "source": patch["source"], "basis": patch.get("basis", "curated")}
        if membrane and members:
            depths = (coords[members] - membrane["center"]) @ np.array(membrane["normal"])
            assignment = membrane.get("side_assignment", {"state": "unknown"})
            signs = set(1 if d > membrane["half_thickness"] else -1 if d < -membrane["half_thickness"] else 0 for d in depths)
            side["compartment_orientation"] = {"state": assignment["state"], "depth_range_A": [float(depths.min()), float(depths.max())],
                "compartment": assignment.get("positive_compartment") if signs == {1} else assignment.get("negative_compartment") if signs == {-1} else None,
                "geometry": "mixed_or_core" if len(signs) > 1 or signs == {0} else "positive_side" if signs == {1} else "negative_side"}
        for record in supporting:
            for item in record["output"].get("side_records", []):
                if item.get("side_id") not in {p["side_id"] for p in patches}:
                    raise ValueError("Evidence names an unknown side")
                if item.get("side_id") == side["side_id"]:
                    for field in ("curvature", "mobility", "dynamic_state", "functional_annotation", "functional_availability"):
                        if field in item:
                            if field == "functional_availability" and (not item.get("source", {}).get("id") or not item.get("source", {}).get("citation") or not item.get("interpretation_limit")):
                                raise ValueError("Functional availability needs scoped source evidence and an interpretation limit")
                            side[field] = {"state": "recorded_method_output", "value": item[field], "engine": record["engine"], "method": record["method"]}
        orientation = side["compartment_orientation"]
        side["availability"] = {"structural": "coordinates_present" if members and not missing else "partial_or_missing",
                                "orientational": orientation["state"] if orientation.get("compartment") else "conflicting" if orientation["state"] == "conflicting" else "unresolved" if membrane else "unknown",
                                "functional": side["functional_availability"]["state"]}
        sides.append(side)
        for obj in objects:
            neighbor = [i for i in obj["indices"] if i not in members]
            if not members or not neighbor:
                continue
            nearest, _ = cKDTree(coords[neighbor]).query(coords[members])
            distances = cKDTree(coords[members]).sparse_distance_matrix(cKDTree(coords[neighbor]), CONTACT_CUTOFF_A, output_type="coo_matrix")
            edge = {"source_side": side["side_id"], "adjacent_object": obj["object_id"], "minimum_distance_A": float(nearest.min()),
                    "contact_atom_pairs": int(distances.nnz), "contact_cutoff_A": CONTACT_CUTOFF_A,
                    "contact_area": {"state": "unevaluated"}, "buried_sasa_A2": None,
                    "steric_occlusion": {"state": "geometry_only"}, "complementarity": {"state": "unevaluated"},
                    "interaction_class": "geometric_contact" if distances.nnz else "no_contact_at_cutoff",
                    "relative_angle_degrees": None, "orientation_state": "unresolved", "dynamic_probability": None,
                    "functional_availability": {"state": "unknown"}}
            # Pairwise burial is recomputed with this neighbor removed; total
            # assembly burial cannot be attributed independently to each object.
            # Preserve exactly the full-assembly surface sample points. A distant
            # neighbor cannot hide any probe-expanded point and needs no rerun.
            if nearest.min() >= float((radii[members].max() + radii[neighbor].max()) + 2 * 1.4):
                edge["buried_sasa_A2"] = 0.0
            else:
                remaining = atom_sasa(sampled_coords, radii, target_indices=members, excluded_indices=neighbor)
                edge["buried_sasa_A2"] = max(0.0, float(remaining[members].sum()) - area)
            if normal:
                direction = np.array(obj["centroid"]) - side["centroid"]
                if np.linalg.norm(direction) > 1e-8:
                    edge["relative_angle_degrees"] = float(np.degrees(np.arccos(np.clip(np.dot(normal, direction / np.linalg.norm(direction)), -1, 1))))
                    edge["orientation_state"] = "measured_direction_to_neighbor"
            edges.append(edge)
        for obj in context_objects:
            edges.append({"source_side": side["side_id"], "adjacent_object": obj["object_id"],
                          "minimum_distance_A": None, "relative_angle_degrees": None,
                          "interaction_class": "declared_context", "dynamic_probability": None,
                          "functional_availability": {"state": "unknown"}})
    # Directed patch normals can also be compared without asserting an interaction.
    side_angles = []
    for i, a in enumerate(sides):
        for b in sides[i + 1:]:
            angle = None
            if a["surface_normal"] and b["surface_normal"]:
                angle = float(np.degrees(np.arccos(np.clip(np.dot(a["surface_normal"], b["surface_normal"]), -1, 1))))
            side_angles.append({"side_a": a["side_id"], "side_b": b["side_id"], "normal_angle_degrees": angle,
                                "state": "measured_directed_normals" if angle is not None else "unresolved",
                                "functional_relationship": "unknown"})
    for edge in edges:
        evidence = []
        for record in supporting:
            for item in record["output"].get("edge_records", []):
                if item.get("source_side") not in patch_by_id or item.get("adjacent_object") not in {o["object_id"] for o in [*objects, *context_objects]}:
                    raise ValueError("Edge evidence names an unknown side or neighbor")
                if sorted(item.get("residue_set", []), key=canonical) != sorted(patch_by_id[item["source_side"]]["residue_set"], key=canonical):
                    raise ValueError("Edge evidence residue scope differs from declared patch")
                if item["source_side"] == edge["source_side"] and item["adjacent_object"] == edge["adjacent_object"]:
                    source = item.get("source", {})
                    if not source.get("id") or not source.get("citation") or not item.get("interpretation_limit"):
                        raise ValueError("Biological edge annotations need a source and interpretation limit")
                    evidence.append({"record": item, "engine": record["engine"], "method": record["method"]})
        if evidence:
            edge["biological_evidence"] = evidence
            functional = [e for e in evidence if "functional_availability" in e["record"]]
            if functional:
                values = {canonical(e["record"]["functional_availability"]) for e in functional}
                edge["functional_availability"] = {"state": "conflicting" if len(values) > 1 else "recorded_evidence", "records": functional}
        edge["availability"] = {"structural": "measured_geometry" if edge.get("minimum_distance_A") is not None else "unknown",
                                "orientational": edge.get("orientation_state", "unresolved"), "functional": edge["functional_availability"]["state"]}
    role_sets = [set(s["biological_roles"]) for s in sides if s["biological_roles"]]
    equivalences = sorted({s["equivalence_group"] for s in sides if s["equivalence_group"]})
    descriptors = {"side_count": len(sides), "multisided": len(sides) > 1, "declared_equivalence_groups": equivalences,
                   "homosided": True if any(sum(s["equivalence_group"] == g for s in sides) > 1 for g in equivalences) else None,
                   "heterosided": True if len({tuple(sorted(r)) for r in role_sets}) > 1 else None,
                   "basis": "declared_and_evidence_mapped; symmetry is not functional equivalence"}
    return {"schema_version": "1.0", "module_id": "bio_orient", "binding": expected,
            "sides": sides, "objects": [{k: v for k, v in obj.items() if k != "indices"} for obj in objects] + context_objects,
            "edges": edges, "side_angles": side_angles, "sidedness_descriptors": descriptors, "membrane": membrane,
            "chirality": chemistry, "atoms": scope["atoms"], "components": scope["components"], "copies": scope["copies"],
            "recorded_state_evidence": [r for r in supporting if r["engine"] in {"state_atlas", "actstate"}],
            "supporting_evidence": supporting, "input_sha256": {"structure": sha(structure_path), "patch_declaration": sha(declaration_path),
                **({"chemical_reference": sha(chemical_reference)} if chemical_reference else {}), **({"topology_evidence": sha(topology_path)} if topology_path else {}),
                **{f"evidence_{i}": sha(p) for i, p in enumerate(evidence_paths)}},
            "software_identity": software_identity(),
            "methods": {"surface": "memorient_shrake_rupley_intrinsic_frame", "probe_radius_A": 1.4, "surface_samples": 240,
                        "minimum_directional_coherence": 0.25, "contacts": "assembly_context_cutoff", "contact_cutoff_A": CONTACT_CUTOFF_A},
            "limitations": LIMITS}


def software_identity():
    from importlib import metadata
    import platform
    import structqc.coordinate_scope, structqc.chirality, memorient.sasa, memorient.geometry, memorient.orientor, memorient.sidedness, memorient.contexts, assembly_context.core
    modules = [structqc.coordinate_scope, structqc.chirality, memorient.sasa, memorient.geometry, memorient.orientor, memorient.sidedness, memorient.contexts, assembly_context.core]
    return {"source_sha256": {"bio_orient.core": sha(__file__), "bio_orient.cli": sha(Path(__file__).with_name("cli.py")), **{m.__name__: sha(m.__file__) for m in modules}},
            "runtime_versions": {"python": platform.python_version(), **{p: metadata.version(p) for p in ("numpy", "scipy", "gemmi", "biopython")}}}


def write_outputs(document, output_dir):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    def write(name, data):
        (out / name).write_bytes(json.dumps(data, indent=2, sort_keys=True, allow_nan=False).encode() + b"\n")
    write("BIO_ORIENT.json", document)
    write("BIO_ORIENT_LAYER.json", {**document, "layer_id": "bio_orient"})
    for filename, rows, fields in [("BIO_ORIENT_SIDES.tsv", document["sides"], ["side_id", "residue_set", "coverage", "centroid", "surface_normal", "availability"]),
                                   ("BIO_ORIENT_EDGES.tsv", document["edges"], ["source_side", "adjacent_object", "minimum_distance_A", "relative_angle_degrees", "buried_sasa_A2", "functional_availability"])]:
        with (out / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows({k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v for k, v in row.items()} for row in rows)
    write("RUN_MANIFEST.json", {"schema_version": "1.0", "module_id": "bio_orient", "version": VERSION,
         "binding": document["binding"], "input_sha256": document["input_sha256"], "methods": document["methods"],
         "software_identity": document["software_identity"],
         "reference_identity": document["chirality"]["reference_identity"], "scientific_readiness": "experimental_nonblocking",
         "limitations": LIMITS, "outputs": {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "RUN_MANIFEST.json"}})
