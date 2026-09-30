"""Exact atom/copy identities shared by local chemistry and Bio-Orient.

This reader never repairs coordinates or applies crystal symmetry implicitly.
"""
from __future__ import annotations
import hashlib
from pathlib import Path
import numpy as np


def proper_rotation(matrix, translation=None):
    matrix = np.asarray(matrix, dtype=float)
    if (matrix.shape != (3, 3) or not np.isfinite(matrix).all()
            or not np.allclose(matrix.T @ matrix, np.eye(3), atol=1e-6, rtol=0)
            or not np.isclose(np.linalg.det(matrix), 1, atol=1e-6, rtol=0)):
        raise ValueError("Biology-preserving transform requires a finite proper rotation")
    if translation is not None:
        translation = np.asarray(translation, dtype=float)
        if translation.shape != (3,) or not np.isfinite(translation).all():
            raise ValueError("Invalid transform translation")
    return matrix


def load_scope(path, model_index=0, assembly_id="asu", chain=None):
    import gemmi
    path = Path(path)
    try:
        structure = gemmi.read_structure(str(path), format=gemmi.CoorFormat.Detect)
    except RuntimeError as exc:
        raise ValueError(f"Coordinate parsing failed: {exc}") from exc
    structure.setup_entities()
    if type(model_index) is not int or not 0 <= model_index < len(structure):
        raise ValueError("Select an existing model index")
    model = structure[model_index]
    assembly_id = str(assembly_id)
    groups = {}
    for source_chain in model:
        for residue in source_chain:
            label = residue.subchain or source_chain.name
            groups.setdefault(label, []).append((source_chain.name, residue))
    copies = []
    if assembly_id == "asu":
        copies = [(label, "identity", np.eye(3), np.zeros(3)) for label in sorted(groups)]
    else:
        assemblies = [a for a in structure.assemblies if a.name == assembly_id]
        if len(assemblies) != 1:
            raise ValueError("Select an explicitly deposited assembly")
        for generator in assemblies[0].generators:
            labels = list(generator.subchains)
            if not labels:
                labels = [label for label, rows in groups.items()
                          if rows[0][0] in generator.chains]
            for operator in generator.operators:
                for label in labels:
                    if label not in groups:
                        raise ValueError("Assembly names an absent chain")
                    copies.append((label, operator.name, np.array(operator.transform.mat.tolist()),
                                   np.array(operator.transform.vec.tolist())))
    atoms, components, copy_records = [], [], []
    seen = set()
    for label, operator, rotation, translation in copies:
        proper_rotation(rotation, translation)
        copy_id = f"{label}:{operator}"
        if copy_id in seen:
            raise ValueError("Duplicate chain-copy identity")
        seen.add(copy_id)
        author = groups[label][0][0]
        if chain and chain not in {author, copy_id}:
            continue
        copy_records.append({"chain_id": copy_id, "auth_asym_id": author,
                             "label_asym_id": label, "operator_id": operator,
                             "rotation": rotation.tolist(), "translation": translation.tolist()})
        component_ids = set()
        for _, residue in groups[label]:
            icode = residue.seqid.icode.strip()
            cid = f"{copy_id}/{residue.seqid.num}{icode}/{residue.name}"
            if cid in component_ids:
                raise ValueError("Ambiguous component identity")
            component_ids.add(cid)
            alt_scores = {}
            for atom in residue:
                alt = atom.altloc.strip("\x00 ")
                if alt:
                    alt_scores.setdefault(alt, []).append(atom.occ)
            selected = sorted(alt_scores, key=lambda a: (-float(np.mean(alt_scores[a])), a))[0] if alt_scores else ""
            component = {"component_id": cid, "chemical_component_id": residue.name,
                         "chain_id": copy_id, "auth_asym_id": author, "label_asym_id": label,
                         "operator_id": operator, "auth_seq_id": residue.seqid.num,
                         "insertion_code": icode, "selected_altloc": selected,
                         "alternate_locations": sorted(alt_scores), "atom_indices": [],
                         "is_polymer": residue.entity_type == gemmi.EntityType.Polymer}
            info = gemmi.find_tabulated_residue(residue.name)
            component["component_kind"] = ("protein" if info.is_amino_acid() else "nucleic_acid" if info.is_nucleic_acid()
                                           else "glycan" if str(info.kind) in {"ResidueKind.PYR", "ResidueKind.KET"}
                                           else "solvent" if info.is_water() else "chemical_component")
            names = set()
            for atom in residue:
                alt = atom.altloc.strip("\x00 ")
                if alt and alt != selected:
                    continue
                if atom.name in names:
                    raise ValueError("Duplicate atom identity within selected conformer")
                names.add(atom.name)
                xyz = rotation @ np.array([atom.pos.x, atom.pos.y, atom.pos.z]) + translation
                if not np.isfinite(xyz).all():
                    raise ValueError("Nonfinite atom coordinates")
                index = len(atoms)
                atoms.append({"index": index, "serial": index, "atom": atom.name,
                              "elem": atom.element.name, "x": float(xyz[0]), "y": float(xyz[1]), "z": float(xyz[2]),
                              "chain": copy_id, "resi": residue.seqid.num, "icode": icode,
                              "occupancy": float(atom.occ),
                              "resn": residue.name, "hetflag": residue.het_flag != "A",
                              "properties": {"component_id": cid}, "altloc": alt})
                component["atom_indices"].append(index)
            components.append(component)
    if not atoms:
        raise ValueError("Selected scope contains no coordinates")
    if len(atoms) > 300_000:
        raise ValueError("Selected scope exceeds atom capacity")
    return {"coordinate_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "model_index": model_index, "model_id": str(getattr(model, "name", getattr(model, "num", ""))), "assembly_id": assembly_id,
            "frame": "deposited_coordinates", "atoms": atoms, "components": components, "copies": copy_records,
            "conformer_policy": "highest_mean_occupancy_per_component; chemistry with alternates remains unevaluated"}
