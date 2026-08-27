from __future__ import annotations

import copy
import csv
import hashlib
import importlib.util
import json
import platform
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any, Iterable, Mapping

import Bio
import numpy as np
import scipy
from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.SASA import ShrakeRupley
from scipy.spatial import cKDTree

SCHEMA_VERSION = "1.0"
RELATIONSHIPS = {"exact_protein", "homolog_assembly", "architecture_analogy", "unresolved"}
CONTACT_CUTOFF_A = 5.0


class InputError(RuntimeError):
    pass


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def read_json(path: str | Path) -> Any:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(f"cannot read JSON {path}: {exc}") from exc


def _load(path: Path):
    try:
        parser = MMCIFParser(QUIET=True) if path.name.lower().endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
        structure = parser.get_structure("assembly", str(path))
        model = next(structure.get_models())
        return structure, model
    except Exception as exc:
        raise InputError(f"cannot parse {path.name}: {type(exc).__name__}: {exc}") from exc


def _expanded_assembly(path: Path, assembly_id: str | None) -> tuple[Any, Any, dict[str, Any]]:
    """Load coordinates and apply a declared mmCIF biological-assembly recipe.

    Deposited biological-assembly files are often already expanded.  We retain
    them when applying the selected Gemmi recipe would not increase the model's
    atom or chain count.  An asymmetric unit with a non-trivial recipe is
    expanded with deterministic ``AddNumber`` copy names (A1, A2, ...).
    """
    structure, model = _load(path)
    result = {
        "operator_backend": "not_required",
        "operator_state": "not_declared",
        "requested_assembly_id": assembly_id,
        "observed_chain_count_before": len(list(model)),
        "observed_atom_count_before": len(list(model.get_atoms())),
        "chain_copies": [
            {"copy_chain_id": str(chain.id), "source_chain_id": str(chain.id), "copy_index": 1}
            for chain in model
        ],
    }
    if not path.name.lower().endswith((".cif", ".mmcif")):
        return structure, model, result
    if importlib.util.find_spec("gemmi") is None:
        metadata = _assembly_metadata(path, assembly_id)
        if metadata["operator_expressions"]:
            raise InputError("mmCIF assembly operators require the declared Gemmi runtime")
        return structure, model, result
    try:
        import gemmi  # type: ignore

        gemmi_structure = gemmi.read_structure(str(path), merge_chain_parts=True)
        assemblies = list(gemmi_structure.assemblies)
        if not assemblies:
            result["operator_backend"] = "gemmi"
            return structure, model, result
        if assembly_id:
            selected = next((item for item in assemblies if str(item.name) == str(assembly_id)), None)
            if selected is None:
                raise InputError(
                    f"assembly_id {assembly_id!r} is not declared; available: "
                    + ", ".join(sorted(str(item.name) for item in assemblies))
                )
        elif len(assemblies) == 1:
            selected = assemblies[0]
        else:
            raise InputError("mmCIF declares multiple biological assemblies; select assembly_id explicitly")
        expanded = gemmi.make_assembly(
            selected, gemmi_structure[0], gemmi.HowToNameCopiedChain.AddNumber
        )
        before_atoms = int(gemmi_structure[0].count_atom_sites())
        after_atoms = int(expanded.count_atom_sites())
        before_chains = len(list(gemmi_structure[0]))
        after_chains = len(list(expanded))
        result.update({
            "operator_backend": "gemmi",
            "selected_assembly_id": str(selected.name),
            "candidate_chain_count": after_chains,
            "candidate_atom_count": after_atoms,
        })
        if after_atoms <= before_atoms and after_chains <= before_chains:
            result["operator_state"] = "already_expanded_or_identity"
            return structure, model, result

        expanded_structure = gemmi.Structure()
        expanded_structure.name = f"assembly_{selected.name}"
        expanded_structure.add_model(expanded)
        with tempfile.TemporaryDirectory(prefix="yauvi-assembly-") as temp_dir:
            expanded_path = Path(temp_dir) / "expanded.cif"
            expanded_structure.make_mmcif_document().write_file(str(expanded_path))
            parsed = MMCIFParser(QUIET=True).get_structure("expanded_assembly", str(expanded_path))
            expanded_model = next(parsed.get_models())
        copy_rows = []
        counters: dict[str, int] = {}
        for chain in expanded_model:
            copy_chain = str(chain.id)
            match = re.match(r"^(.*?)([1-9][0-9]*)$", copy_chain)
            source_chain = match.group(1) if match and match.group(1) else copy_chain
            counters[source_chain] = counters.get(source_chain, 0) + 1
            copy_rows.append({
                "copy_chain_id": copy_chain,
                "source_chain_id": source_chain,
                "copy_index": counters[source_chain],
            })
        result.update({
            "operator_state": "applied",
            "observed_chain_count_after": len(copy_rows),
            "observed_atom_count_after": len(list(expanded_model.get_atoms())),
            "chain_copies": copy_rows,
        })
        return parsed, expanded_model, result
    except InputError:
        raise
    except Exception as exc:
        raise InputError(f"Gemmi assembly expansion failed: {type(exc).__name__}: {exc}") from exc


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    return [str(value)] if isinstance(value, str) else [str(x) for x in value]


def _assembly_metadata(path: Path, assembly_id: str | None) -> dict[str, Any]:
    if not path.name.lower().endswith((".cif", ".mmcif")):
        return {"format": "pdb", "assembly_id": assembly_id, "operator_expressions": [], "asym_ids": []}
    try:
        doc = MMCIF2Dict(str(path))
    except Exception:
        return {"format": "mmcif", "assembly_id": assembly_id, "operator_expressions": [], "asym_ids": []}
    ids = _as_list(doc.get("_pdbx_struct_assembly_gen.assembly_id"))
    ops = _as_list(doc.get("_pdbx_struct_assembly_gen.oper_expression"))
    asym = _as_list(doc.get("_pdbx_struct_assembly_gen.asym_id_list"))
    rows = [
        {"assembly_id": ids[i] if i < len(ids) else "", "operator_expression": ops[i] if i < len(ops) else "",
         "asym_ids": sorted(x.strip() for x in (asym[i] if i < len(asym) else "").split(",") if x.strip())}
        for i in range(max(len(ids), len(ops), len(asym)))
    ]
    selected = [r for r in rows if not assembly_id or r["assembly_id"] == assembly_id]
    return {
        "format": "mmcif", "assembly_id": assembly_id,
        "operator_expressions": sorted({r["operator_expression"] for r in selected if r["operator_expression"]}),
        "asym_ids": sorted({c for r in selected for c in r["asym_ids"]}),
    }


def _heavy_atoms(chain) -> list[Any]:
    return [a for a in chain.get_atoms() if str(getattr(a, "element", "")).upper() not in {"H", "D"}]


def _residue_key(atom) -> tuple[str, int, str, str]:
    residue = atom.get_parent()
    chain = residue.get_parent()
    return str(chain.id), int(residue.id[1]), str(residue.id[2]).strip(), str(residue.resname)


def _contacts(left: list[Any], right: list[Any], cutoff: float) -> tuple[list[dict[str, Any]], str]:
    # In an installed Protein Platform, use the canonical geometry primitive.
    # Standalone packages must not require that sibling, so the independently
    # named cKDTree implementation provides the same Euclidean cutoff semantics.
    try:
        from yauvi_platform.modules.native import interface_geometry as ig  # type: ignore
        la = [ig.Atom(_residue_key(a)[0], _residue_key(a)[1], _residue_key(a)[2],
                      _residue_key(a)[3], a.name, str(a.element), *map(float, a.coord)) for a in left]
        ra = [ig.Atom(_residue_key(a)[0], _residue_key(a)[1], _residue_key(a)[2],
                      _residue_key(a)[3], a.name, str(a.element), *map(float, a.coord)) for a in right]
        pairs = ig.close_pairs(la, ra, cutoff)
        records = [
            {"subject": [a.chain, a.resseq, a.icode, a.resname],
             "partner": [b.chain, b.resseq, b.icode, b.resname], "distance_A": round(float(d), 6)}
            for a, b, d in pairs
        ]
        return records, "yauvi_interface_geometry"
    except (ImportError, AttributeError):
        pass
    if not left or not right:
        return [], "scipy_ckdtree"
    left_xyz = np.asarray([a.coord for a in left], dtype=float)
    right_xyz = np.asarray([a.coord for a in right], dtype=float)
    sparse = cKDTree(left_xyz).sparse_distance_matrix(cKDTree(right_xyz), cutoff, output_type="coo_matrix")
    records = [
        {"subject": list(_residue_key(left[int(i)])), "partner": list(_residue_key(right[int(j)])),
         "distance_A": round(float(d), 6)}
        for i, j, d in sorted(zip(sparse.row, sparse.col, sparse.data), key=lambda item: (int(item[0]), int(item[1])))
    ]
    return records, "scipy_ckdtree"


def _canonical_sasa_frame(model) -> tuple[np.ndarray, np.ndarray]:
    atoms = sorted(
        model.get_atoms(),
        key=lambda atom: (
            str(atom.get_parent().get_parent().id),
            str(atom.get_parent().id[0]), int(atom.get_parent().id[1]), str(atom.get_parent().id[2]),
            str(atom.name), str(atom.altloc),
        ),
    )
    if not atoms:
        raise InputError("SASA calculation requires at least one atom")
    coordinates = np.asarray([atom.coord for atom in atoms], dtype=float)
    origin = coordinates[0]
    offsets = coordinates - origin
    x_axis = next((row for row in offsets if np.linalg.norm(row) > 1e-8), np.asarray([1.0, 0.0, 0.0]))
    x_axis = x_axis / np.linalg.norm(x_axis)
    y_axis = None
    for row in offsets:
        candidate = row - np.dot(row, x_axis) * x_axis
        if np.linalg.norm(candidate) > 1e-8:
            y_axis = candidate / np.linalg.norm(candidate)
            break
    if y_axis is None:
        seed = np.eye(3)[int(np.argmin(np.abs(x_axis)))]
        candidate = seed - np.dot(seed, x_axis) * x_axis
        y_axis = candidate / np.linalg.norm(candidate)
    z_axis = np.cross(x_axis, y_axis)
    z_axis = z_axis / np.linalg.norm(z_axis)
    return origin, np.column_stack((x_axis, y_axis, z_axis))


def _apply_frame(entity, origin: np.ndarray, basis: np.ndarray) -> None:
    for atom in entity.get_atoms():
        atom.coord = np.asarray(atom.coord - origin, dtype=float) @ basis


def _freesasa_document(entity: Any, path: Path) -> dict[str, Any]:
    binary = shutil.which("freesasa")
    if not binary:
        raise InputError("FreeSASA runtime is unavailable")
    chains = [entity] if getattr(entity, "level", "") == "C" else list(entity.get_chains())
    chain_ids = {str(chain.id) for chain in chains}
    use_cif = any(len(chain_id) > 1 for chain_id in chain_ids)
    actual_path = path.with_suffix(".cif") if use_cif else path.with_suffix(".pdb")
    writer = MMCIFIO() if use_cif else PDBIO()
    writer.set_structure(entity)
    writer.save(str(actual_path))
    try:
        arguments = [binary, "--format=json", "--depth=chain", "--n-threads=1"]
        if use_cif:
            arguments.append("--cif")
        arguments.append(str(actual_path))
        process = subprocess.run(
            arguments,
            capture_output=True, text=True, timeout=120, check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise InputError(f"FreeSASA execution failed: {type(exc).__name__}: {exc}") from exc
    if process.returncode != 0:
        detail = (process.stderr or process.stdout).strip().splitlines()
        raise InputError("FreeSASA failed: " + (detail[-1][:240] if detail else "unknown error"))
    try:
        document = json.loads(process.stdout)
        result = document["results"][0]
        # FreeSASA 2.1.x emits ``structure``.  Retain the plural spelling as a
        # compatibility read for older adapter fixtures rather than pinning the
        # parser to one undocumented variation.
        structures = result.get("structure", result.get("structures"))
    except (KeyError, IndexError, TypeError, json.JSONDecodeError) as exc:
        raise InputError("FreeSASA did not return the expected JSON structure") from exc
    if not structures:
        raise InputError("FreeSASA returned no structure result")
    return document


def _freesasa_structures(document: Mapping[str, Any]) -> Any:
    result = document["results"][0]
    return result.get("structure", result.get("structures"))


def _freesasa_total(document: Mapping[str, Any]) -> float:
    try:
        return float(_freesasa_structures(document)[0]["area"]["total"])
    except (KeyError, IndexError, TypeError, ValueError) as exc:
        raise InputError("FreeSASA JSON lacks a total surface area") from exc


def _freesasa_chain(document: Mapping[str, Any], chain_id: str) -> float:
    try:
        chains = _freesasa_structures(document)[0]["chains"]
        row = next(item for item in chains if str(item["label"]) == chain_id)
        return float(row["area"]["total"])
    except (KeyError, StopIteration, IndexError, TypeError, ValueError) as exc:
        raise InputError(f"FreeSASA JSON lacks subject chain {chain_id!r}") from exc


def _sasa(chain, model) -> tuple[float, float, str]:
    if shutil.which("freesasa"):
        with tempfile.TemporaryDirectory(prefix="yauvi-freesasa-") as temp_dir:
            isolated_doc = _freesasa_document(copy.deepcopy(chain), Path(temp_dir) / "isolated.pdb")
            assembly_doc = _freesasa_document(copy.deepcopy(model), Path(temp_dir) / "assembly.pdb")
        return (
            _freesasa_total(isolated_doc),
            _freesasa_chain(assembly_doc, str(chain.id)),
            "freesasa_lee_richards_default_single_thread",
        )
    try:
        isolated = copy.deepcopy(chain)
        assembly = copy.deepcopy(model)
        origin, basis = _canonical_sasa_frame(assembly)
        _apply_frame(isolated, origin, basis)
        _apply_frame(assembly, origin, basis)
        sr = ShrakeRupley(n_points=240)
        sr.compute(isolated, level="R")
        isolated_sasa = sum(float(getattr(r, "sasa", 0.0)) for r in isolated.get_residues())
        sr.compute(assembly, level="R")
        assembly_chain = assembly[chain.id]
        assembly_sasa = sum(float(getattr(r, "sasa", 0.0)) for r in assembly_chain.get_residues())
        return isolated_sasa, assembly_sasa, "biopython_shrake_rupley_240_canonical_frame"
    except Exception as exc:
        raise InputError(f"SASA calculation failed: {type(exc).__name__}: {exc}") from exc


def analyze(
    manifest: Mapping[str, Any], isolated_path: str | Path, assembly_path: str | Path, *,
    subject_chain: str, relationship: str, reference_id: str = "", assembly_id: str | None = None,
    expected_chains: Iterable[str] = (), contact_cutoff_A: float = CONTACT_CUTOFF_A,
) -> dict[str, Any]:
    if relationship not in RELATIONSHIPS:
        raise InputError(f"relationship must be one of {sorted(RELATIONSHIPS)}")
    isolated, assembly = Path(isolated_path), Path(assembly_path)
    if not isolated.is_file() or not assembly.is_file():
        raise InputError("isolated and assembly coordinate files are both required")
    expected_digest = str(manifest.get("coordinate", {}).get("sha256", ""))
    actual_digest = sha256(isolated)
    if not expected_digest or expected_digest != actual_digest:
        raise InputError("StructQC manifest checksum does not match isolated coordinates")
    _isolated_structure, _isolated_model = _load(isolated)
    _assembly_structure, model, operator_record = _expanded_assembly(assembly, assembly_id)
    chains = sorted(str(c.id) for c in model)
    resolved_subject_chain = subject_chain
    if resolved_subject_chain not in chains and operator_record["operator_state"] == "applied":
        matching = [
            row["copy_chain_id"] for row in operator_record["chain_copies"]
            if row["source_chain_id"] == subject_chain
        ]
        if matching:
            resolved_subject_chain = sorted(matching)[0]
    if resolved_subject_chain not in chains:
        raise InputError(f"subject chain {subject_chain!r} not found in assembly; available: {chains}")
    metadata = _assembly_metadata(assembly, assembly_id)
    metadata["operator_application"] = operator_record
    subject = model[resolved_subject_chain]
    partners = [c for c in model if str(c.id) != resolved_subject_chain]
    pairs, geometry_method = _contacts(_heavy_atoms(subject), [a for c in partners for a in _heavy_atoms(c)], contact_cutoff_A)
    minimum: dict[tuple, float] = {}
    for pair in pairs:
        key = (tuple(pair["subject"]), tuple(pair["partner"]))
        minimum[key] = min(minimum.get(key, float("inf")), float(pair["distance_A"]))
    isolated_sasa, assembly_sasa, sasa_method = _sasa(subject, model)
    expected = sorted(set(str(c) for c in expected_chains))
    observed_identities = set(chains) | {
        row["source_chain_id"] for row in operator_record.get("chain_copies", [])
    }
    complete = None if not expected else set(expected).issubset(observed_identities)
    lower_bound = complete is not True
    residue_contacts: dict[tuple, set[str]] = {}
    for subject_key, partner_key in minimum:
        residue_contacts.setdefault(subject_key, set()).add(str(partner_key[0]))
    return {
        "schema_version": SCHEMA_VERSION, "module_id": "assembly_context",
        "subject": manifest["subject"],
        "coordinate_sha256": expected_digest, "assembly_sha256": sha256(assembly),
        "input_sha256": {
            "structure_evidence_manifest": _json_sha256(manifest),
            "isolated": expected_digest,
            "assembly": sha256(assembly),
        },
        "reference": {"id": reference_id, "relationship": relationship, "assembly_id": assembly_id},
        "assembly": {
            "file_name": assembly.name, "chains_observed": chains, "chains_expected": expected,
            "complete": complete, "lower_bound": lower_bound, "metadata": metadata,
        },
        "methods": {"contact_geometry": geometry_method, "contact_cutoff_A": contact_cutoff_A, "sasa": sasa_method},
        "config": {
            "subject_chain": subject_chain,
            "resolved_subject_chain": resolved_subject_chain,
            "relationship": relationship,
            "assembly_id": assembly_id,
            "chains_expected": expected,
            "contact_cutoff_A": contact_cutoff_A,
        },
        "surface": {
            "subject_isolated_sasa_A2": round(isolated_sasa, 6),
            "subject_assembly_sasa_A2": round(assembly_sasa, 6),
            "buried_sasa_A2": round(max(0.0, isolated_sasa - assembly_sasa), 6),
        },
        "interfaces": [
            {"subject_chain": key[0][0], "subject_resseq": key[0][1], "subject_icode": key[0][2],
             "subject_resname": key[0][3], "partner_chain": key[1][0], "partner_resseq": key[1][1],
             "partner_icode": key[1][2], "partner_resname": key[1][3], "minimum_distance_A": round(distance, 6)}
            for key, distance in sorted(minimum.items())
        ],
        "residue_contacts": [
            {"chain_id": key[0], "auth_seq_id": key[1], "insertion_code": key[2], "resname": key[3],
             "partner_chains": sorted(partners)}
            for key, partners in sorted(residue_contacts.items())
        ],
        "limitations": [
            "Assembly geometry from one coordinate state is not native surface exposure.",
            "A homolog assembly transfers architecture, not residue identity.",
            "Incomplete assemblies make contacts and burial lower bounds.",
        ],
    }


def _json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_outputs(out_dir: str | Path, document: Mapping[str, Any]) -> list[Path]:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    summary, table, layer, run = (out / name for name in (
        "ASSEMBLY_CONTEXT.json", "INTERFACES.tsv", "ASSEMBLY_LAYER.json", "RUN_MANIFEST.json"))
    _json(summary, document)
    fields = ["subject_chain", "subject_resseq", "subject_icode", "subject_resname",
              "partner_chain", "partner_resseq", "partner_icode", "partner_resname", "minimum_distance_A"]
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(document["interfaces"])
    _json(layer, {
        "schema_version": SCHEMA_VERSION, "contract_id": "structure_layer_bundle",
        "subject": document["subject"], "coordinate_sha256": document["coordinate_sha256"],
        "layer_id": "assembly_interface",
        "records": [
            {"chain_id": r["chain_id"], "auth_seq_id": r["auth_seq_id"], "insertion_code": r["insertion_code"],
             "metric": "assembly_contact", "state": "observed", "value": len(r["partner_chains"]),
             "detail": "contacts partner chain(s) " + ",".join(r["partner_chains"]),
             "source_digest": document["assembly_sha256"],
             "evidence_class": document["reference"]["relationship"]}
            for r in document["residue_contacts"]
        ],
    })
    _json(run, {"schema_version": SCHEMA_VERSION, "module_id": "assembly_context", "version": "0.1.0",
                "input_sha256": document["input_sha256"],
                "parameters": document["config"],
                "runtime_versions": {"python": platform.python_version(), "biopython": Bio.__version__,
                                     "numpy": np.__version__, "scipy": scipy.__version__},
                "optional_runtimes": {
                    "gemmi": "available_not_required" if importlib.util.find_spec("gemmi") else "not_available",
                    "freesasa": (
                        "available_invoked" if str(document["methods"]["sasa"]).startswith("freesasa_")
                        else ("available_not_invoked" if shutil.which("freesasa") else "not_available")
                    ),
                },
                "outputs": [summary.name, table.name, layer.name],
                "missing_evidence": (
                    (["complete_assembly"] if document["assembly"]["lower_bound"] else [])
                    + (["publication_grade_freesasa"] if not str(document["methods"]["sasa"]).startswith("freesasa_") else [])
                )})
    return [summary, table, layer, run]
