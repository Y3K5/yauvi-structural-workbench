from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import platform
from pathlib import Path
import re
import shutil
from typing import Any, Mapping
import xml.etree.ElementTree as ET

import Bio
import numpy as np
from Bio import Align
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

SCHEMA_VERSION = "1.0"
PROVENANCE = {"experimental", "predicted", "unknown"}
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}
BACKBONE = {"N", "CA", "C", "O"}


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


def read_json(path: str | Path | None) -> Any:
    if path is None:
        return None
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(f"cannot read JSON {path}: {exc}") from exc


def read_validation_report(path: str | Path | None) -> dict[str, Any] | None:
    """Import a checksum-bound validation summary without reinterpreting its scores."""
    if path is None:
        return None
    source = Path(path)
    if not source.is_file():
        raise InputError(f"validation report not found: {source}")
    raw = source.read_bytes()
    suffix = source.suffix.lower()
    if suffix not in {".json", ".xml"}:
        suffix = ".xml" if raw.lstrip().startswith(b"<") else ".json"
    candidates: dict[str, Any] = {}
    if suffix == ".json":
        try:
            document = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise InputError(f"validation JSON cannot be parsed: {exc}") from exc

        def walk(value: Any, prefix: str = "") -> None:
            if isinstance(value, Mapping):
                for key, child in value.items():
                    walk(child, f"{prefix}_{key}" if prefix else str(key))
            elif isinstance(value, (int, float)) and not isinstance(value, bool):
                candidates[prefix.lower()] = value
        walk(document)
    elif suffix == ".xml":
        try:
            root = ET.fromstring(raw)
        except ET.ParseError as exc:
            raise InputError(f"validation XML cannot be parsed: {exc}") from exc
        for element in root.iter():
            tag = element.tag.rsplit("}", 1)[-1].lower()
            if element.text and element.text.strip():
                try: candidates[tag] = float(element.text.strip())
                except ValueError: pass
            for key, value in element.attrib.items():
                try: candidates[f"{tag}_{key.lower()}"] = float(value)
                except ValueError: pass
    else:
        raise InputError("validation report must be JSON or XML")
    normalized = {
        re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_"): value
        for key, value in candidates.items()
    }
    aliases = {
        "clashscore": ("clashscore", "clash_score"),
        "ramachandran_outliers_percent": (
            "percent_rama_outliers", "rama_outlier_percent", "ramachandran_outlier_percent",
        ),
        "rotamer_outliers_percent": ("percent_rota_outliers", "rotamer_outlier_percent"),
        "r_free": ("pdb_rfree", "r_free", "rfree"),
        "resolution_angstrom": ("pdb_resolution", "resolution"),
        "rsrz_outliers_percent": ("percent_rsrz_outliers", "rsrz_outlier_percent"),
    }
    metrics: dict[str, float] = {}
    for metric, needles in aliases.items():
        # Attribute names such as ``absolute-percentile-clashscore`` contain the
        # word clashscore but are not the metric. Prefer an exact suffix and the
        # shortest qualified key (e.g. Entry.clashscore) to avoid importing a
        # percentile in place of the raw value.
        matches = [
            (key, value) for key, value in normalized.items()
            if any(key == needle or key.endswith("_" + needle) for needle in needles)
        ]
        if matches:
            _key, value = min(matches, key=lambda item: (len(item[0]), item[0]))
            metrics[metric] = float(value)
    return {
        "state": "imported", "format": suffix.removeprefix("."),
        "file_name": source.name, "sha256": hashlib.sha256(raw).hexdigest(),
        "metrics": metrics,
        "limitations": ["Imported values retain the producing validator's definitions and are not recomputed by StructQC."],
    }


def read_fasta(path: str | Path | None) -> tuple[str | None, str | None]:
    if path is None:
        return None, None
    identifier: str | None = None
    chunks: list[str] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if identifier is not None:
                raise InputError("reference FASTA must contain exactly one sequence")
            header = line[1:].strip()
            parts = header.split("|")
            identifier = parts[1] if len(parts) >= 3 else header.split()[0]
        elif identifier is not None:
            chunks.append(line.strip())
    sequence = "".join(chunks).upper()
    if not identifier or not sequence:
        raise InputError("reference FASTA contains no sequence")
    return identifier, sequence


def _parser(path: Path):
    lower = path.name.lower()
    if lower.endswith((".cif", ".mmcif")):
        return MMCIFParser(QUIET=True).get_structure("query", str(path)), "mmcif"
    return PDBParser(QUIET=True).get_structure("query", str(path)), "pdb"


def _gemmi_metadata(path: Path) -> dict[str, Any]:
    if importlib.util.find_spec("gemmi") is None:
        return {"state": "runtime_missing", "backend": "biopython", "residue_identity": {}}
    try:
        import gemmi  # type: ignore
        parsed = gemmi.read_structure(str(path))
        if len(parsed) == 0:
            raise ValueError("Gemmi found no coordinate models")
        identities: dict[tuple[str, int, str], dict[str, Any]] = {}
        if path.name.lower().endswith((".cif", ".mmcif")):
            block = gemmi.cif.read_file(str(path)).sole_block()
            columns = {name: list(block.find_values(name)) for name in (
                "_atom_site.auth_asym_id", "_atom_site.auth_seq_id", "_atom_site.pdbx_PDB_ins_code",
                "_atom_site.label_asym_id", "_atom_site.label_seq_id", "_atom_site.label_entity_id",
            )}
            count = max((len(values) for values in columns.values()), default=0)
            for index in range(count):
                def value(name: str) -> str:
                    values = columns[name]; return str(values[index]) if index < len(values) else ""
                try: auth_seq = int(float(value("_atom_site.auth_seq_id")))
                except ValueError: continue
                insertion = value("_atom_site.pdbx_PDB_ins_code").replace("?", "").replace(".", "")
                label_raw = value("_atom_site.label_seq_id")
                try: label_seq: int | None = int(float(label_raw))
                except ValueError: label_seq = None
                key = (value("_atom_site.auth_asym_id"), auth_seq, insertion)
                identities.setdefault(key, {
                    "label_asym_id": value("_atom_site.label_asym_id"),
                    "label_seq_id": label_seq,
                    "entity_id": value("_atom_site.label_entity_id"),
                })
        return {
            "state": "validated", "backend": "gemmi", "gemmi_version": getattr(gemmi, "__version__", "unknown"),
            "model_count": len(parsed), "residue_identity": identities,
        }
    except Exception as exc:
        raise InputError(f"Gemmi coordinate validation failed: {type(exc).__name__}: {exc}") from exc


def _assemblies(path: Path, fmt: str) -> list[dict[str, Any]]:
    if fmt != "mmcif":
        return []
    try:
        doc = MMCIF2Dict(str(path))
    except Exception:
        return []
    ids = doc.get("_pdbx_struct_assembly.id", [])
    details = doc.get("_pdbx_struct_assembly.details", [])
    if isinstance(ids, str):
        ids = [ids]
    if isinstance(details, str):
        details = [details]
    return [
        {"assembly_id": str(aid), "detail": str(details[i]) if i < len(details) else ""}
        for i, aid in enumerate(ids)
    ]


def _prediction_hint(path: Path) -> str | None:
    upper = path.read_text(encoding="utf-8", errors="ignore")[:8000].upper()
    for name in ("ALPHAFOLD", "COLABFOLD", "OPENFOLD", "CHAI", "BOLTZ", "ESMFOLD", "ROSETTAFOLD"):
        if name in upper:
            return f"header mentions {name}"
    if path.name.upper().startswith("AF-"):
        return "filename follows the AlphaFold DB convention"
    return None


def _provenance(raw: Mapping[str, Any] | None, structure: Path) -> tuple[dict[str, Any], list[str]]:
    warnings: list[str] = []
    document = dict(raw or {})
    cls = str(document.get("class", "unknown")).strip().lower()
    if cls not in PROVENANCE:
        raise InputError(f"provenance.class must be one of {sorted(PROVENANCE)}")
    if cls == "experimental" and not str(document.get("method", "")).strip():
        raise InputError("experimental provenance requires method")
    hint = _prediction_hint(structure)
    if cls == "unknown" and hint:
        warnings.append(f"prediction hint observed but not promoted to declared provenance: {hint}")
    if cls == "experimental" and hint:
        warnings.append(f"declared experimental provenance conflicts with prediction hint: {hint}")
    return {
        "class": cls,
        "method": str(document.get("method", "")),
        "source_id": str(document.get("source_id", "")),
        "confidence_encoding": str(document.get("confidence_encoding", "")),
        **({"resolution_angstrom": float(document["resolution_angstrom"])}
           if document.get("resolution_angstrom") not in (None, "") else {}),
    }, warnings


def _sequence_map(observed: str, reference: str | None) -> tuple[dict[int, int], dict[str, Any]]:
    if reference is None:
        return {}, {"state": "unevaluated", "detail": "no reference sequence supplied"}
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference, observed)[0]
    mapping: dict[int, int] = {}
    for (rs, re), (os, oe) in zip(alignment.aligned[0], alignment.aligned[1]):
        for offset in range(min(re - rs, oe - os)):
            mapping[int(os + offset)] = int(rs + offset + 1)
    identity = sum(
        1 for oi, ri in mapping.items() if oi < len(observed) and reference[ri - 1] == observed[oi]
    )
    return mapping, {
        "state": "evaluated",
        "reference_length": len(reference),
        "coordinate_residues": len(observed),
        "mapped_residues": len(mapping),
        "identity_fraction": round(identity / len(mapping), 6) if mapping else 0.0,
        "coverage_fraction": round(len(set(mapping.values())) / len(reference), 6),
    }


def _pae_summary(raw: Any) -> dict[str, Any]:
    if raw is None:
        return {"state": "unevaluated", "detail": "no PAE supplied"}
    matrix = raw
    if isinstance(raw, list) and raw and isinstance(raw[0], dict):
        matrix = raw[0].get("predicted_aligned_error")
    elif isinstance(raw, dict):
        matrix = raw.get("predicted_aligned_error", raw.get("pae"))
    try:
        array = np.asarray(matrix, dtype=float)
    except (TypeError, ValueError) as exc:
        raise InputError(f"PAE is not numeric: {exc}") from exc
    if array.ndim != 2 or array.shape[0] != array.shape[1] or array.size == 0:
        raise InputError("PAE must be a non-empty square matrix")
    return {
        "state": "evaluated",
        "size": int(array.shape[0]),
        "mean_angstrom": round(float(array.mean()), 6),
        "max_angstrom": round(float(array.max()), 6),
    }


def analyze(
    structure_path: str | Path,
    *,
    subject_id: str | None = None,
    provenance: Mapping[str, Any] | None = None,
    reference_sequence: str | None = None,
    reference_id: str | None = None,
    pae: Any = None,
    validation_report: Mapping[str, Any] | None = None,
    model_index: int = 0,
    chain: str | None = None,
) -> dict[str, Any]:
    path = Path(structure_path)
    if not path.is_file():
        raise InputError(f"structure not found: {path}")
    gemmi_metadata = _gemmi_metadata(path)
    try:
        structure, fmt = _parser(path)
    except Exception as exc:
        raise InputError(f"structure cannot be parsed: {type(exc).__name__}: {exc}") from exc
    models = list(structure.get_models())
    if not models:
        raise InputError("structure contains no models")
    if model_index < 0 or model_index >= len(models):
        raise InputError(f"model index {model_index} is outside 0..{len(models)-1}")

    prov, warnings = _provenance(provenance, path)
    model = models[model_index]
    available_chains = sorted(str(c.id) for c in model)
    if chain and chain not in available_chains:
        raise InputError(f"chain {chain!r} not found; available: {available_chains}")

    residues: list[dict[str, Any]] = []
    chain_summaries: list[dict[str, Any]] = []
    sequence_offset = 0
    observed_all = ""
    chain_entries: list[tuple[str, list[Any], str]] = []
    for ch in sorted(model, key=lambda item: str(item.id)):
        if chain and str(ch.id) != chain:
            continue
        amino = [r for r in ch if r.id[0] == " "]
        observed = "".join(AA3.get(r.resname.upper(), "X") for r in amino)
        observed_all += observed
        chain_entries.append((str(ch.id), amino, observed))

    mapping, completeness = _sequence_map(observed_all, reference_sequence)
    for chain_id, amino, observed in chain_entries:
        missing_backbone = 0
        nonstandard = 0
        chain_breaks = 0
        previous_c = None
        for local_index, residue in enumerate(amino):
            atoms = list(residue.get_atoms())
            names = {a.name for a in atoms}
            missing = sorted(BACKBONE - names)
            missing_backbone += bool(missing)
            aa = AA3.get(residue.resname.upper(), "X")
            nonstandard += aa == "X"
            ca = residue["CA"] if "CA" in residue else None
            if previous_c is not None and "N" in residue:
                distance = float(np.linalg.norm(previous_c.coord - residue["N"].coord))
                if distance > 2.5:
                    chain_breaks += 1
            previous_c = residue["C"] if "C" in residue else None
            b_factors = [float(a.bfactor) for a in atoms]
            altlocs = sorted({str(a.altloc).strip() for a in atoms if str(a.altloc).strip()})
            global_index = sequence_offset + local_index
            key = residue.id
            record = {
                "chain_id": chain_id,
                "auth_seq_id": int(key[1]),
                "insertion_code": str(key[2]).strip(),
                "resname": residue.resname,
                "one_letter": aa,
                "sequence_index": mapping.get(global_index),
                "missing_backbone_atoms": missing,
                "altlocs": altlocs,
                "mean_b_factor": round(float(np.mean(b_factors)), 6) if b_factors else None,
                "ca_xyz": [round(float(v), 6) for v in ca.coord] if ca is not None else None,
            }
            identity = gemmi_metadata["residue_identity"].get(
                (chain_id, int(key[1]), str(key[2]).strip()), {}
            )
            record.update({
                "label_asym_id": identity.get("label_asym_id", chain_id if fmt == "pdb" else ""),
                "label_seq_id": identity.get("label_seq_id"),
                "entity_id": identity.get("entity_id", ""),
            })
            if prov["class"] == "predicted" and (
                prov.get("confidence_encoding") == "plddt_in_bfactor"
                or any(name in prov.get("method", "").upper() for name in ("ALPHAFOLD", "COLABFOLD", "OPENFOLD"))
            ):
                record["plddt"] = record["mean_b_factor"]
            residues.append(record)
        chain_summaries.append({
            "chain_id": chain_id,
            "residues": len(amino),
            "sequence": observed,
            "missing_backbone_residues": missing_backbone,
            "nonstandard_residues": nonstandard,
            "chain_breaks": chain_breaks,
        })
        sequence_offset += len(amino)

    if not residues:
        raise InputError("selected model/chain contains no amino-acid residues")
    return {
        "schema_version": SCHEMA_VERSION,
        "module_id": "structure_quality",
        "subject": {"id": subject_id or reference_id or path.stem},
        "coordinate": {
            "file_name": path.name,
            "sha256": sha256(path),
            "format": fmt,
            "selected_model": model_index,
            "selected_chain": chain,
            "model_count": len(models),
            "chains": available_chains,
            "assemblies": _assemblies(path, fmt),
            "parser": {
                "coordinate_analysis": "biopython",
                "gemmi_validation": gemmi_metadata["state"],
                **({"gemmi_version": gemmi_metadata.get("gemmi_version")} if gemmi_metadata.get("gemmi_version") else {}),
            },
        },
        "input_sha256": {
            "structure": sha256(path),
            **({"provenance": _json_sha256(dict(provenance))} if provenance is not None else {}),
            **({"reference_sequence": hashlib.sha256(reference_sequence.encode("ascii")).hexdigest()}
               if reference_sequence is not None else {}),
            **({"pae": _json_sha256(pae)} if pae is not None else {}),
            **({"external_validation": str(validation_report["sha256"])}
               if validation_report is not None and validation_report.get("sha256") else {}),
        },
        "provenance": prov,
        "reference": {"id": reference_id or "", "sequence_supplied": reference_sequence is not None},
        "completeness": completeness,
        "pae": _pae_summary(pae),
        "external_validation": dict(validation_report) if validation_report is not None else {
            "state": "missing", "metrics": {},
            "limitations": ["No wwPDB, MolProbity, or Phenix validation report was supplied."],
        },
        "chain_summaries": chain_summaries,
        "residues": residues,
        "warnings": warnings,
        "limitations": [
            "Coordinate quality does not establish native conformation or biological function.",
            "Unknown provenance remains unknown; filename/header hints cannot certify experiment.",
            "Completeness is unevaluated without a reference sequence.",
        ],
    }


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_outputs(out_dir: str | Path, document: Mapping[str, Any]) -> list[Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    evidence = out / "STRUCTURE_EVIDENCE.json"
    table = out / "RESIDUE_QUALITY.tsv"
    layer = out / "STRUCTURE_LAYER.json"
    run = out / "RUN_MANIFEST.json"
    _write_json(evidence, document)
    fields = ["chain_id", "auth_seq_id", "insertion_code", "label_asym_id", "label_seq_id", "entity_id",
              "resname", "one_letter", "sequence_index",
              "missing_backbone_atoms", "altlocs", "mean_b_factor", "plddt", "ca_xyz"]
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for item in document["residues"]:
            row = dict(item)
            for key in ("missing_backbone_atoms", "altlocs", "ca_xyz"):
                row[key] = ";".join(map(str, row.get(key) or []))
            writer.writerow(row)
    layer_doc = {
        "schema_version": SCHEMA_VERSION,
        "contract_id": "structure_layer_bundle",
        "subject": document["subject"],
        "coordinate_sha256": document["coordinate"]["sha256"],
        "layer_id": "structure_quality",
        "records": [
            {
                "chain_id": r["chain_id"], "auth_seq_id": r["auth_seq_id"],
                "insertion_code": r["insertion_code"], "metric": "coordinate_quality",
                "state": "limited" if r["missing_backbone_atoms"] else "observed",
                "value": r.get("plddt"),
                "detail": ("missing " + ",".join(r["missing_backbone_atoms"])) if r["missing_backbone_atoms"] else "coordinates observed",
                "source_digest": document["coordinate"]["sha256"],
                "evidence_class": document["provenance"]["class"],
            }
            for r in document["residues"]
        ],
    }
    _write_json(layer, layer_doc)
    _write_json(run, {
        "schema_version": SCHEMA_VERSION,
        "module_id": "structure_quality",
        "version": "0.1.0",
        "input_sha256": document["input_sha256"],
        "parameters": {
            "selected_model": document["coordinate"]["selected_model"],
            "selected_chain": document["coordinate"]["selected_chain"],
            "provenance_class": document["provenance"]["class"],
            "confidence_encoding": document["provenance"]["confidence_encoding"],
            "reference_sequence_supplied": document["reference"]["sequence_supplied"],
            "pae_state": document["pae"]["state"],
            "external_validation_state": document["external_validation"]["state"],
        },
        "runtime_versions": {
            "python": platform.python_version(),
            "biopython": Bio.__version__,
            "numpy": np.__version__,
        },
        "optional_runtimes": {
            "gemmi": "available_invoked" if importlib.util.find_spec("gemmi") else "not_available",
            "mkdssp": "available_not_invoked" if shutil.which("mkdssp") else "not_available",
        },
        "outputs": [p.name for p in (evidence, table, layer)],
        "missing_evidence": [
            name for name, state in (("provenance", document["provenance"]["class"]),
                                     ("reference_sequence", document["completeness"]["state"]),
                                     ("community_geometry_validation", document["external_validation"]["state"]),
                                     ("gemmi_coordinate_validation", document["coordinate"]["parser"]["gemmi_validation"]))
            if state in {"unknown", "unevaluated"}
            or (name == "community_geometry_validation" and state == "missing")
            or (name == "gemmi_coordinate_validation" and state == "runtime_missing")
        ],
    })
    return [evidence, table, layer, run]
