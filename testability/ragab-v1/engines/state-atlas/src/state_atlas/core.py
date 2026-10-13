from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

import Bio
import numpy as np
import scipy
from Bio import Align
from Bio.PDB import MMCIFParser, PDBParser
from scipy.cluster.hierarchy import fcluster, linkage

SCHEMA_VERSION = "1.1"
REFERENCE_SET_V2 = "2.0"
DEFAULT_SELECTION = "protein and name CA"
ABL_DOMAIN_START = 242
ABL_DOMAIN_END = 495
ABL_MAX_RMSD_A = 2.5
ABL_MIN_MARGIN_A = 0.25
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
STATES = {"active", "inactive"}
LABELS = {"active_like", "inactive_like", "mixed", "unresolved"}
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O",
}


class InputError(RuntimeError):
    pass


def sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def read_json(path: str | Path) -> Any:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(f"cannot read JSON {path}: {exc}") from exc


def _exact_map_rows(rows: Any, *, label: str, errors: list[str]) -> dict[int, tuple[str, int, str]]:
    if not isinstance(rows, list):
        errors.append(f"{label}: residue mappings must be a list")
        return {}
    exact: dict[int, tuple[str, int, str]] = {}
    seen: set[int] = set()
    for row_number, row in enumerate(rows, 1):
        if not isinstance(row, Mapping):
            errors.append(f"{label}: mapping row {row_number} must be an object")
            continue
        try:
            position = int(row["uniprot_position"])
        except (KeyError, TypeError, ValueError):
            errors.append(f"{label}: mapping row {row_number} has no integer UniProt position")
            continue
        if position in seen:
            errors.append(f"{label}: duplicate UniProt position {position}")
            continue
        seen.add(position)
        state = str(row.get("mapping_state", "exact"))
        if state != "exact":
            continue
        chain = str(row.get("chain_id", ""))
        try:
            auth_seq_id = int(row["auth_seq_id"])
        except (KeyError, TypeError, ValueError):
            errors.append(f"{label}: exact position {position} has no integer auth_seq_id")
            continue
        if not chain:
            errors.append(f"{label}: exact position {position} has no chain_id")
            continue
        exact[position] = (chain, auth_seq_id, str(row.get("insertion_code", "")))
    return exact


def _mapping_states(rows: Any) -> dict[int, str]:
    """Retain declared non-exact states for per-comparison diagnostics."""
    states: dict[int, str] = {}
    if not isinstance(rows, list):
        return states
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        try:
            position = int(row["uniprot_position"])
        except (KeyError, TypeError, ValueError):
            continue
        states[position] = str(row.get("mapping_state", "exact"))
    return states


def validate_alignment_map(
    document: Mapping[str, Any], reference_set: Mapping[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    """Validate the exact residue-equivalence map required by Reference Set v2."""
    errors: list[str] = []
    mask = reference_set.get("alignment_mask", {})
    if not isinstance(mask, Mapping):
        return ["alignment_mask must be an object"], {}
    if str(document.get("coordinate_system", "")) != "uniprot":
        errors.append("alignment map coordinate_system must be uniprot")
    source = document.get("source", {})
    if not isinstance(source, Mapping) or not str(source.get("id", "")).strip() or not str(source.get("citation", "")).strip():
        errors.append("alignment map source requires id and citation")
    if not SHA256_PATTERN.fullmatch(str(source.get("sha256", ""))):
        errors.append("alignment map source.sha256 must be a SHA-256 digest")
    try:
        start = int(mask["uniprot_start"]); end = int(mask["uniprot_end"])
        minimum = float(mask.get("minimum_coverage", 0.9))
    except (KeyError, TypeError, ValueError):
        return [*errors, "alignment_mask requires numeric uniprot_start, uniprot_end, and minimum_coverage"], {}
    if start > end:
        errors.append("alignment_mask start must not exceed end")
    if minimum < 0.9 or minimum > 1.0:
        errors.append("Reference Set v2 minimum_coverage must be between 0.9 and 1.0")
    map_domain = document.get("domain", {})
    try:
        map_start = int(map_domain.get("uniprot_start", -1)) if isinstance(map_domain, Mapping) else -1
        map_end = int(map_domain.get("uniprot_end", -1)) if isinstance(map_domain, Mapping) else -1
    except (TypeError, ValueError):
        map_start, map_end = -1, -1
    if map_start != start or map_end != end:
        errors.append("alignment map domain does not match the reference-set alignment mask")
    expected = max(end - start + 1, 1)
    query = _exact_map_rows(document.get("query"), label="query", errors=errors)
    query_states = _mapping_states(document.get("query"))
    references_raw = document.get("references", {})
    if not isinstance(references_raw, Mapping):
        errors.append("alignment map references must be an object keyed by reference_id")
        references_raw = {}
    reference_maps: dict[str, dict[int, tuple[str, int, str]]] = {}
    reference_states: dict[str, dict[int, str]] = {}
    for reference in reference_set.get("references", []):
        reference_id = str(reference.get("reference_id", ""))
        reference_maps[reference_id] = _exact_map_rows(
            references_raw.get(reference_id), label=f"reference {reference_id}", errors=errors,
        )
        reference_states[reference_id] = _mapping_states(references_raw.get(reference_id))
    for label, mapping in [("query", query), *[(f"reference {key}", value) for key, value in reference_maps.items()]]:
        within = {position for position in mapping if start <= position <= end}
        coverage = len(within) / expected
        if coverage < minimum:
            errors.append(f"{label}: exact mapped coverage {coverage:.3f} is below {minimum:.3f}")
    return errors, {
        "start": start, "end": end, "minimum_coverage": minimum,
        "query": query, "references": reference_maps,
        "query_states": query_states, "reference_states": reference_states,
    }


def validate_reference_set(
    document: Mapping[str, Any], *, base: str | Path = ".",
    alignment_map: Mapping[str, Any] | None = None, alignment_map_sha256: str | None = None,
) -> list[str]:
    errors: list[str] = []
    refs = document.get("references", [])
    if not isinstance(refs, list) or not refs:
        return ["reference set has no references"]
    seen: set[str] = set(); experimental_states: set[str] = set(); pdb_entries: set[str] = set()
    is_v2 = str(document.get("schema_version", "1.0")) == REFERENCE_SET_V2
    if is_v2:
        if not str(document.get("qualification_scope", "")).strip(): errors.append("Reference Set v2 qualification_scope is required")
        if not str(document.get("subject_family", "")).strip(): errors.append("Reference Set v2 subject_family is required")
        expected_map = str(document.get("alignment_map_sha256", ""))
        if not SHA256_PATTERN.fullmatch(expected_map): errors.append("Reference Set v2 alignment_map_sha256 is required")
        elif alignment_map_sha256 != expected_map: errors.append("alignment map checksum does not match Reference Set v2")
        if alignment_map is None:
            errors.append("Reference Set v2 requires an alignment map")
        if str(document.get("qualification_scope", "")) == "abl_family":
            if str(document.get("subject_family", "")) != "ABL1":
                errors.append("Mark 1 ABL-family Reference Set v2 subject_family must be ABL1")
            mask = document.get("alignment_mask", {})
            if not isinstance(mask, Mapping) or (
                mask.get("uniprot_start") != ABL_DOMAIN_START
                or mask.get("uniprot_end") != ABL_DOMAIN_END
            ):
                errors.append(
                    f"Mark 1 ABL-family alignment mask must be UniProt ABL1 "
                    f"{ABL_DOMAIN_START}-{ABL_DOMAIN_END}"
                )
    for ref in refs:
        rid = str(ref.get("reference_id", ""))
        state = str(ref.get("state", ""))
        provenance = ref.get("provenance", {})
        if not rid: errors.append("reference_id is required")
        elif rid in seen: errors.append(f"duplicate reference_id: {rid}")
        seen.add(rid)
        if state not in STATES: errors.append(f"{rid}: state must be active or inactive")
        if str(provenance.get("class", "")) != "experimental":
            errors.append(f"{rid}: only experimental references can carry known state")
        else: experimental_states.add(state)
        if not str(provenance.get("method", "")).strip(): errors.append(f"{rid}: experimental method is required")
        basis = str(ref.get("state_evidence", {}).get("basis", "")).strip()
        if not basis or basis.lower() in {state, f"{state} state"}: errors.append(f"{rid}: state evidence must state its basis")
        location = Path(base) / str(ref.get("structure", ""))
        if not ref.get("structure") or not location.is_file(): errors.append(f"{rid}: reference structure is missing")
        if is_v2:
            chain = str(ref.get("chain", ""))
            pdb_entry = str(ref.get("pdb_entry_id", "")).upper()
            if not chain: errors.append(f"{rid}: exact reference chain is required")
            if not pdb_entry: errors.append(f"{rid}: pdb_entry_id is required")
            elif pdb_entry in pdb_entries: errors.append(f"{rid}: duplicate PDB entry in reference set: {pdb_entry}")
            pdb_entries.add(pdb_entry)
            expected_structure = str(ref.get("structure_sha256", ""))
            if not SHA256_PATTERN.fullmatch(expected_structure):
                errors.append(f"{rid}: structure_sha256 is required")
            elif location.is_file() and sha256(location) != expected_structure:
                errors.append(f"{rid}: reference structure checksum mismatch")
            citation = str(ref.get("state_evidence", {}).get("citation", "")).strip()
            if not citation: errors.append(f"{rid}: state evidence citation is required")
    if experimental_states != STATES:
        errors.append("reference set requires at least one experimental active and inactive reference")
    if is_v2:
        state_counts = {
            state: sum(str(ref.get("state", "")) == state for ref in refs)
            for state in STATES
        }
        for state in sorted(STATES):
            if state_counts[state] < 2:
                errors.append(f"Reference Set v2 requires multiple {state} references (at least two)")
    rules = document.get("decision_rules", {})
    try:
        if float(rules.get("max_rmsd_A", 0)) <= 0: errors.append("max_rmsd_A must be positive")
        if float(rules.get("min_margin_A", -1)) < 0: errors.append("min_margin_A must be non-negative")
        if is_v2 and str(document.get("qualification_scope", "")) == "abl_family":
            if float(rules.get("max_rmsd_A", 0)) != ABL_MAX_RMSD_A:
                errors.append(f"Mark 1 ABL max_rmsd_A must remain frozen at {ABL_MAX_RMSD_A}")
            if float(rules.get("min_margin_A", -1)) != ABL_MIN_MARGIN_A:
                errors.append(f"Mark 1 ABL min_margin_A must remain frozen at {ABL_MIN_MARGIN_A}")
    except (TypeError, ValueError):
        errors.append("decision-rule thresholds must be numeric")
    if is_v2 and alignment_map is not None:
        map_errors, _details = validate_alignment_map(alignment_map, document)
        errors.extend(map_errors)
    return errors


def _bio_models(path: Path, chain: str | None = None) -> tuple[list[np.ndarray], list[tuple[str, int, str]], str]:
    try:
        parser = MMCIFParser(QUIET=True) if path.name.lower().endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
        structure = parser.get_structure("state", str(path))
    except Exception as exc:
        raise InputError(f"cannot parse structure {path.name}: {type(exc).__name__}: {exc}") from exc
    frames: list[np.ndarray] = []; keys: list[tuple[str, int, str]] = []; sequence = ""
    for model_index, model in enumerate(structure):
        coords = []; model_keys = []; letters = []
        for ch in model:
            if chain and str(ch.id) != chain: continue
            for residue in ch:
                if residue.id[0] != " " or "CA" not in residue: continue
                coords.append(np.asarray(residue["CA"].coord, dtype=float))
                model_keys.append((str(ch.id), int(residue.id[1]), str(residue.id[2]).strip()))
                letters.append(AA3.get(residue.resname.upper(), "X"))
        if not coords: continue
        if model_index == 0:
            keys, sequence = model_keys, "".join(letters)
        elif model_keys != keys:
            raise InputError("all static ensemble models must carry the same ordered CA residue keys")
        frames.append(np.asarray(coords, dtype=float))
    if not frames: raise InputError(f"no CA coordinates found in {path.name}")
    return frames, keys, sequence


def _trajectory(topology: Path, trajectory: Path, *, selection: str, stride: int, pbc: str):
    if stride < 1: raise InputError("trajectory stride must be a positive integer")
    if pbc not in {"none", "unwrap"}: raise InputError("trajectory pbc must be explicitly none or unwrap")
    try:
        import MDAnalysis as mda  # type: ignore
    except ImportError as exc:
        raise InputError("trajectory analysis requires pip install 'yauvi-state-atlas[md]'") from exc
    try:
        universe = mda.Universe(str(topology), str(trajectory))
        atoms = universe.select_atoms(selection)
        if not len(atoms): raise InputError(f"MDAnalysis selection matched no atoms: {selection}")
        if pbc == "unwrap":
            try:
                from MDAnalysis.transformations import unwrap
                universe.trajectory.add_transformations(unwrap(universe.atoms))
            except Exception as exc:
                raise InputError(f"trajectory unwrap requested but unavailable: {exc}") from exc
        frames = [np.asarray(atoms.positions, dtype=float).copy() for _ in universe.trajectory[::stride]]
        keys = [(str(a.chainID or a.segid or "_"), int(a.resid), str(getattr(a, "icode", "") or "")) for a in atoms]
        sequence = "".join(AA3.get(str(a.resname).upper(), "X") for a in atoms)
        if not frames: raise InputError("trajectory selection produced no frames")
        return frames, keys, sequence
    except InputError:
        raise
    except Exception as exc:
        raise InputError(f"trajectory cannot be read: {type(exc).__name__}: {exc}") from exc


def _sequence_indices(query: str, reference: str) -> tuple[np.ndarray, np.ndarray]:
    aligner = Align.PairwiseAligner(); aligner.mode = "global"; aligner.match_score = 2
    aligner.mismatch_score = -1; aligner.open_gap_score = -5; aligner.extend_gap_score = -0.5
    alignment = aligner.align(query, reference)[0]
    q: list[int] = []; r: list[int] = []
    for (qs, qe), (rs, re) in zip(alignment.aligned[0], alignment.aligned[1]):
        n = min(qe - qs, re - rs); q.extend(range(int(qs), int(qs + n))); r.extend(range(int(rs), int(rs + n)))
    if len(q) < 3: raise InputError("query/reference mapping has fewer than three shared residues")
    return np.asarray(q, dtype=int), np.asarray(r, dtype=int)


def _coordinate_index(keys: list[tuple[str, int, str]], *, label: str) -> dict[tuple[str, int, str], int]:
    index: dict[tuple[str, int, str], int] = {}
    for position, key in enumerate(keys):
        if key in index:
            raise InputError(f"{label} has duplicate coordinate residue {key[0]}:{key[1]}{key[2]}")
        index[key] = position
    return index


def _mapped_indices(
    query_keys: list[tuple[str, int, str]], reference_keys: list[tuple[str, int, str]],
    query_map: Mapping[int, tuple[str, int, str]], reference_map: Mapping[int, tuple[str, int, str]],
    *, start: int, end: int, minimum_coverage: float, reference_id: str,
    query_states: Mapping[int, str] | None = None,
    reference_states: Mapping[int, str] | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    query_index = _coordinate_index(query_keys, label="query")
    reference_index = _coordinate_index(reference_keys, label=f"reference {reference_id}")
    declared = list(range(start, end + 1))
    shared = [position for position in declared if position in query_map and position in reference_map]
    missing_coordinate: list[int] = []
    query_positions: list[int] = []
    reference_positions: list[int] = []
    for uniprot_position in shared:
        query_key = query_map[uniprot_position]
        reference_key = reference_map[uniprot_position]
        if query_key not in query_index or reference_key not in reference_index:
            missing_coordinate.append(uniprot_position)
            continue
        query_positions.append(query_index[query_key])
        reference_positions.append(reference_index[reference_key])
    coverage = len(query_positions) / max(len(declared), 1)
    if missing_coordinate:
        raise InputError(
            f"alignment map declares exact residues absent from coordinates for {reference_id}: "
            + ", ".join(map(str, missing_coordinate[:12]))
        )
    if coverage < minimum_coverage:
        raise InputError(
            f"mapped alignment coverage for {reference_id} is {coverage:.3f}, below {minimum_coverage:.3f}"
        )
    if len(query_positions) < 3:
        raise InputError(f"mapped alignment for {reference_id} has fewer than three residues")
    included = [declared[index] for index in range(len(declared)) if declared[index] in shared]
    query_states = query_states or {}; reference_states = reference_states or {}
    ambiguous = [
        position for position in declared
        if query_states.get(position) == "ambiguous" or reference_states.get(position) == "ambiguous"
    ]
    declared_missing = [
        position for position in declared
        if query_states.get(position) in {"missing", "missing_coordinate"}
        or reference_states.get(position) in {"missing", "missing_coordinate"}
    ]
    excluded = [
        position for position in declared
        if position not in shared and position not in ambiguous and position not in declared_missing
    ]
    return (
        np.asarray(query_positions, dtype=int), np.asarray(reference_positions, dtype=int),
        {
            "coordinate_system": "uniprot", "uniprot_start": start, "uniprot_end": end,
            "positions_declared": len(declared), "positions_mapped": len(query_positions),
            "coverage": round(coverage, 6), "included_uniprot_positions": included,
            "excluded_uniprot_positions": excluded,
            "missing_coordinate_positions": declared_missing,
            "ambiguous_uniprot_positions": ambiguous,
        },
    )


def _superpose(moving: np.ndarray, fixed: np.ndarray) -> tuple[np.ndarray, float]:
    moving_center = moving.mean(axis=0); fixed_center = fixed.mean(axis=0)
    x, y = moving - moving_center, fixed - fixed_center
    u, _s, vt = np.linalg.svd(x.T @ y)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1; rotation = u @ vt
    placed = x @ rotation + fixed_center
    rmsd = float(np.sqrt(np.mean(np.sum((placed - fixed) ** 2, axis=1))))
    return placed, rmsd


def _aligned_frames(frames: list[np.ndarray], indices: np.ndarray | None = None) -> list[np.ndarray]:
    """Align complete frames using only the declared comparison residues."""
    reference = frames[0]
    selected = np.arange(len(reference), dtype=int) if indices is None else np.asarray(indices, dtype=int)
    fixed = reference[selected]
    aligned: list[np.ndarray] = []
    for frame in frames:
        moving = frame[selected]
        moving_center = moving.mean(axis=0); fixed_center = fixed.mean(axis=0)
        x, y = moving - moving_center, fixed - fixed_center
        u, _singular, vt = np.linalg.svd(x.T @ y)
        rotation = u @ vt
        if np.linalg.det(rotation) < 0:
            u[:, -1] *= -1; rotation = u @ vt
        aligned.append((frame - moving_center) @ rotation + fixed_center)
    return aligned


def _clusters(frames: list[np.ndarray], cutoff: float) -> list[int]:
    if cutoff <= 0: raise InputError("RMSD clustering cutoff must be positive")
    if len(frames) == 1: return [1]
    distances = []
    for i, left in enumerate(frames):
        for right in frames[i + 1:]:
            distances.append(float(np.sqrt(np.mean(np.sum((left - right) ** 2, axis=1)))))
    tree = linkage(np.asarray(distances), method="average")
    labels = fcluster(tree, t=cutoff, criterion="distance")
    unique = {old: new for new, old in enumerate(sorted(set(map(int, labels))), start=1)}
    return [unique[int(x)] for x in labels]


def _collective_variables(frames: list[np.ndarray], keys: list[tuple[str, int, str]], config: Iterable[Mapping[str, Any]]) -> dict[str, list[float]]:
    index = {key: i for i, key in enumerate(keys)}; out: dict[str, list[float]] = {}
    for cv in config:
        if str(cv.get("type")) != "distance": raise InputError("version 1 collective variables support type=distance")
        def key(name):
            part = cv.get(name, {}); return str(part.get("chain_id", "_")), int(part["auth_seq_id"]), str(part.get("insertion_code", ""))
        a, b = key("a"), key("b")
        if a not in index or b not in index: raise InputError(f"collective variable {cv.get('name')} references an unmapped residue")
        out[str(cv.get("name", "distance"))] = [round(float(np.linalg.norm(frame[index[a]] - frame[index[b]])), 6) for frame in frames]
    return out


def analyze(
    manifest: Mapping[str, Any], reference_set: Mapping[str, Any], *, reference_base: str | Path,
    structure_path: str | Path | None = None, topology_path: str | Path | None = None,
    trajectory_path: str | Path | None = None, chain: str | None = None,
    selection: str = "protein and name CA", stride: int = 1, pbc: str | None = None,
    cluster_cutoff_A: float = 2.0, collective_variables: Iterable[Mapping[str, Any]] = (),
    alignment_map: Mapping[str, Any] | None = None, alignment_map_sha256: str | None = None,
) -> dict[str, Any]:
    errors = validate_reference_set(
        reference_set, base=reference_base, alignment_map=alignment_map,
        alignment_map_sha256=alignment_map_sha256,
    )
    if errors: raise InputError("invalid reference set: " + "; ".join(errors))
    if trajectory_path:
        if not topology_path or pbc is None: raise InputError("trajectory analysis requires topology and explicit --pbc")
        if sha256(topology_path) != str(manifest.get("coordinate", {}).get("sha256", "")):
            raise InputError("StructQC manifest checksum does not match trajectory topology")
        frames, keys, sequence = _trajectory(Path(topology_path), Path(trajectory_path), selection=selection, stride=stride, pbc=pbc)
        input_digest = {"topology": sha256(topology_path), "trajectory": sha256(trajectory_path)}
        input_kind = "trajectory"
    else:
        if not structure_path: raise InputError("structure is required when no trajectory is supplied")
        if selection != DEFAULT_SELECTION:
            raise InputError("--selection applies only to trajectories; static structures require --alignment-map")
        if sha256(structure_path) != str(manifest.get("coordinate", {}).get("sha256", "")):
            raise InputError("StructQC manifest checksum does not match structure")
        frames, keys, sequence = _bio_models(Path(structure_path), chain=chain)
        input_digest = {"structure": sha256(structure_path)}; input_kind = "static_ensemble"

    input_digest.update({
        "structure_evidence_manifest": _json_sha256(manifest),
        "reference_set": _json_sha256(reference_set),
    })
    if alignment_map is not None:
        input_digest["alignment_map"] = alignment_map_sha256 or _json_sha256(alignment_map)

    base = Path(reference_base); references = []
    is_v2 = str(reference_set.get("schema_version", "1.0")) == REFERENCE_SET_V2
    map_details: dict[str, Any] | None = None
    if is_v2:
        _map_errors, map_details = validate_alignment_map(alignment_map or {}, reference_set)
    for raw in reference_set["references"]:
        ref_path = base / str(raw["structure"])
        ref_frames, ref_keys, ref_sequence = _bio_models(ref_path, chain=raw.get("chain"))
        if is_v2 and map_details is not None:
            q_idx, r_idx, mapping = _mapped_indices(
                keys, ref_keys, map_details["query"], map_details["references"][str(raw["reference_id"])],
                start=map_details["start"], end=map_details["end"],
                minimum_coverage=map_details["minimum_coverage"],
                reference_id=str(raw["reference_id"]),
                query_states=map_details["query_states"],
                reference_states=map_details["reference_states"][str(raw["reference_id"])],
            )
        else:
            q_idx, r_idx = _sequence_indices(sequence, ref_sequence)
            mapping = {
                "coordinate_system": "coordinate_sequence_alignment",
                "positions_declared": len(q_idx), "positions_mapped": len(q_idx),
                "coverage": 1.0, "included_uniprot_positions": [],
                "excluded_uniprot_positions": [], "missing_coordinate_positions": [],
                "ambiguous_uniprot_positions": [],
            }
        references.append({"id": str(raw["reference_id"]), "state": str(raw["state"]),
                           "coords": ref_frames[0][r_idx], "query_indices": q_idx,
                           "sha256": sha256(ref_path), "mapping": mapping})
    references.sort(key=lambda item: (item["state"], item["id"]))
    input_digest.update({f"reference_{ref['id']}": ref["sha256"] for ref in references})
    rules = reference_set["decision_rules"]; max_rmsd = float(rules["max_rmsd_A"]); min_margin = float(rules["min_margin_A"])
    frame_rows = []
    for frame_index, frame in enumerate(frames):
        scores = {}
        for ref in references:
            _placed, score = _superpose(frame[ref["query_indices"]], ref["coords"])
            scores[ref["id"]] = score
        by_state = {state: min((scores[r["id"]], r["id"]) for r in references if r["state"] == state) for state in STATES}
        best_state = min(STATES, key=lambda state: by_state[state][0]); other_state = (STATES - {best_state}).pop()
        best_rmsd, best_ref = by_state[best_state]; margin = by_state[other_state][0] - best_rmsd
        call = f"{best_state}_like" if best_rmsd <= max_rmsd and margin >= min_margin else "unresolved"
        frame_rows.append({"frame": frame_index, "call": call, "best_reference": best_ref,
                           "best_rmsd_A": round(best_rmsd, 6), "margin_A": round(margin, 6),
                           **{f"rmsd_{key}_A": round(value, 6) for key, value in sorted(scores.items())}})
    if is_v2 and map_details is not None:
        query_index = _coordinate_index(keys, label="query")
        ensemble_alignment_indices = np.asarray([
            query_index[map_details["query"][position]]
            for position in range(map_details["start"], map_details["end"] + 1)
            if position in map_details["query"] and map_details["query"][position] in query_index
        ], dtype=int)
        declared_count = map_details["end"] - map_details["start"] + 1
        if len(ensemble_alignment_indices) / declared_count < map_details["minimum_coverage"]:
            raise InputError("query coordinates do not meet the declared ensemble-alignment coverage")
    else:
        ensemble_alignment_indices = np.arange(len(keys), dtype=int)
    aligned = _aligned_frames(frames, ensemble_alignment_indices)
    cluster_ids = _clusters(
        [frame[ensemble_alignment_indices] for frame in aligned], cluster_cutoff_A,
    )
    for row, cluster in zip(frame_rows, cluster_ids): row["cluster"] = cluster
    cvs = _collective_variables(frames, keys, collective_variables)
    for name, values in cvs.items():
        for row, value in zip(frame_rows, values): row[f"cv_{name}"] = value
    stack = np.stack(aligned); mean = stack.mean(axis=0)
    rmsf = np.sqrt(np.mean(np.sum((stack - mean) ** 2, axis=2), axis=0))
    counts = {label: sum(row["call"] == label for row in frame_rows) for label in sorted(LABELS)}
    interpreted = counts["active_like"] + counts["inactive_like"]
    if not interpreted: overall = "unresolved"
    elif counts["active_like"] and counts["inactive_like"]: overall = "mixed"
    elif counts["active_like"]: overall = "active_like"
    else: overall = "inactive_like"
    cluster_rows = []
    for cluster in sorted(set(cluster_ids)):
        rows = [r for r in frame_rows if r["cluster"] == cluster]
        calls = {label: sum(r["call"] == label for r in rows) for label in sorted(LABELS)}
        majority = sorted(calls, key=lambda label: (-calls[label], label))[0]
        cluster_rows.append({"cluster": cluster, "frames": len(rows), "fraction_total": round(len(rows) / len(frame_rows), 6),
                             "majority_call": majority, "calls": calls})
    return {
        "schema_version": SCHEMA_VERSION, "module_id": "conformational_state", "subject": manifest["subject"],
        "coordinate_sha256": manifest["coordinate"]["sha256"], "input_kind": input_kind, "input_sha256": input_digest,
        "reference_set": {"id": reference_set.get("reference_set_id", ""),
                          "schema_version": reference_set.get("schema_version", "1.0"),
                          "qualification_scope": reference_set.get("qualification_scope", "unscoped_prototype"),
                          "subject_family": reference_set.get("subject_family", "unresolved"),
                          "references": [{k: r[k] for k in ("id", "state", "sha256", "mapping")} for r in references],
                          "decision_rules": {"max_rmsd_A": max_rmsd, "min_margin_A": min_margin}},
        "config": {"chain": chain, "selection": selection, "stride": stride, "pbc": pbc,
                   "cluster_cutoff_A": cluster_cutoff_A,
                   "ensemble_alignment_residue_count": len(ensemble_alignment_indices),
                   "alignment_map_sha256": input_digest.get("alignment_map")},
        "scientific_scope": {
            "scope_id": reference_set.get("qualification_scope", "unscoped_prototype"),
            "scientific_state": "prototype",
            "supported_subject_class": "ABL-family conformational resemblance" if reference_set.get("qualification_scope") == "abl_family" else "user-declared reference family",
            "release_blocking": bool(reference_set.get("qualification_scope") == "abl_family"),
            "qualification_gate": "qualification_v2_held_out_panel_pending" if is_v2 and reference_set.get("qualification_scope") == "abl_family" else "outside_mark_1_qualified_scope",
            "claim_boundary": "Conformational resemblance is not biochemical activity.",
        },
        "overall_label": overall, "frames_total": len(frame_rows), "frames_interpretable": interpreted,
        "populations": {label: {"count": counts[label], "fraction_total": round(counts[label] / len(frame_rows), 6),
                                "fraction_interpretable": round(counts[label] / interpreted, 6) if interpreted and label != "unresolved" else None}
                        for label in sorted(LABELS)},
        "frame_metrics": frame_rows, "clusters": cluster_rows,
        "rmsf": [{"chain_id": key[0], "auth_seq_id": key[1], "insertion_code": key[2], "rmsf_A": round(float(value), 6)}
                 for key, value in zip(keys, rmsf)],
        "limitations": [
            "State labels describe structural resemblance, not biochemical activity.",
            "Frame populations depend on trajectory preparation, selection, stride, references, and thresholds.",
            "Unresolved frames remain in the total-frame denominator.",
            "Reference Set v1 sequence alignment is retained for compatibility but is not a Mark 1 qualified scope." if not is_v2 else "The exact-mapped ABL-family path remains prototype until the Qualification v2 held-out gate passes.",
        ],
    }


def _json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_outputs(out_dir: str | Path, document: Mapping[str, Any]) -> list[Path]:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    summary, frames, clusters, layer, run = (out / n for n in (
        "STATE_ENSEMBLE.json", "FRAME_METRICS.tsv", "CLUSTERS.tsv", "STATE_LAYER.json", "RUN_MANIFEST.json"))
    _json(summary, document)
    frame_fields = sorted({key for row in document["frame_metrics"] for key in row}, key=lambda k: (k not in {"frame", "call", "cluster", "best_reference", "best_rmsd_A", "margin_A"}, k))
    with frames.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=frame_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(document["frame_metrics"])
    with clusters.open("w", encoding="utf-8", newline="") as handle:
        fields = ["cluster", "frames", "fraction_total", "majority_call"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader(); writer.writerows(document["clusters"])
    _json(layer, {"schema_version": SCHEMA_VERSION, "contract_id": "structure_layer_bundle",
                  "subject": document["subject"], "coordinate_sha256": document["coordinate_sha256"],
                  "layer_id": "conformational_rmsf", "records": [
                      {"chain_id": r["chain_id"], "auth_seq_id": r["auth_seq_id"], "insertion_code": r["insertion_code"],
                       "metric": "rmsf_A", "state": "measured", "value": r["rmsf_A"],
                       "detail": "RMSF after alignment to the first evaluated frame",
                       "source_digest": document["input_sha256"].get("trajectory", document["coordinate_sha256"]),
                       "evidence_class": document["input_kind"]}
                      for r in document["rmsf"]]})
    mda_version = None
    if importlib.util.find_spec("MDAnalysis"):
        try:
            mda_version = importlib.metadata.version("MDAnalysis")
        except importlib.metadata.PackageNotFoundError:
            mda_version = "available_unknown_version"
    _json(run, {"schema_version": SCHEMA_VERSION, "module_id": "conformational_state", "version": "0.2.0",
                "input_sha256": document["input_sha256"],
                "parameters": document["config"],
                "runtime_versions": {"python": platform.python_version(), "biopython": Bio.__version__,
                                     "numpy": np.__version__, "scipy": scipy.__version__},
                "optional_runtimes": {"mdanalysis": (mda_version if document["input_kind"] == "trajectory"
                                                     else "available_not_invoked" if mda_version else "not_available")},
                "outputs": [summary.name, frames.name, clusters.name, layer.name],
                "missing_evidence": [
                    *(["interpretable_state_frames"] if not document["frames_interpretable"] else []),
                    *(["reference_set_v2_exact_alignment_map"] if document["reference_set"].get("schema_version") != REFERENCE_SET_V2 else []),
                    *(["qualification_v2_held_out_panel"]
                      if document.get("scientific_scope", {}).get("qualification_gate") == "qualification_v2_held_out_panel_pending" else []),
                ],
                "limitations": [document.get("scientific_scope", {}).get("claim_boundary", "State resemblance is not activity.")]})
    return [summary, frames, clusters, layer, run]
