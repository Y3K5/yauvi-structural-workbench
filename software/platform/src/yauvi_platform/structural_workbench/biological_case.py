"""Identity-bound biological cases; read-only views over explicit local evidence.

This module joins records and applies deposited assembly operators. Scientific
measurements belong to the existing engines, not to the viewer or this adapter.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import re
from pathlib import Path
from urllib.parse import urlsplit

import gemmi
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,95}\Z")
STATES = {"supported", "contradicted", "disputed", "unknown", "not_applicable"}
BASES = {"experimental", "curated", "predicted", "derived"}
MAX_ARTIFACT = 32 * 1024 * 1024
MAX_ATOMS = 300_000


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def identifier(value):
    if not isinstance(value, str) or not IDENTIFIER.fullmatch(value):
        raise ValueError("Invalid biological case identifier")
    return value


def _required_text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Missing {label}")
    return value


def _rows(doc, prefix, fields):
    columns = [doc.get(prefix + field, []) for field in fields]
    if len({len(c) for c in columns}) > 1:
        raise ValueError(f"Incomplete mmCIF category: {prefix}")
    return [dict(zip(fields, values)) for values in zip(*columns)]


def _blank(value):
    return "" if value in (None, ".", "?") else str(value)


def _integer(value, label, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError(f"Invalid {label}")
    return value


class BiologicalCase:
    """Open one self-contained local case, verifying every declared artifact."""

    def __init__(self, directory):
        self.directory = Path(directory).resolve()
        manifest = self.directory / "case.json"
        if manifest.is_symlink() or manifest.stat().st_size > 2 * 1024 * 1024:
            raise ValueError("Invalid biological case manifest")
        self.document = json.loads(manifest.read_text())
        d = self.document
        if d.get("schema_version") != "1.0":
            raise ValueError("Unsupported biological case schema")
        identifier(d.get("id"))
        _required_text(d.get("title"), "case title")
        self.sources = {}
        self.payloads = {}
        for source in d.get("sources", []):
            sid = identifier(source.get("id"))
            if sid in self.sources:
                raise ValueError("Duplicate source identity")
            relative = Path(source["path"])
            path = (self.directory / relative).resolve()
            if relative.is_absolute() or ".." in relative.parts or not path.is_relative_to(self.directory):
                raise ValueError("Case source escapes its directory")
            if path.stat().st_size > MAX_ARTIFACT:
                raise ValueError("Case source is too large")
            data = path.read_bytes()
            if sha256(data) != source.get("sha256"):
                raise ValueError(f"Source checksum mismatch: {sid}")
            _required_text(source.get("label"), "source label")
            _required_text(source.get("retrieved_at"), "source acquisition date")
            url = source.get("url")
            if url:
                parts = urlsplit(url)
                if parts.scheme != "https" or not parts.hostname or parts.username or parts.password:
                    raise ValueError("Source citation must be an HTTPS URL")
            self.sources[sid], self.payloads[sid] = source, data
        self.protein = d["protein"]
        sequence = self.protein["sequence"]
        if not isinstance(sequence, str) or not re.fullmatch("[ACDEFGHIKLMNPQRSTVWYXBZUO]+", sequence):
            raise ValueError("Invalid reference protein sequence")
        if sha256(sequence.encode("ascii")) != self.protein.get("sequence_sha256"):
            raise ValueError("Reference sequence checksum mismatch")
        _required_text(self.protein.get("accession"), "protein accession")
        _required_text(self.protein.get("organism"), "protein organism")
        protein_source = json.loads(self._source(self.protein["source_id"]))
        if (protein_source.get("primaryAccession") != self.protein["accession"]
                or protein_source.get("sequence", {}).get("value") != sequence
                or protein_source.get("organism", {}).get("scientificName") != self.protein["organism"]):
            raise ValueError("Protein identity differs from its locked source record")
        self.structures = {}
        for structure in d.get("structures", []):
            sid = identifier(structure["id"])
            if sid in self.structures:
                raise ValueError("Duplicate structure identity")
            data = self._source(structure["source_id"])
            if not self.sources[structure["source_id"]]["path"].lower().endswith((".cif", ".mmcif")):
                raise ValueError("Biological cases require deposited mmCIF coordinates")
            doc = MMCIF2Dict(io.StringIO(data.decode("utf-8")))
            if doc.get("_entry.id", [""])[0].upper() != structure["entry_id"].upper():
                raise ValueError("Coordinate entry identity mismatch")
            parsed = gemmi.make_structure_from_block(gemmi.cif.read_string(data.decode()).sole_block())
            if set(doc.get("_pdbx_struct_assembly.id", [])) != {a.name for a in parsed.assemblies}:
                raise ValueError("Assembly declarations could not be completely parsed")
            mapping = self._sequence_mapping(doc, structure.get("sequence_mappings", []))
            models = sorted(set(doc.get("_atom_site.pdbx_PDB_model_num", [])), key=int)
            if not models:
                raise ValueError("No coordinate models")
            self.structures[sid] = (structure, doc, parsed, mapping, models)
        if not self.structures:
            raise ValueError("No structures declared")
        seen = set()
        for assertion in d.get("assertions", []):
            aid = identifier(assertion["id"])
            if aid in seen:
                raise ValueError("Duplicate assertion identity")
            seen.add(aid)
            if assertion.get("state") not in STATES or assertion.get("basis") not in BASES:
                raise ValueError("Invalid assertion evidence state or basis")
            _required_text(assertion.get("property"), "assertion property")
            _required_text(assertion.get("interpretation_limit"), "assertion interpretation limit")
            self._source(assertion["source_id"])
            scope = assertion.get("scope", {})
            if scope.get("sequence_sha256") != self.protein["sequence_sha256"]:
                raise ValueError("Assertion has a different sequence identity")
            if "start" in scope or "end" in scope:
                lo = _integer(scope.get("start"), "assertion start")
                hi = _integer(scope.get("end"), "assertion end")
                if lo > hi or hi > len(sequence):
                    raise ValueError("Assertion residue range outside reference")
            if "structure_id" in scope:
                self._binding({key: scope.get(key) for key in (
                    "structure_id", "coordinate_sha256", "sequence_sha256", "model_id", "assembly_id")})
        self.evidence = [self._adapt(item) for item in d.get("evidence", [])]
        if len({item["id"] for item in self.evidence}) != len(self.evidence):
            raise ValueError("Duplicate engine evidence identity")
        for membrane in d.get("membranes", []):
            self._binding(membrane["binding"])
            self._membrane(membrane["binding"])

    def _source(self, source_id):
        if source_id not in self.payloads:
            raise ValueError("Evidence refers to an undeclared source")
        return self.payloads[source_id]

    def _sequence_mapping(self, doc, records):
        entities = {r["entity_id"]: "".join(r["pdbx_seq_one_letter_code_can"].split())
                    for r in _rows(doc, "_entity_poly.", ["entity_id", "pdbx_seq_one_letter_code_can"])}
        mapping = {}
        used = set()
        for row in records:
            entity = str(row["entity_id"])
            if entity not in entities:
                raise ValueError("Mapping names an unknown polymer entity")
            start = _integer(row["label_start"], "polymer start")
            reference = _integer(row["reference_start"], "reference start")
            length = _integer(row["length"], "mapping length")
            observed = entities[entity][start - 1:start - 1 + length]
            expected = self.protein["sequence"][reference - 1:reference - 1 + length]
            if len(observed) != length or observed != expected:
                raise ValueError("Declared polymer mapping does not match the exact reference sequence")
            for offset in range(length):
                key, target = (entity, start + offset), (entity, reference + offset)
                if key in mapping or target in used:
                    raise ValueError("Overlapping or ambiguous polymer mapping")
                mapping[key] = reference + offset
                used.add(target)
        return mapping

    def _binding(self, binding):
        sid = binding.get("structure_id")
        if sid not in self.structures:
            raise ValueError("Evidence names an unknown structure")
        structure, _, parsed, _, models = self.structures[sid]
        if binding.get("coordinate_sha256") != self.sources[structure["source_id"]]["sha256"]:
            raise ValueError("Evidence has a different coordinate identity")
        if binding.get("sequence_sha256") != self.protein["sequence_sha256"]:
            raise ValueError("Evidence has a different sequence identity")
        if binding.get("model_id") not in models:
            raise ValueError("Evidence names an unknown model")
        if binding.get("assembly_id") not in {"asu", *(a.name for a in parsed.assemblies)}:
            raise ValueError("Evidence names an unknown assembly")
        return structure

    def _adapt(self, item):
        """Verify a run/output binding before exposing raw engine measurements.

        Older engines omit complete input identity. A separately locked run
        record is therefore mandatory; it binds input scope and exact output.
        The adapter never upgrades an engine measurement to native exposure.
        """
        identifier(item["id"])
        binding = item["binding"]
        self._binding(binding)
        raw_bytes = self._source(item["source_id"])
        raw = json.loads(raw_bytes)
        run = json.loads(self._source(item["run_source_id"]))
        if run.get("binding") != binding or run.get("output_sha256") != sha256(raw_bytes):
            raise ValueError("Engine run/output binding mismatch")
        if not run.get("method") or not re.fullmatch(r"[0-9a-f]{64}", str(run.get("method_source_sha256", ""))):
            raise ValueError("Engine method provenance missing")
        engine = item["engine"]
        if engine not in {"structqc", "assembly_context", "memorient", "site_context", "actstate", "state_atlas", "sf_csa"}:
            raise ValueError("Unsupported biological evidence adapter")
        internal_hash = raw.get("input_sha256", {}).get("structure")
        if internal_hash and internal_hash != binding["coordinate_sha256"]:
            raise ValueError("Engine's internal coordinate identity conflicts with case")
        if engine == "structqc":
            coordinate = raw.get("coordinate", {})
            if coordinate.get("sha256") != binding["coordinate_sha256"]:
                raise ValueError("StructQC coordinate identity mismatch")
            models = self.structures[binding["structure_id"]][4]
            if coordinate.get("selected_model") != models.index(binding["model_id"]):
                raise ValueError("StructQC model identity mismatch")
            if binding["assembly_id"] != "asu":
                raise ValueError("Deposited-coordinate StructQC cannot be relabeled as expanded assembly evidence")
            summary = raw.get("completeness", {})
        elif engine == "assembly_context":
            if raw.get("assembly_sha256") != binding["coordinate_sha256"]:
                raise ValueError("AssemblyContext coordinate identity mismatch")
            if raw.get("reference", {}).get("assembly_id") != binding["assembly_id"]:
                raise ValueError("AssemblyContext assembly identity mismatch")
            summary = raw.get("summary", raw.get("assembly", {}))
        else:
            summary = raw.get("summary", {})
        return {"id": item["id"], "engine": engine, "source_id": item["source_id"],
                "binding": binding, "state": "recorded_method_output", "summary": summary,
                "native_accessibility": "unknown",
                "interpretation_limit": "A recorded calculation is not evidence of native accessibility, biological activity, or qualification."}

    def summary(self):
        structures = []
        for sid, (record, doc, parsed, _, models) in self.structures.items():
            descriptions = dict(zip(doc.get("_pdbx_struct_assembly.id", []), doc.get("_pdbx_struct_assembly.details", [])))
            structures.append({"id": sid, "entry_id": record["entry_id"],
                "sha256": self.sources[record["source_id"]]["sha256"], "models": models,
                "assemblies": [{"id": "asu", "label": "Deposited asymmetric unit"}] + [
                    {"id": a.name, "label": f"Assembly {a.name} · {descriptions.get(a.name, 'deposited recipe')}"}
                    for a in parsed.assemblies]})
        return {"id": self.document["id"], "title": self.document["title"], "protein": self.protein,
                "structures": structures, "assertions": self.document.get("assertions", []),
                "sources": [{k: v for k, v in source.items() if k != "path"} for source in self.sources.values()],
                "limitations": self.document.get("limitations", [])}

    def view(self, structure_id, model_id, assembly_id):
        if structure_id not in self.structures:
            raise ValueError("Unknown structure selection")
        structure, doc, parsed, mapping, models = self.structures[structure_id]
        if model_id not in models:
            raise ValueError("Select a deposited model explicitly")
        fields = ["id", "group_PDB", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
                  "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z",
                  "occupancy", "auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num"]
        atoms_in = [r for r in _rows(doc, "_atom_site.", fields) if r["pdbx_PDB_model_num"] == model_id]
        asym_ids = sorted({r["label_asym_id"] for r in atoms_in})
        copies = []
        if assembly_id == "asu":
            copies = [(label, "identity", np.eye(3), np.zeros(3)) for label in asym_ids]
        else:
            assembly = next((a for a in parsed.assemblies if a.name == assembly_id), None)
            if assembly is None:
                raise ValueError("Select a deposited assembly explicitly")
            for gen in assembly.generators:
                if not gen.subchains:
                    raise ValueError("Assembly lacks label-asym identities")
                for operator in gen.operators:
                    for label in gen.subchains:
                        if label not in asym_ids:
                            raise ValueError("Assembly operator references a missing chain")
                        copies.append((label, operator.name, np.array(operator.transform.mat.tolist()),
                                       np.array(operator.transform.vec.tolist())))
        if len(copies) > 1000 or len(copies) * len(atoms_in) > MAX_ATOMS * max(len(asym_ids), 1):
            raise ValueError("Selected assembly exceeds viewer capacity")
        scheme = _rows(doc, "_pdbx_poly_seq_scheme.", ["asym_id", "entity_id", "seq_id", "mon_id", "auth_seq_num", "pdb_seq_num", "pdb_ins_code", "pdb_strand_id"])
        atoms, chains, warnings = [], [], []
        seen_copies = set()
        for label, operator, rotation, translation in copies:
            chain_id = f"{label}:{operator}"
            if chain_id in seen_copies:
                raise ValueError("Duplicate generated-chain identity")
            seen_copies.add(chain_id)
            source_atoms = [a for a in atoms_in if a["label_asym_id"] == label]
            author = source_atoms[0]["auth_asym_id"]
            residues = {}
            for row in scheme:
                if row["asym_id"] != label:
                    continue
                label_seq = int(row["seq_id"])
                position = mapping.get((row["entity_id"], label_seq))
                rid = f"{chain_id}/{label_seq}"
                residues[str(label_seq)] = {"id": rid, "label_seq_id": label_seq,
                    "sequence_position": position, "auth_seq_id": _blank(row["auth_seq_num"]) or _blank(row["pdb_seq_num"]),
                    "insertion_code": _blank(row["pdb_ins_code"]), "comp_id": row["mon_id"],
                    "one_letter": self.protein["sequence"][position - 1] if position else "X",
                    "observed": False, "mapping_state": "exact_polymer_sequence" if position else "unmapped",
                    "atom_indices": [], "alternate_locations": []}
            # Select one residue conformer by mean occupancy; retain blank atoms.
            alt_scores = {}
            for a in source_atoms:
                residue_key = a["label_seq_id"] if _blank(a["label_seq_id"]) else f"het:{a['auth_seq_id']}:{a['pdbx_PDB_ins_code']}:{a['label_comp_id']}"
                alt = _blank(a["label_alt_id"])
                if alt:
                    alt_scores.setdefault(residue_key, {}).setdefault(alt, []).append(float(a["occupancy"]))
            chosen = {key: sorted(scores, key=lambda alt: (-sum(scores[alt])/len(scores[alt]), alt))[0]
                      for key, scores in alt_scores.items()}
            for a in source_atoms:
                key = a["label_seq_id"] if _blank(a["label_seq_id"]) else f"het:{a['auth_seq_id']}:{a['pdbx_PDB_ins_code']}:{a['label_comp_id']}"
                alt = _blank(a["label_alt_id"])
                if alt and chosen[key] != alt:
                    continue
                if key not in residues:
                    if _blank(a["label_seq_id"]):
                        raise ValueError("Coordinate polymer residue absent from deposited sequence scheme")
                    residues[key] = {"id": f"{chain_id}/{key}", "sequence_position": None, "label_seq_id": None,
                        "auth_seq_id": a["auth_seq_id"], "insertion_code": _blank(a["pdbx_PDB_ins_code"]),
                        "comp_id": a["label_comp_id"], "one_letter": "X", "observed": False,
                        "mapping_state": "nonpolymer", "atom_indices": [], "alternate_locations": []}
                residue = residues[key]
                if residue["comp_id"] != a["label_comp_id"]:
                    raise ValueError("Coordinate and polymer scheme residue identities disagree")
                residue.update({"auth_seq_id": a["auth_seq_id"], "insertion_code": _blank(a["pdbx_PDB_ins_code"]), "observed": True})
                residue["alternate_locations"] = sorted(alt_scores.get(key, {}))
                residue["selected_altloc"] = chosen.get(key, "")
                xyz = rotation @ np.array([float(a["Cartn_x"]), float(a["Cartn_y"]), float(a["Cartn_z"])]) + translation
                if not np.isfinite(xyz).all():
                    raise ValueError("Nonfinite atom coordinates")
                index = len(atoms)
                atoms.append({"serial": index, "index": index, "x": float(xyz[0]), "y": float(xyz[1]), "z": float(xyz[2]),
                    "chain": chain_id, "resi": int(a["auth_seq_id"]), "icode": residue["insertion_code"],
                    "resn": a["label_comp_id"], "atom": a["label_atom_id"], "elem": a["type_symbol"],
                    "hetflag": a["group_PDB"] == "HETATM", "properties": {"residue_id": residue["id"]}})
                residue["atom_indices"].append(index)
            chains.append({"id": chain_id, "label": f"{author} · label {label} · operator {operator}",
                "source_label_asym_id": label, "auth_asym_id": author, "operator_ids": [operator],
                "rotation": rotation.tolist(), "translation": translation.tolist(), "residues": list(residues.values())})
        if len(atoms) > MAX_ATOMS:
            raise ValueError("Selected assembly exceeds viewer capacity")
        if any(r["alternate_locations"] for c in chains for r in c["residues"]):
            warnings.append("One conformer per residue is displayed, selected by mean occupancy; alternate-location IDs remain recorded.")
        binding = {"structure_id": structure_id, "coordinate_sha256": self.sources[structure["source_id"]]["sha256"],
                   "sequence_sha256": self.protein["sequence_sha256"], "model_id": model_id, "assembly_id": assembly_id}
        evidence = [e for e in self.evidence if e["binding"] == binding]
        membrane = self._membrane(binding)
        return {"case_id": self.document["id"], **binding, "chains": chains, "atoms": atoms,
                "membrane": membrane, "evidence": evidence, "warnings": warnings}

    def _membrane(self, binding):
        # A membrane overlay needs its own input-frame record. A legacy normal
        # alone cannot reconstruct the center or an assembly-wide bilayer.
        matches = [m for m in self.document.get("membranes", []) if m.get("binding") == binding]
        if len(matches) > 1:
            raise ValueError("Conflicting membrane placements require separate case alternatives")
        if not matches:
            return None
        m = matches[0]
        self._binding(m["binding"])
        raw = json.loads(self._source(m["source_id"]))
        if raw.get("binding") != binding or raw.get("frame") != "deposited_coordinates":
            raise ValueError("Membrane frame/identity mismatch")
        if any(a["property"] == "membrane_applicability" and a["state"] == "not_applicable" for a in self.document.get("assertions", [])):
            raise ValueError("A soluble context cannot carry a membrane overlay")
        normal, center = np.asarray(raw["normal"], dtype=float), np.asarray(raw["center"], dtype=float)
        half = float(raw["half_thickness"])
        if normal.shape != (3,) or center.shape != (3,) or not np.isfinite(normal).all() or not np.isfinite(center).all() or not math.isfinite(half) or half <= 0:
            raise ValueError("Invalid membrane geometry")
        if not np.isclose(np.linalg.norm(normal), 1, atol=1e-6):
            raise ValueError("Membrane normal is not a unit vector")
        if raw.get("sidedness") != "unknown":
            raise ValueError("Signed membrane overlays require a supported sidedness adapter")
        return {"state": "research_only", "normal": normal.tolist(), "center": center.tolist(),
                "half_thickness": half, "sidedness": "unknown", "source_id": m["source_id"]}


class BiologicalCaseStore:
    def __init__(self, workspace):
        self.root = Path(workspace).resolve() / "biological-cases"

    def load(self, case_id):
        directory = self.root / identifier(case_id)
        if directory.is_symlink() or not directory.resolve().is_relative_to(self.root.resolve()):
            raise ValueError("Case directory escapes workspace")
        case = BiologicalCase(directory)
        if case.document["id"] != case_id:
            raise ValueError("Case directory identity mismatch")
        return case

    def list(self):
        if not self.root.exists():
            return []
        result = []
        for directory in sorted(self.root.iterdir()):
            if directory.is_dir() and not directory.is_symlink() and (directory / "case.json").is_file():
                try:
                    case = self.load(directory.name)
                    result.append({"id": case.document["id"], "title": case.document["title"], "state": "ready"})
                except (OSError, ValueError, KeyError, TypeError) as exc:
                    result.append({"id": directory.name, "title": directory.name, "state": "invalid", "error": str(exc)})
        return result
