#!/usr/bin/env python3
"""Build the nominated public integration cases from explicit local public inputs.

No downloads, sequence submission, publication, or overwrite of existing cases.
This is an integration example builder, not a membrane qualification runner.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from yauvi_platform.structural_workbench.biological_case import BiologicalCase


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def build(inputs, workspace):
    results = []
    for case_id, entry, accession, title, is_membrane in [
        ("aqp1", "1J4N", "P47865", "AQP1 · a water channel in biological context", True),
        ("carbonic-anhydrase", "2POW", "P00918", "Carbonic anhydrase II · soluble control", False),
    ]:
        directory = workspace / "biological-cases" / case_id
        if directory.exists():
            raise ValueError(f"Case already exists; use a new output workspace: {case_id}")
        protein = json.loads((inputs / (accession + ".json")).read_text())
        if protein["primaryAccession"] != accession:
            raise ValueError("Input protein accession mismatch")
        sequence = protein["sequence"]["value"]
        seq_sha = digest(sequence.encode())
        doc = MMCIF2Dict(str(inputs / (entry + ".cif")))
        if doc["_entry.id"][0].upper() != entry:
            raise ValueError("Input structure identity mismatch")
        mapping = [{"entity_id": entity, "label_start": 1, "reference_start": 1, "length": len(sequence)}
                   for entity, value in zip(doc["_entity_poly.entity_id"], doc["_entity_poly.pdbx_seq_one_letter_code_can"])
                   if "".join(value.split()) == sequence]
        if not mapping:
            raise ValueError("Nominated example no longer has an exact polymer/reference match")
        directory.mkdir(parents=True)
        sources = []
        now = datetime.now(timezone.utc).isoformat()

        def add_source(sid, filename, label, url=None):
            row = {"id": sid, "path": filename, "label": label, "sha256": digest((directory / filename).read_bytes()),
                   "retrieved_at": now, "acquisition_note": "Imported from explicitly selected local public inputs; original retrieval time is not asserted."}
            if url:
                row["url"] = url
            sources.append(row)
            return row

        shutil.copyfile(inputs / (accession + ".json"), directory / "protein.json")
        shutil.copyfile(inputs / (entry + ".cif"), directory / "structure.cif")
        add_source("protein", "protein.json", f"UniProt {accession}", f"https://www.uniprot.org/uniprotkb/{accession}/entry")
        coordinate = add_source("coordinates", "structure.cif", f"wwPDB {entry} deposited coordinates", f"https://www.rcsb.org/structure/{entry}")
        assertions = [{"id": "membrane-applicability", "property": "membrane_applicability", "value": is_membrane,
            "state": "supported" if is_membrane else "not_applicable", "basis": "curated", "source_id": "protein",
            "scope": {"sequence_sha256": seq_sha},
            "interpretation_limit": "Curated setting for this nominated example; not inferred from coordinate geometry."}]
        for index, feature in enumerate(protein.get("features", [])):
            if feature["type"] not in {"Transmembrane", "Topological domain", "Motif", "Active site", "Binding site"}:
                continue
            loc = feature.get("location", {})
            if any(loc.get(k, {}).get("modifier") != "EXACT" for k in ("start", "end")):
                continue
            assertions.append({"id": f"uniprot-feature-{index}", "property": feature["type"],
                "value": {"description": feature.get("description", ""), "evidences": feature.get("evidences", [])},
                "state": "supported", "basis": "curated", "source_id": "protein",
                "scope": {"sequence_sha256": seq_sha, "start": loc["start"]["value"], "end": loc["end"]["value"]},
                "interpretation_limit": "UniProt annotation, including its evidence codes; not independently measured by this workbench."})
        assertions.append({"id": "native-accessibility", "property": "native_accessibility", "value": None,
            "state": "unknown", "basis": "curated", "source_id": "protein", "scope": {"sequence_sha256": seq_sha},
            "interpretation_limit": "Neither solvent exposure nor a membrane placement establishes probe-specific native accessibility."})
        manifest = {"schema_version": "1.0", "id": case_id, "title": title,
            "protein": {"accession": accession, "sequence": sequence, "sequence_sha256": seq_sha,
                "organism": protein["organism"]["scientificName"], "source_id": "protein", "entry_audit": protein.get("entryAudit", {})},
            "structures": [{"id": entry.lower(), "entry_id": entry, "source_id": "coordinates", "sequence_mappings": mapping}],
            "assertions": assertions, "sources": sources, "evidence": [], "membranes": [],
            "limitations": ["Integration example; not membrane algorithm qualification.",
                "Assembly alternatives are deposited proposals, not proof of native oligomeric state.",
                "Missing coordinates and unmapped annotations remain unresolved.",
                "Native accessibility, transport rate, and biological activity are not established by this view."]}
        # Reuse the existing StructQC engine, preserving its actual scope and gaps.
        qc = importlib.import_module("structqc.core")
        raw = qc.analyze(directory / "structure.cif", subject_id=accession,
                         reference_id=accession, reference_sequence=sequence)
        write(directory / "structqc.json", raw)
        output = add_source("structqc", "structqc.json", "StructQC result for deposited model 1")
        binding = {"structure_id": entry.lower(), "coordinate_sha256": coordinate["sha256"],
                   "sequence_sha256": seq_sha, "model_id": "1", "assembly_id": "asu"}
        write(directory / "structqc-run.json", {"binding": binding, "output_sha256": output["sha256"],
            "method": "structqc.core.analyze", "method_source_sha256": digest(Path(qc.__file__).read_bytes()),
            "parameters": {"model_index": 0, "chain": None, "reference_id": accession}, "scope": "integration example"})
        add_source("structqc-run", "structqc-run.json", "StructQC input/output binding")
        manifest["evidence"].append({"id": "structqc", "engine": "structqc", "binding": binding,
                                     "source_id": "structqc", "run_source_id": "structqc-run"})
        def record_engine(engine, eid, result, engine_binding, module, parameters):
            write(directory / (eid + ".json"), result)
            artifact = add_source(eid, eid + ".json", engine + " recorded calculation")
            package = Path(module.__file__).parent
            code = {p.relative_to(package).as_posix(): digest(p.read_bytes()) for p in sorted(package.rglob("*.py"))}
            write(directory / (eid + "-run.json"), {"binding": engine_binding, "output_sha256": artifact["sha256"],
                "method": module.__name__, "method_source_sha256": digest(json.dumps(code, sort_keys=True).encode()),
                "source_files": code, "parameters": parameters, "scope": "integration example; no qualification claim"})
            add_source(eid + "-run", eid + "-run.json", engine + " input/output/method binding")
            manifest["evidence"].append({"id": eid, "engine": engine, "binding": engine_binding,
                "source_id": eid, "run_source_id": eid + "-run"})

        assembly_module = importlib.import_module("assembly_context.core")
        for assembly_id in doc.get("_pdbx_struct_assembly.id", []):
            assembly_binding = {**binding, "assembly_id": assembly_id}
            result = assembly_module.analyze(raw, directory / "structure.cif", directory / "structure.cif",
                subject_chain="A", relationship="exact_protein", reference_id=entry, assembly_id=assembly_id, sasa_backend="biopython")
            record_engine("assembly_context", "assembly-" + assembly_id, result, assembly_binding, assembly_module,
                          {"subject_chain": "A", "relationship": "exact_protein", "assembly_id": assembly_id, "sasa_backend": "biopython"})

        if is_membrane:
            # Derive explicit residue keys from independent UniProt annotations.
            # The OPM reference orientation is not used to choose these spans.
            sequence_rows = list(zip(doc["_pdbx_poly_seq_scheme.asym_id"], doc["_pdbx_poly_seq_scheme.seq_id"],
                doc["_pdbx_poly_seq_scheme.auth_seq_num"], doc["_pdbx_poly_seq_scheme.pdb_ins_code"],
                doc["_pdbx_poly_seq_scheme.pdb_strand_id"]))
            spans = []
            for feature in assertions:
                if feature["property"] != "Transmembrane":
                    continue
                scope = feature["scope"]
                residues = [{"chain_id": author, "auth_seq_id": int(auth), "insertion_code": "" if ins in {".", "?"} else ins}
                    for label, label_seq, auth, ins, author in sequence_rows
                    if label == "A" and scope["start"] <= int(label_seq) <= scope["end"] and auth not in {"?", "."}]
                if len(residues) != scope["end"] - scope["start"] + 1:
                    raise ValueError("A nominated transmembrane span lacks exact coordinate mapping")
                spans.append({"residues": residues})
            topology = {"coordinate_sha256": coordinate["sha256"], "spans": spans,
                "source": {"id": accession, "citation": sources[0]["url"], "sha256": sources[0]["sha256"]}}
            write(directory / "topology.json", topology)
            add_source("topology", "topology.json", "UniProt transmembrane spans mapped to exact author residues")
            orientation = importlib.import_module("memorient.orientor")
            geometry = importlib.import_module("memorient.geometry")
            contexts = importlib.import_module("memorient.contexts")
            result = orientation.orient_structure(geometry.load_structure(str(directory / "structure.cif"), chain="A"),
                contexts.get_context("tm_receptor"), topology_evidence=topology, validate=False, n_points=240)
            raw_membrane = result.to_dict()
            record_engine("memorient", "membrane", raw_membrane, binding, orientation,
                          {"context": "tm_receptor", "chain": "A", "validate": False, "n_points": 240,
                           "topology_sha256": digest((directory / "topology.json").read_bytes())})
            pose = raw_membrane["input_coordinate_membrane"]
            write(directory / "membrane-pose.json", {**pose, "frame": "deposited_coordinates", "binding": binding,
                "derived_from_sha256": digest((directory / "membrane.json").read_bytes())})
            add_source("membrane-pose", "membrane-pose.json", "Unsigned modeled membrane in the deposited coordinate frame")
            manifest["membranes"].append({"binding": binding, "source_id": "membrane-pose"})
        write(directory / "case.json", manifest)
        case = BiologicalCase(directory)
        for item in case.summary()["structures"]:
            for assembly in item["assemblies"]:
                view = case.view(item["id"], "1", assembly["id"])
                results.append({"case": case_id, "assembly": assembly["id"], "atoms": len(view["atoms"]),
                    "polymer_copies": sum(any(r["sequence_position"] for r in c["residues"]) for c in view["chains"]),
                    "source_sha256": coordinate["sha256"]})
    write(workspace / "BIOLOGICAL_CASE_BUILD.json", {"scope": "integration only", "views": results})
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--workspace", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.inputs.resolve(), args.workspace.resolve()), indent=2))
