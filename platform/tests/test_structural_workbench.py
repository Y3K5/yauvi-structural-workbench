from __future__ import annotations

import hashlib
import json
from pathlib import Path
import zipfile

import pytest

from yauvi_platform.structural_workbench import (
    AnalysisError, StructuralAnalysisStore, analysis_definitions, metric_definitions,
    structural_source_descriptors, template_artifact, tool_readiness,
)
from yauvi_platform.structural_workbench.store import _tree_sha


PDB = """\
ATOM      1  N   ALA A   1      -1.200   0.000   0.000  1.00 91.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 91.00           C
ATOM      3  C   ALA A   1       1.300   0.000   0.000  1.00 91.00           C
ATOM      4  O   ALA A   1       2.100   0.300   0.000  1.00 91.00           O
ATOM      5  N   GLY A   2       1.500   1.200   0.000  1.00 82.00           N
ATOM      6  CA  GLY A   2       2.700   1.500   0.000  1.00 82.00           C
ATOM      7  C   GLY A   2       3.600   0.400   0.000  1.00 82.00           C
ATOM      8  O   GLY A   2       4.800   0.600   0.000  1.00 82.00           O
TER
END
"""


def inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    structure = tmp_path / "model.pdb"
    provenance = tmp_path / "provenance.json"
    validation = tmp_path / "validation.json"
    structure.write_text(PDB, encoding="utf-8")
    provenance.write_text('{"class":"predicted","method":"AlphaFold"}\n', encoding="utf-8")
    validation.write_text(json.dumps({"clashscore": 3.2, "rama_outlier_percent": 0.4}) + "\n", encoding="utf-8")
    return structure, provenance, validation


def test_definitions_freeze_six_core_tasks_and_display_metrics():
    definitions = analysis_definitions()
    assert [item["analysis_type"] for item in definitions] == [
        "structure_qc", "membrane_orientation", "conformational_state",
        "functional_site_state", "assembly_interface", "sf_csa",
    ]
    assert metric_definitions()["rmsd_A"]["decimals"] == 3
    assert all(item["claim_ceiling"] for item in definitions)
    for definition in definitions:
        roles = [item["role"] for item in definition["inputs"]]
        assert len(roles) == len(set(roles))
        assert "validation_report" in roles
        for input_role in definition["inputs"]:
            assert {
                "description", "why_needed", "absence_effect", "accepted_artifact_types",
                "accepted_extensions", "format_guide", "source_ids", "template_id",
                "validator_id", "sensitivity",
            }.issubset(input_role)
        assert definition["use_when"] and definition["measures"] and definition["receives"]
        assert definition["non_claim"]
    artifact_types = {
        artifact["artifact_type"]
        for source in structural_source_descriptors()
        for artifact in source["artifacts"]
    }
    assert {"pdb.coordinates", "pdb.biological_assembly", "wwpdb.validation", "alphafold.pae", "uniprot.sequence"}.issubset(artifact_types)
    assert definitions[-1]["module_ids"] == ["structure_quality", "sf_csa"]
    by_type = {row["analysis_type"]: row for row in definitions}
    assert "topology_evidence" in {row["role"] for row in by_type["membrane_orientation"]["inputs"]}
    assert "topology_evidence" not in {row["role"] for row in by_type["structure_qc"]["inputs"]}
    assert "alignment_map" in {row["role"] for row in by_type["conformational_state"]["inputs"]}
    sources = {source["source_id"]: source for source in structural_source_descriptors()}
    artifacts = {
        artifact["artifact_type"]: artifact
        for source in sources.values() for artifact in source["artifacts"]
    }
    for definition in definitions:
        for role in definition["inputs"]:
            assert set(role["source_ids"]).issubset(sources)
            for artifact_type in role["accepted_artifact_types"]:
                artifact = artifacts[artifact_type]
                if artifact["accepted_extensions"]:
                    assert set(role["accepted_extensions"]) & set(artifact["accepted_extensions"])


def test_source_tree_digest_excludes_local_build_artifacts(tmp_path):
    package = tmp_path / "package"
    package.mkdir()
    (package / "module.py").write_text("VALUE = 1\n", encoding="utf-8")
    expected = _tree_sha(package)
    (package / "build").mkdir()
    (package / "build" / "module.py").write_text("generated\n", encoding="utf-8")
    (package / "package.egg-info").mkdir()
    (package / "package.egg-info" / "PKG-INFO").write_text("generated\n", encoding="utf-8")
    assert _tree_sha(package) == expected


def test_scope_specific_readiness_does_not_promote_experimental_methods(tmp_path):
    rows = {row["analysis_type"]: row for row in tool_readiness(tmp_path)}
    membrane = {row["scope_id"]: row for row in rows["membrane_orientation"]["scientific_scopes"]}
    assert membrane["beta_barrel"]["scientific_state"] == "conditionally_qualified"
    assert membrane["beta_barrel"]["release_blocking"] is True
    assert membrane["alpha_helical"]["scientific_state"] == "prototype"
    assert membrane["alpha_helical"]["release_blocking"] is False
    abl = {row["scope_id"]: row for row in rows["conformational_state"]["scientific_scopes"]}
    assert abl["abl_family"]["scientific_state"] == "prototype"
    assert abl["other_proteins"]["release_blocking"] is False


def test_new_evidence_templates_are_downloadable_and_explicit():
    topology_name, _mime, topology_bytes = template_artifact("membrane_topology_evidence")
    alignment_name, _mime, alignment_bytes = template_artifact("state_alignment_map_v2")
    assert topology_name.endswith(".json") and alignment_name.endswith(".json")
    topology = json.loads(topology_bytes)
    alignment = json.loads(alignment_bytes)
    assert topology["coordinate_sha256"].startswith("replace-")
    assert alignment["domain"] == {"uniprot_start": 242, "uniprot_end": 495}
    assert set(alignment["reference_metadata"]) == {
        "ACTIVE_1", "ACTIVE_2", "INACTIVE_1", "INACTIVE_2",
    }

def test_supporting_sources_do_not_masquerade_as_directly_accepted_files():
    rows = structural_source_descriptors(
        artifact_types=["pdb.coordinates"], source_ids=["pdb", "sifts"],
    )
    by_id = {row["source_id"]: row for row in rows}
    assert [item["artifact_type"] for item in by_id["pdb"]["artifacts"]] == ["pdb.coordinates"]
    assert by_id["sifts"]["artifacts"] == []
    assert structural_source_descriptors(artifact_types=[], source_ids=[]) == []


def test_bounded_ingest_rejects_bad_hash_and_duplicate_chunks(tmp_path):
    store = StructuralAnalysisStore(tmp_path)
    store.create("qc-one", analysis_type="structure_qc", question="Is it inspectable?")
    upload = store.begin_ingest(
        "qc-one", role="structure", file_name="model.pdb", size=3,
        expected_sha256="0" * 64,
    )
    store.ingest_chunk(upload["upload_id"], 0, b"ABC")
    with pytest.raises(AnalysisError, match="match the declared checksum"):
        store.finalize_ingest(upload["upload_id"])
    with pytest.raises(AnalysisError, match="once and in order"):
        store.ingest_chunk(upload["upload_id"], 0, b"ABC")


def test_server_finalized_checksum_supports_streaming_browsers_without_web_crypto(tmp_path):
    store = StructuralAnalysisStore(tmp_path)
    store.create("qc-stream", analysis_type="structure_qc", question="Is it inspectable?")
    payload = b"ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n"
    upload = store.begin_ingest(
        "qc-stream", role="structure", file_name="model.pdb", size=len(payload),
        expected_sha256="",
    )
    assert upload["checksum_policy"] == "server_finalized"
    store.ingest_chunk(upload["upload_id"], 0, payload)
    finalized = store.finalize_ingest(upload["upload_id"])
    assert finalized["sha256"] == hashlib.sha256(payload).hexdigest()


def test_structqc_case_runs_cli_and_writes_deterministic_report_bundle(tmp_path):
    structure, provenance, validation = inputs(tmp_path)
    store = StructuralAnalysisStore(tmp_path)
    store.create(
        "qc-case", analysis_type="structure_qc",
        question="Are the supplied coordinates suitable for interpretation?",
        subject_id="SYNTHETIC_1",
    )
    store.add_file("qc-case", role="structure", path=structure)
    store.add_file("qc-case", role="provenance", path=provenance)
    store.add_file("qc-case", role="validation_report", path=validation)
    preflight = store.preflight("qc-case")
    assert preflight["valid"], checks

    first = store.run("qc-case")
    assert first["status"] == "completed"
    run_dir = tmp_path / "structural_analyses" / "cases" / "qc-case" / "runs" / first["run_id"]
    required = {"REPORT_DATA.json", "REPORT.html", "RAW_EVIDENCE.zip", "CHECKSUMS.json", "RUN_MANIFEST.json"}
    assert required.issubset({path.name for path in run_dir.iterdir()})
    before = {name: (run_dir / name).read_bytes() for name in required}
    second = store.run("qc-case")
    after = {name: (run_dir / name).read_bytes() for name in required}
    assert second["run_id"] == first["run_id"]
    assert before == after
    assert str(tmp_path).encode() not in (run_dir / "RUN_MANIFEST.json").read_bytes()
    report = json.loads((run_dir / "REPORT_DATA.json").read_text(encoding="utf-8"))
    assert report["platform_identity"] == {
        "display_name": "YAUVI Structural Biology Platform — Mark 1",
        "edition": "Mark 1",
        "platform_id": "yauvi_structural_biology_platform_mark_1",
        "scientific_suite_name": "YAUVI Structural Workbench",
    }
    assert "YAUVI Structural Biology Platform — Mark 1" in (run_dir / "REPORT.html").read_text(encoding="utf-8")
    validation_doc = next(item["document"] for item in report["documents"] if item["path"].endswith("STRUCTURE_EVIDENCE.json"))
    assert validation_doc["external_validation"]["metrics"]["clashscore"] == 3.2
    with zipfile.ZipFile(run_dir / "RAW_EVIDENCE.zip") as archive:
        assert archive.namelist() == sorted(archive.namelist())
        assert all(info.date_time == (1980, 1, 1, 0, 0, 0) for info in archive.infolist())


def test_input_artifact_access_is_bound_to_the_analysis(tmp_path):
    structure, _provenance, _validation = inputs(tmp_path)
    store = StructuralAnalysisStore(tmp_path)
    store.create("bound-one", analysis_type="structure_qc", question="q")
    item = store.add_file("bound-one", role="structure", path=structure)
    assert store.input_path("bound-one", item["sha256"]).read_bytes() == structure.read_bytes()
    store.create("bound-two", analysis_type="structure_qc", question="q")
    with pytest.raises(AnalysisError, match="not attached"):
        store.input_path("bound-two", hashlib.sha256(structure.read_bytes()).hexdigest())


def test_sf_csa_preflight_requires_checksum_pinned_pack_and_organism_neutral_tables(tmp_path):
    structure, _provenance, _validation = inputs(tmp_path)
    fasta = tmp_path / "query.faa"; fasta.write_text(">SYN1\nAG\n", encoding="utf-8")
    proteome = tmp_path / "proteome.faa"; proteome.write_text(">SYN1\nAG\n", encoding="utf-8")
    tables = tmp_path / "interpretation.json"
    tables.write_text(json.dumps({
        "mechanism_families": [{"group": "synthetic_group", "pattern": "synthetic"}],
        "contested_groups": [], "divergence_sets": [],
        "classification_vocabulary": ["structural_analogy_only", "unresolved_or_conflicted"],
    }), encoding="utf-8")
    pack = tmp_path / ".yauvi-cache" / "sources" / "packs" / "synthetic-pack"
    pack.mkdir(parents=True)
    database = pack / "database_manifest.json"
    database.write_text('{"schema_version":1,"path_base":"."}\n', encoding="utf-8")
    (pack / "PACK_MANIFEST.json").write_text(json.dumps({
        "sf_csa_database_manifest": database.name,
        "database_manifest_sha256": hashlib.sha256(database.read_bytes()).hexdigest(),
    }), encoding="utf-8")
    store = StructuralAnalysisStore(tmp_path)
    store.create("sf-case", analysis_type="sf_csa", question="Which relationships are bounded?", subject_id="SYN1")
    for role, path in (("query_structure", structure), ("query_fasta", fasta),
                       ("source_proteome", proteome), ("interpretation_tables", tables)):
        store.add_file("sf-case", role=role, path=path)
    store.update_parameters("sf-case", {
        "accession": "SYN1", "organism": "Synthetic organism", "mechanism_group": "synthetic_group",
        "protein_specific_boundary": "A synthetic fold does not establish function.",
        "database_pack": "synthetic-pack", "chain": "A",
    })
    checks = {item["name"]: item for item in store.preflight("sf-case")["checks"]}
    assert checks["reference:sf_csa_database_pack"]["ok"]
    assert checks["evidence:sf_csa_interpretation_tables"]["ok"]


def test_alpha_helical_preflight_requires_coordinate_bound_topology(tmp_path):
    structure = tmp_path / "helix.pdb"
    structure.write_text("\n".join([
        f"ATOM  {index:5d}  CA  ALA A{index:4d}    {float(index):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
        for index in range(1, 7)
    ] + ["END", ""]), encoding="utf-8")
    store = StructuralAnalysisStore(tmp_path)
    store.create("alpha-case", analysis_type="membrane_orientation", question="Where is the membrane axis?")
    item = store.add_file("alpha-case", role="structure", path=structure)
    store.update_parameters("alpha-case", {"context": "eukaryotic_pm"})
    checks = {row["name"]: row for row in store.preflight("alpha-case")["checks"]}
    assert checks["evidence:membrane_topology"]["ok"] is False

    bad_topology = tmp_path / "bad-topology.json"
    bad_topology.write_text(json.dumps({
        "coordinate_sha256": item["sha256"],
        "source": {"id": "synthetic-topology", "citation": "synthetic test evidence"},
        "spans": [{"chain_id": "A", "start_auth_seq_id": 2, "end_auth_seq_id": 7}],
    }), encoding="utf-8")
    bad_store = StructuralAnalysisStore(tmp_path / "bad-workspace")
    bad_store.create("alpha-bad", analysis_type="membrane_orientation", question="Where is the membrane axis?")
    bad_store.add_file("alpha-bad", role="structure", path=structure)
    bad_store.add_file("alpha-bad", role="topology_evidence", path=bad_topology)
    bad_store.update_parameters("alpha-bad", {"context": "eukaryotic_pm"})
    bad_checks = {row["name"]: row for row in bad_store.preflight("alpha-bad")["checks"]}
    assert bad_checks["evidence:membrane_topology"]["ok"] is False
    assert "absent from the coordinate model" in bad_checks["evidence:membrane_topology"]["detail"]

    topology = tmp_path / "topology.json"
    topology.write_text(json.dumps({
        "coordinate_sha256": item["sha256"],
        "source": {"id": "synthetic-topology", "citation": "synthetic test evidence"},
        "spans": [{"chain_id": "A", "start_auth_seq_id": 1, "end_auth_seq_id": 6}],
    }), encoding="utf-8")
    store.add_file("alpha-case", role="topology_evidence", path=topology)
    checks = {row["name"]: row for row in store.preflight("alpha-case")["checks"]}
    assert checks["evidence:membrane_topology"]["ok"] is True, checks["evidence:membrane_topology"]


def test_abl_preflight_requires_complete_exact_v2_map(tmp_path):
    def long_pdb(path: Path) -> None:
        lines = []
        for index in range(254):
            residue = index + 1
            lines.append(
                f"ATOM  {residue:5d}  CA  ALA A{residue:4d}    {float(index):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
            )
        path.write_text("\n".join([*lines, "END", ""]), encoding="utf-8")

    query = tmp_path / "query.pdb"
    active = tmp_path / "active.pdb"; active_two = tmp_path / "active-two.pdb"
    inactive = tmp_path / "inactive.pdb"; inactive_two = tmp_path / "inactive-two.pdb"
    for path in (query, active, active_two, inactive, inactive_two):
        long_pdb(path)
    rows = [
        {"uniprot_position": 242 + index, "chain_id": "A", "auth_seq_id": index + 1,
         "insertion_code": "", "mapping_state": "exact"}
        for index in range(254)
    ]
    alignment = tmp_path / "alignment.json"
    alignment.write_text(json.dumps({
        "schema_version": "2.0", "coordinate_system": "uniprot",
        "source": {"id": "SIFTS:synthetic", "citation": "synthetic exact map", "sha256": "1" * 64},
        "domain": {"uniprot_start": 242, "uniprot_end": 495},
        "query": rows,
        "reference_metadata": {
            "ACTIVE_1": {"pdb_entry_id": "SYN1", "chain_id": "A"},
            "ACTIVE_2": {"pdb_entry_id": "SYN2", "chain_id": "A"},
            "INACTIVE_1": {"pdb_entry_id": "SYN3", "chain_id": "A"},
            "INACTIVE_2": {"pdb_entry_id": "SYN4", "chain_id": "A"},
        },
        "references": {
            "ACTIVE_1": rows, "ACTIVE_2": rows,
            "INACTIVE_1": rows, "INACTIVE_2": rows,
        },
    }, sort_keys=True), encoding="utf-8")
    store = StructuralAnalysisStore(tmp_path)
    store.create("abl-case", analysis_type="conformational_state", question="Which ABL conformation is resembled?")
    for role, path in (("structure", query), ("active_reference", active),
                       ("active_reference", active_two), ("inactive_reference", inactive),
                       ("inactive_reference", inactive_two), ("alignment_map", alignment)):
        store.add_file("abl-case", role=role, path=path)
    store.update_parameters("abl-case", {
        "active_state_evidence": "independent synthetic active-state assignment",
        "inactive_state_evidence": "independent synthetic inactive-state assignment",
        "active_reference_citation": "synthetic active citation",
        "inactive_reference_citation": "synthetic inactive citation",
        "reference_method": "synthetic coordinates",
        "subject_family": "ABL1",
    })
    preflight = store.preflight("abl-case")
    checks = {row["name"]: row for row in preflight["checks"]}
    assert preflight["valid"], checks
    assert checks["evidence:state_alignment_map_v2"]["ok"]
