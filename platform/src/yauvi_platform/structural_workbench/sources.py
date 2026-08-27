"""Registered public artifacts for the local structural workbench.

This module is intentionally an identity and policy layer.  It never accepts an
arbitrary URL and it never derives a query from an uploaded structure.  A user
must choose a declared artifact type and enter a public accession.  Retrieval is
performed by :mod:`yauvi_sources`, while this module records the result in a
content-addressed cache that can later be adopted into an analysis revision.
"""
from __future__ import annotations

import hashlib
import json
import re
import secrets
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping


SCHEMA_VERSION = "1.0"
SAFE_ACQUISITION_ID = re.compile(r"^source_[0-9a-f]{24}$")


class StructuralSourceError(RuntimeError):
    """A source request is unknown, unsafe, unavailable, or unverifiable."""


def _canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(_canonical(value))


def _descriptor(
    source_id: str,
    provider: str,
    homepage: str,
    documentation_url: str,
    license_name: str,
    artifacts: list[dict[str, Any]],
    *,
    limitations: Iterable[str] = (),
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "contract_id": "structural_source_descriptor",
        "source_id": source_id,
        "provider": provider,
        "homepage": homepage,
        "documentation_url": documentation_url,
        "license": license_name,
        "artifacts": artifacts,
        "limitations": list(limitations),
    }


def _artifact(
    artifact_type: str,
    label: str,
    description: str,
    extensions: Iterable[str],
    identifier_kind: str,
    identifier_example: str,
    *,
    fetchable: bool,
    generated_locally: bool = False,
) -> dict[str, Any]:
    return {
        "contract_id": "source_artifact_type",
        "artifact_type": artifact_type,
        "label": label,
        "description": description,
        "accepted_extensions": list(extensions),
        "identifier_kind": identifier_kind,
        "identifier_example": identifier_example,
        "fetchable": fetchable,
        "generated_locally": generated_locally,
    }


SOURCE_DESCRIPTORS: tuple[dict[str, Any], ...] = (
    _descriptor(
        "pdb", "RCSB Protein Data Bank", "https://www.rcsb.org/",
        "https://www.rcsb.org/docs/programmatic-access/file-download-services", "CC0",
        [
            _artifact("pdb.coordinates", "Experimental coordinates", "Deposited asymmetric-unit coordinates in PDBx/mmCIF.", [".cif", ".mmcif"], "PDB ID", "1CRN", fetchable=True),
            _artifact("pdb.biological_assembly", "Biological assembly", "A deposited biological assembly selected by PDB ID and assembly number.", [".cif", ".mmcif"], "PDB ID:assembly", "4HHB:1", fetchable=True),
        ],
        limitations=("A deposited assembly is an author/provider model and is not proof of the dominant assembly in the tested biological context.",),
    ),
    _descriptor(
        "wwpdb_validation", "wwPDB validation archive", "https://www.wwpdb.org/validation/validation-reports",
        "https://www.rcsb.org/docs/general-help/validation-reports", "CC0",
        [_artifact("wwpdb.validation", "wwPDB validation XML", "Community validation metrics for an experimental PDB entry.", [".xml"], "PDB ID", "1CRN", fetchable=True)],
        limitations=("Validation metrics assess the deposited model and experiment; they do not establish biological function.",),
    ),
    _descriptor(
        "alphafold_db", "AlphaFold Protein Structure Database", "https://alphafold.ebi.ac.uk/",
        "https://alphafold.ebi.ac.uk/api-docs", "CC BY 4.0",
        [
            _artifact("alphafold.model", "AlphaFold coordinate model", "Predicted single-chain coordinates resolved from the current AlphaFold API record.", [".pdb", ".cif", ".mmcif"], "UniProt accession", "P69905", fetchable=True),
            _artifact("alphafold.pae", "Predicted aligned error", "Pairwise confidence values for an AlphaFold prediction.", [".json"], "UniProt accession", "P69905", fetchable=True),
        ],
        limitations=("Predicted coordinates and confidence are not experimental structure or biochemical evidence.",),
    ),
    _descriptor(
        "uniprot_proteomes", "UniProtKB", "https://www.uniprot.org/",
        "https://www.uniprot.org/help/api_queries", "CC BY 4.0",
        [
            _artifact("uniprot.sequence", "Protein sequence FASTA", "The current UniProtKB sequence for one public accession.", [".fasta", ".fa", ".faa"], "UniProt accession", "P69905", fetchable=True),
            _artifact("uniprot.annotations", "UniProt feature table", "A TSV containing catalytic, binding, cofactor, function, and cross-reference fields.", [".tsv"], "UniProt accession", "P69905", fetchable=True),
            _artifact("uniprot.proteome", "Reference proteome FASTA", "A UniProt reference proteome used for a declared sequence-comparison universe.", [".fasta", ".fa", ".faa"], "Proteome ID", "UP000005640", fetchable=True),
        ],
        limitations=("Annotations have heterogeneous evidence levels; an annotation is not an observation in the submitted structure.",),
    ),
    _descriptor(
        "sifts", "SIFTS", "https://www.ebi.ac.uk/pdbe/docs/sifts/",
        "https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html", "CC BY 4.0",
        [_artifact("sifts.mapping", "Residue mapping", "PDB-to-UniProt and related residue mappings, normalized before use.", [".json", ".xml", ".tsv"], "PDB ID", "1CRN", fetchable=False)],
        limitations=("Ambiguous or incomplete mappings cannot drive residue-level conclusions.",),
    ),
    _descriptor(
        "mcsa", "Mechanism and Catalytic Site Atlas", "https://www.ebi.ac.uk/thornton-srv/m-csa/",
        "https://www.ebi.ac.uk/thornton-srv/m-csa/download/", "CC BY 4.0",
        [_artifact("mcsa.sites", "Catalytic-site records", "Curated catalytic residues, roles, cofactors, and mechanism records.", [".json", ".csv", ".tsv"], "M-CSA entry or PDB ID", "1B73", fetchable=False)],
        limitations=("A curated catalytic annotation does not prove activity in the submitted coordinate state.",),
    ),
    _descriptor(
        "pdb_ccd", "PDB Chemical Component Dictionary", "https://www.wwpdb.org/data/ccd",
        "https://www.wwpdb.org/data/ccd", "CC0",
        [_artifact("pdb_ccd.component", "Chemical component definition", "Declared connectivity and identifiers for a PDB chemical component.", [".cif", ".mmcif"], "CCD component ID", "ATP", fetchable=False)],
    ),
    _descriptor(
        "chebi", "ChEBI", "https://www.ebi.ac.uk/chebi/",
        "https://www.ebi.ac.uk/chebi/downloads", "CC BY 4.0",
        [_artifact("chebi.mapping", "ChEBI identifier mapping", "Chemical identifiers and exact synonyms used to normalize declared cofactors.", [".json", ".tsv", ".obo"], "ChEBI ID", "CHEBI:15422", fetchable=False)],
    ),
    _descriptor(
        "opm_ppm", "OPM / PPM", "https://opm.phar.umich.edu/",
        "https://opm.phar.umich.edu/ppm_server3", "Provider terms apply",
        [_artifact("opm_ppm.reference", "Membrane-positioning reference", "A curator-frozen comparison record for membrane placement benchmarks.", [".pdb", ".json"], "PDB ID", "1BXW", fetchable=False)],
        limitations=("OPM/PPM is an external comparison standard, not proof of native intact-cell exposure.",),
    ),
    _descriptor(
        "foldseek", "Foldseek", "https://github.com/steineggerlab/foldseek",
        "https://github.com/steineggerlab/foldseek/wiki", "GPL-3.0",
        [_artifact("foldseek.database_pack", "Frozen Foldseek database pack", "A checksum-pinned database snapshot built outside the browser and adopted as a registered pack.", [".json"], "Registered pack ID", "sf-csa-public-mini-v1", fetchable=False, generated_locally=True)],
        limitations=("Fold similarity is not exact protein identity or functional transfer.",),
    ),
    _descriptor(
        "freesasa", "FreeSASA", "https://freesasa.github.io/",
        "https://freesasa.github.io/doxygen/CLI.html", "MIT",
        [_artifact("freesasa.runtime", "FreeSASA runtime", "Local publication-grade SASA executable resolved through readiness checks.", [], "Runtime", "freesasa --version", fetchable=False, generated_locally=True)],
    ),
)


TEMPLATES: Mapping[str, tuple[str, str, bytes]] = {
    "explicit_structure_provenance": (
        "structure-provenance.template.json", "application/json",
        _canonical({
            "provenance_class": "unknown",
            "source": {"namespace": "local", "identifier": "replace-me", "release": "unknown"},
            "prediction_format": None,
            "confidence_encoding": "none_declared",
            "limitations": ["Replace every placeholder; unknown provenance remains unknown."],
        }),
    ),
    "functional_site_annotations": (
        "functional-sites.template.json", "application/json",
        _canonical({
            "sites": [{"site_id": "site-1", "chain_id": "A", "auth_seq_id": 1, "role": "unspecified", "source_id": "replace-me", "evidence_class": "annotated"}],
            "limitations": ["Roles and residue coordinates must come from an explicit source."],
        }),
    ),
    "membrane_topology_evidence": (
        "membrane-topology-evidence.template.json", "application/json",
        _canonical({
            "schema_version": "1.0",
            "coordinate_sha256": "replace-with-the-exact-coordinate-sha256",
            "source": {
                "id": "replace-with-source-record-id",
                "citation": "replace-with-source-citation",
            },
            "spans": [
                {"chain_id": "A", "start_auth_seq_id": 10, "end_auth_seq_id": 30},
            ],
            "sidedness": {
                "extracellular_residue": {
                    "chain_id": "A", "auth_seq_id": 31, "insertion_code": "",
                },
            },
            "limitations": [
                "Replace all placeholders. Omit sidedness when no external topology evidence supports it.",
            ],
        }),
    ),
    "state_alignment_map_v2": (
        "state-alignment-map-v2.template.json", "application/json",
        _canonical({
            "schema_version": "2.0",
            "coordinate_system": "uniprot",
            "source": {
                "id": "SIFTS-or-explicit-equivalence-record",
                "citation": "replace-with-mapping-source-citation",
                "sha256": "replace-with-source-record-sha256",
            },
            "domain": {"uniprot_start": 242, "uniprot_end": 495},
            "query": [
                {"uniprot_position": 242, "chain_id": "A", "auth_seq_id": 242, "insertion_code": "", "mapping_state": "exact"},
            ],
            "reference_metadata": {
                "ACTIVE_1": {"pdb_entry_id": "replace-active-1", "chain_id": "A"},
                "ACTIVE_2": {"pdb_entry_id": "replace-active-2", "chain_id": "A"},
                "INACTIVE_1": {"pdb_entry_id": "replace-inactive-1", "chain_id": "A"},
                "INACTIVE_2": {"pdb_entry_id": "replace-inactive-2", "chain_id": "A"},
            },
            "references": {
                "ACTIVE_1": [
                    {"uniprot_position": 242, "chain_id": "A", "auth_seq_id": 242, "insertion_code": "", "mapping_state": "exact"},
                ],
                "ACTIVE_2": [
                    {"uniprot_position": 242, "chain_id": "A", "auth_seq_id": 242, "insertion_code": "", "mapping_state": "exact"},
                ],
                "INACTIVE_1": [
                    {"uniprot_position": 242, "chain_id": "A", "auth_seq_id": 242, "insertion_code": "", "mapping_state": "exact"},
                ],
                "INACTIVE_2": [
                    {"uniprot_position": 242, "chain_id": "A", "auth_seq_id": 242, "insertion_code": "", "mapping_state": "exact"},
                ],
            },
            "limitations": [
                "This is a format template, not a valid one-residue map. Analysis requires at least two independently supported references per state and at least 90 percent exact coverage of UniProt ABL1 residues 242-495.",
            ],
        }),
    ),
    "component_identifier_map": (
        "component-map.template.json", "application/json",
        _canonical({"components": [{"observed_id": "ATP", "ccd_id": "ATP", "chebi_id": "CHEBI:15422", "mapping_state": "exact"}]}),
    ),
    "sf_csa_interpretation_tables": (
        "sf-csa-interpretation.template.json", "application/json",
        _canonical({"mechanism_families": [], "contested_groups": [], "divergence_sets": [], "classification_vocabulary": []}),
    ),
}


def structural_source_descriptors(
    *, artifact_types: Iterable[str] | None = None, source_ids: Iterable[str] | None = None,
) -> list[dict[str, Any]]:
    """Return providers and distinguish direct artifacts from supporting links.

    With no filters this returns the complete catalog. A role-specific request
    may additionally name supporting providers. Supporting providers are
    returned with an empty ``artifacts`` list, so the UI offers official links
    without suggesting a file that the role's validator cannot ingest.
    """
    wanted = None if artifact_types is None else set(artifact_types)
    supporting = None if source_ids is None else set(source_ids)
    rows: list[dict[str, Any]] = []
    for source in SOURCE_DESCRIPTORS:
        artifacts = [
            dict(item) for item in source["artifacts"]
            if wanted is None or item["artifact_type"] in wanted
        ]
        if artifacts or (supporting is not None and source["source_id"] in supporting):
            row = dict(source)
            row["artifacts"] = artifacts
            rows.append(row)
    return rows


def template_artifact(template_id: str) -> tuple[str, str, bytes]:
    try:
        return TEMPLATES[template_id]
    except KeyError:
        raise StructuralSourceError(f"unknown structural template: {template_id}") from None


class StructuralSourceStore:
    """Persistent acquisition records over the shared content-addressed cache."""

    def __init__(self, workspace: str | Path, *, fetcher: Callable[[str, str], Any] | None = None) -> None:
        self.workspace = Path(workspace).resolve()
        self.root = self.workspace / ".yauvi-cache" / "sources" / "structural"
        self.jobs_root = self.root / "acquisitions"
        self.jobs_root.mkdir(parents=True, exist_ok=True)
        try:
            from yauvi_sources import SourceCache
        except ImportError as exc:
            raise StructuralSourceError("yauvi-sources is required for verified source acquisition") from exc
        self.cache = SourceCache(self.workspace / ".yauvi-cache" / "sources")
        if fetcher is None:
            try:
                from yauvi_sources.fetchers import fetch_structural_artifact
            except ImportError as exc:
                raise StructuralSourceError("the structural source fetcher is unavailable") from exc
            fetcher = fetch_structural_artifact
        self.fetcher = fetcher

    @staticmethod
    def _artifact_descriptor(artifact_type: str) -> tuple[dict[str, Any], dict[str, Any]]:
        for source in SOURCE_DESCRIPTORS:
            for artifact in source["artifacts"]:
                if artifact["artifact_type"] == artifact_type:
                    return source, artifact
        raise StructuralSourceError(f"unknown structural artifact type: {artifact_type}")

    def acquire(self, artifact_type: str, identifier: str) -> dict[str, Any]:
        source, artifact = self._artifact_descriptor(artifact_type)
        if not artifact["fetchable"]:
            raise StructuralSourceError(f"{artifact_type} is link-only or locally generated and cannot be fetched by the workbench")
        identifier = identifier.strip()
        if not identifier or len(identifier) > 64 or not re.fullmatch(r"[A-Za-z0-9:_-]+", identifier):
            raise StructuralSourceError("public identifier contains unsupported characters")
        acquisition_id = "source_" + secrets.token_hex(12)
        outcome = self.fetcher(artifact_type, identifier)
        if not getattr(outcome, "ok", False):
            record = {
                "schema_version": SCHEMA_VERSION, "contract_id": "source_acquisition_request",
                "acquisition_id": acquisition_id, "artifact_type": artifact_type,
                "identifier": identifier, "state": "failed", "reason": str(getattr(outcome, "reason", "source_failed")),
            }
            _write_json(self.jobs_root / f"{acquisition_id}.json", record)
            raise StructuralSourceError(record["reason"])
        cache_source_id = f"{source['source_id']}__{artifact_type.replace('.', '_')}"
        entry = self.cache.store(
            cache_source_id, outcome.payload, filename=outcome.filename,
            origin=outcome.origin, version=outcome.version,
            note=f"artifact_type={artifact_type}; identifier={identifier}",
        )
        artifact_manifest = {
            "schema_version": SCHEMA_VERSION, "contract_id": "source_artifact_manifest",
            "source_id": source["source_id"], "provider": source["provider"],
            "artifact_type": artifact_type, "identifier": identifier,
            "file_name": entry.filename, "sha256": entry.sha256, "bytes": entry.bytes,
            "origin": entry.origin, "release": entry.version,
            "license": source["license"], "format_validation": "passed",
            "cache_source_id": cache_source_id,
        }
        artifact_manifest["manifest_sha256"] = hashlib.sha256(_canonical(artifact_manifest)).hexdigest()
        record = {
            "schema_version": SCHEMA_VERSION, "contract_id": "source_acquisition_request",
            "acquisition_id": acquisition_id, "artifact_type": artifact_type,
            "identifier": identifier, "state": "completed", "artifact": artifact_manifest,
        }
        _write_json(self.jobs_root / f"{acquisition_id}.json", record)
        return record

    def load(self, acquisition_id: str) -> dict[str, Any]:
        if not SAFE_ACQUISITION_ID.fullmatch(acquisition_id):
            raise StructuralSourceError("invalid acquisition id")
        try:
            return json.loads((self.jobs_root / f"{acquisition_id}.json").read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise StructuralSourceError("source acquisition is unavailable") from exc

    def artifact_path(self, acquisition_id: str) -> Path:
        record = self.load(acquisition_id)
        artifact = record.get("artifact")
        if record.get("state") != "completed" or not isinstance(artifact, Mapping):
            raise StructuralSourceError("source acquisition did not produce an adoptable artifact")
        entries = self.cache.entries(str(artifact["cache_source_id"]))
        entry = next((item for item in entries if item.sha256 == artifact["sha256"] and item.filename == artifact["file_name"]), None)
        if entry is None:
            raise StructuralSourceError("cached artifact manifest is unavailable")
        path = self.cache.path_for(entry)
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != artifact["sha256"]:
            raise StructuralSourceError("cached artifact checksum is mismatched")
        return path
