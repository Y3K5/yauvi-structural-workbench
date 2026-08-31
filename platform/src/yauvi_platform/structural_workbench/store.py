"""Content-addressed inputs, registered structural CLIs, and deterministic reports.

The workbench is deliberately an orchestration layer.  Scientific calculations
remain in the standalone packages; this module validates typed inputs, invokes
only registered commands, records exact bytes, and renders their evidence without
inventing a score or scientific conclusion.
"""
from __future__ import annotations

import csv
import hashlib
import html
import importlib.util
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import threading
import zipfile
from pathlib import Path
from typing import Any, Iterable, Mapping


SCHEMA_VERSION = "1.0"
PLATFORM_ID = "yauvi_structural_biology_platform_mark_1"
PLATFORM_DISPLAY_NAME = "YAUVI Structural Biology Platform — Mark 1"
PLATFORM_SCIENTIFIC_SUITE = "YAUVI Structural Workbench"
SAFE_ID = re.compile(r"^[a-z0-9][a-z0-9_-]{0,63}$")
SAFE_UPLOAD_ID = re.compile(r"^upload_[0-9a-f]{24}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
MAX_CHUNK_BYTES = 8 * 1024 * 1024
MAX_INPUT_BYTES = 20 * 1024 * 1024 * 1024
FIXED_ZIP_TIME = (1980, 1, 1, 0, 0, 0)


class AnalysisError(RuntimeError):
    """A structural analysis cannot safely proceed."""


def _canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(_canonical(value))


def _sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _tree_sha(path: Path) -> str | None:
    if not path.is_dir():
        return None
    digest = hashlib.sha256()
    def scientific_source_file(candidate: Path) -> bool:
        excluded = {"__pycache__", ".pytest_cache", ".hypothesis", "build", "dist"}
        return (
            candidate.is_file()
            and not excluded.intersection(candidate.parts)
            and not any(part.endswith(".egg-info") for part in candidate.parts)
            and candidate.suffix not in {".pyc", ".pyo"}
        )
    for item in sorted(p for p in path.rglob("*") if scientific_source_file(p)):
        digest.update(item.relative_to(path).as_posix().encode())
        digest.update(b"\0")
        digest.update(item.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _checkout_source_root(workspace: Path) -> Path | None:
    if (workspace / "structqc").is_dir() and (workspace / "state-atlas").is_dir():
        return workspace
    candidate = Path(__file__).resolve().parents[4]
    if (candidate / "structqc").is_dir() and (candidate / "state-atlas").is_dir():
        return candidate
    return None


def _installed_package_root(import_name: str) -> Path | None:
    spec = importlib.util.find_spec(import_name)
    if spec is None:
        return None
    if spec.submodule_search_locations:
        return Path(next(iter(spec.submodule_search_locations))).resolve()
    return Path(spec.origin).resolve().parent if spec.origin else None


def metric_definitions() -> dict[str, dict[str, Any]]:
    return {
        "distance_A": {"label": "Distance", "unit": "Å", "decimals": 2, "kind": "distance", "help": "Euclidean coordinate distance."},
        "rmsd_A": {"label": "Root-mean-square deviation", "unit": "Å", "decimals": 3, "kind": "distance", "help": "Coordinate resemblance after the declared alignment."},
        "rmsf_A": {"label": "Root-mean-square fluctuation", "unit": "Å", "decimals": 3, "kind": "distance", "help": "Per-residue positional variation across an interpretable ensemble."},
        "sasa_A2": {"label": "Solvent-accessible surface area", "unit": "Å²", "decimals": 1, "kind": "area", "help": "Method-specific solvent-accessible area."},
        "buried_sasa_A2": {"label": "Buried surface area", "unit": "Å²", "decimals": 1, "kind": "area", "help": "Difference between isolated and assembly-accessible area using the recorded method."},
        "fraction": {"label": "Fraction", "unit": "", "decimals": 3, "kind": "fraction", "help": "A value from 0 to 1; its denominator is reported with the result."},
        "percentage": {"label": "Percentage", "unit": "%", "decimals": 1, "kind": "percentage", "help": "Presentation of a recorded fraction; raw values are retained."},
        "plddt": {"label": "pLDDT", "unit": "", "decimals": 1, "kind": "confidence", "help": "Prediction confidence only when provenance declares a compatible format."},
        "mapping_confidence": {"label": "Mapping confidence", "unit": "", "decimals": 3, "kind": "confidence", "help": "Confidence in the declared sequence-to-coordinate identity mapping."},
        "count": {"label": "Count", "unit": "", "decimals": 0, "kind": "count", "help": "An exact item count."},
    }


def _base_analysis_definitions() -> list[dict[str, Any]]:
    """Public task definitions consumed by the CLI and browser UI."""
    return [
        {
            "analysis_type": "structure_qc", "title": "Structure QC",
            "question": "Are these coordinates identity-bound, provenance-declared, and suitable for interpretation?",
            "module_ids": ["structure_quality"], "readiness": "prototype",
            "claim_ceiling": "Coordinate and provenance quality only; not native structure or function.",
            "inputs": [
                {"role": "structure", "label": "PDB or mmCIF coordinates", "required": True, "multiple": False, "extensions": [".pdb", ".cif", ".mmcif"]},
                {"role": "reference_fasta", "label": "Reference sequence", "required": False, "multiple": False, "extensions": [".fasta", ".fa", ".faa"]},
                {"role": "provenance", "label": "Provenance declaration", "required": False, "multiple": False, "extensions": [".json"]},
                {"role": "pae", "label": "Predicted aligned error", "required": False, "multiple": False, "extensions": [".json"]},
                {"role": "validation_report", "label": "wwPDB or local validation report", "required": False, "multiple": False, "extensions": [".xml", ".json"]},
            ],
            "parameters": [
                {"name": "model", "label": "Model index", "type": "integer", "default": 0, "advanced": True},
                {"name": "chain", "label": "Chain", "type": "text", "default": "", "advanced": True},
            ],
        },
        {
            "analysis_type": "membrane_orientation", "title": "Membrane orientation",
            "question": "How does this protein sit in its declared membrane or surface context?",
            "module_ids": ["structure_quality", "membrane_orientation"], "readiness": "conditionally_qualified",
            "claim_ceiling": "Modeled orientation and accessibility; not native intact-cell exposure.",
            "inputs": [
                {"role": "structure", "label": "PDB or mmCIF coordinates", "required": True, "multiple": False, "extensions": [".pdb", ".cif", ".mmcif"]},
                {"role": "topology_evidence", "label": "Coordinate-bound transmembrane spans", "required": False, "multiple": False, "extensions": [".json"]},
                {"role": "reference_fasta", "label": "Reference sequence", "required": False, "multiple": False, "extensions": [".fasta", ".fa", ".faa"]},
                {"role": "provenance", "label": "Provenance declaration", "required": False, "multiple": False, "extensions": [".json"]},
                {"role": "validation_report", "label": "wwPDB or local validation report", "required": False, "multiple": False, "extensions": [".xml", ".json"]},
            ],
            "parameters": [
                {"name": "context", "label": "Biological membrane context", "type": "select", "default": "gram_negative_om",
                 "choices": ["gram_negative_om", "eukaryotic_pm", "tm_receptor", "gram_positive_surface", "soluble_secreted"]},
                {"name": "chain", "label": "Chain", "type": "text", "default": "", "advanced": True},
            ],
        },
        {
            "analysis_type": "conformational_state", "title": "Conformational resemblance",
            "question": "Which experimentally bounded conformations does this structure or ensemble resemble?",
            "module_ids": ["structure_quality", "conformational_state"], "readiness": "prototype",
            "claim_ceiling": "Active-like or inactive-like structural resemblance; not biochemical activity.",
            "inputs": [
                {"role": "structure", "label": "Query structure or trajectory topology", "required": True, "multiple": False, "extensions": [".pdb", ".cif", ".mmcif"]},
                {"role": "trajectory", "label": "Optional trajectory", "required": False, "multiple": False, "extensions": [".xtc", ".dcd", ".trr"]},
                {"role": "active_reference", "label": "Experimental active-state references", "required": True, "multiple": True, "minimum_files": 2, "extensions": [".pdb", ".cif", ".mmcif"]},
                {"role": "inactive_reference", "label": "Experimental inactive-state references", "required": True, "multiple": True, "minimum_files": 2, "extensions": [".pdb", ".cif", ".mmcif"]},
                {"role": "alignment_map", "label": "Exact SIFTS or residue-equivalence map", "required": True, "multiple": False, "extensions": [".json"]},
                {"role": "reference_fasta", "label": "Reference sequence", "required": False, "multiple": False, "extensions": [".fasta", ".fa", ".faa"]},
                {"role": "provenance", "label": "Query provenance declaration", "required": False, "multiple": False, "extensions": [".json"]},
                {"role": "validation_report", "label": "wwPDB or local validation report", "required": False, "multiple": False, "extensions": [".xml", ".json"]},
            ],
            "parameters": [
                {"name": "active_state_evidence", "label": "Basis for active reference state", "type": "text", "required": True},
                {"name": "inactive_state_evidence", "label": "Basis for inactive reference state", "type": "text", "required": True},
                {"name": "active_reference_citation", "label": "Active reference citation or source ID", "type": "text", "required": True},
                {"name": "inactive_reference_citation", "label": "Inactive reference citation or source ID", "type": "text", "required": True},
                {"name": "reference_method", "label": "Experimental method", "type": "text", "required": True},
                {"name": "subject_family", "label": "Qualified reference family", "type": "select", "default": "ABL1", "choices": ["ABL1"]},
                {"name": "max_rmsd_A", "label": "Maximum interpretable RMSD (Å)", "type": "number", "default": 2.5, "advanced": True},
                {"name": "min_margin_A", "label": "Minimum between-state margin (Å)", "type": "number", "default": 0.25, "advanced": True},
                {"name": "cluster_cutoff_A", "label": "Clustering cutoff (Å)", "type": "number", "default": 2.0, "advanced": True},
                {"name": "chain", "label": "Query chain", "type": "text", "default": "", "advanced": True},
                {"name": "stride", "label": "Trajectory stride", "type": "integer", "default": 1, "advanced": True},
                {"name": "pbc", "label": "Periodic boundary handling", "type": "select", "default": "none", "choices": ["none", "unwrap"], "advanced": True},
            ],
        },
        {
            "analysis_type": "functional_site_state", "title": "Functional-site evidence",
            "question": "Are declared functional residues mapped, chemically plausible, and observed in this coordinate state?",
            "module_ids": ["structure_quality", "site_context", "activity_state"], "readiness": "prototype",
            "claim_ceiling": "Site completeness and catalytic-competence evidence; not observed activity.",
            "inputs": [
                {"role": "structure", "label": "PDB or mmCIF coordinates", "required": True, "multiple": False, "extensions": [".pdb", ".cif", ".mmcif"]},
                {"role": "reference_fasta", "label": "Reference sequence", "required": True, "multiple": False, "extensions": [".fasta", ".fa", ".faa"]},
                {"role": "site_annotations", "label": "Site annotations (JSON, TSV, or CSV)", "required": True, "multiple": False, "extensions": [".json", ".tsv", ".csv"]},
                {"role": "uniprot_annotations", "label": "Optional UniProt feature export for ActState", "required": False, "multiple": False, "extensions": [".tsv", ".csv"]},
                {"role": "component_map", "label": "Exact CCD/ChEBI component map", "required": False, "multiple": False, "extensions": [".json"]},
                {"role": "pocket_result", "label": "Method-declared pocket result", "required": False, "multiple": True, "extensions": [".json"]},
                {"role": "provenance", "label": "Provenance declaration", "required": False, "multiple": False, "extensions": [".json"]},
                {"role": "validation_report", "label": "wwPDB or local validation report", "required": False, "multiple": False, "extensions": [".xml", ".json"]},
            ],
            "parameters": [{"name": "chain", "label": "Chain", "type": "text", "default": "", "advanced": True}],
        },
        {
            "analysis_type": "assembly_interface", "title": "Assembly and interfaces",
            "question": "Which residues contact or become buried in a declared biological assembly?",
            "module_ids": ["structure_quality", "assembly_context"], "readiness": "prototype",
            "claim_ceiling": "One coordinate assembly and interface geometry; not native exposure or binding.",
            "inputs": [
                {"role": "structure", "label": "Isolated subject coordinates", "required": True, "multiple": False, "extensions": [".pdb", ".cif", ".mmcif"]},
                {"role": "assembly", "label": "Expanded biological assembly", "required": True, "multiple": False, "extensions": [".pdb", ".cif", ".mmcif"]},
                {"role": "reference_fasta", "label": "Reference sequence", "required": False, "multiple": False, "extensions": [".fasta", ".fa", ".faa"]},
                {"role": "provenance", "label": "Provenance declaration", "required": False, "multiple": False, "extensions": [".json"]},
                {"role": "validation_report", "label": "wwPDB or local validation report", "required": False, "multiple": False, "extensions": [".xml", ".json"]},
            ],
            "parameters": [
                {"name": "subject_chain", "label": "Subject chain", "type": "text", "required": True},
                {"name": "relationship", "label": "Assembly relationship", "type": "select", "default": "exact_protein", "choices": ["exact_protein", "homolog_assembly", "architecture_analogy", "unresolved"]},
                {"name": "assembly_id", "label": "Assembly identifier", "type": "text", "default": "", "advanced": True},
                {"name": "expected_chains", "label": "Expected chain IDs (comma-separated)", "type": "text", "default": "", "advanced": True},
            ],
        },
        {
            "analysis_type": "sf_csa", "title": "SF-CSA comparison",
            "question": "How do structural similarity and sequence homology compare without collapsing them into one claim?",
            "module_ids": ["structure_quality", "sf_csa"], "readiness": "prototype",
            "claim_ceiling": "Typed similarity, homology, analogy, and divergence evidence; not exact functional transfer.",
            "inputs": [
                {"role": "query_structure", "label": "Query structure", "required": True, "multiple": False, "extensions": [".pdb", ".cif", ".mmcif"]},
                {"role": "query_fasta", "label": "Query sequence FASTA", "required": True, "multiple": False, "extensions": [".fasta", ".fa", ".faa"]},
                {"role": "source_proteome", "label": "Local comparison proteome", "required": True, "multiple": True, "extensions": [".fasta", ".fa", ".faa"]},
                {"role": "interpretation_tables", "label": "Organism-appropriate interpretation tables", "required": True, "multiple": False, "extensions": [".json"]},
                {"role": "provenance", "label": "Query provenance declaration", "required": False, "multiple": False, "extensions": [".json"]},
                {"role": "validation_report", "label": "wwPDB or local validation report", "required": False, "multiple": False, "extensions": [".xml", ".json"]},
            ],
            "parameters": [
                {"name": "accession", "label": "Query accession", "type": "text", "required": True},
                {"name": "organism", "label": "Organism and strain", "type": "text", "required": True},
                {"name": "mechanism_group", "label": "Declared mechanism group", "type": "text", "required": True},
                {"name": "protein_specific_boundary", "label": "Protein-specific interpretation boundary", "type": "text", "required": True},
                {"name": "database_pack", "label": "Pinned Foldseek database pack ID", "type": "text", "required": True},
                {"name": "chain", "label": "Query chain", "type": "text", "default": "A", "advanced": True},
            ],
        },
    ]


_WORKFLOW_GUIDANCE: Mapping[str, Mapping[str, Any]] = {
    "structure_qc": {
        "use_when": "Start here whenever you need to inspect a PDB or mmCIF model before making residue-level claims.",
        "measures": "Coordinate identity, models, chains, numbering, missing backbone atoms, provenance, confidence encoding, reference-sequence mapping, PAE, and imported validation evidence.",
        "receives": ["Structure evidence manifest", "Residue quality table", "Coordinate-bound quality layer", "Reproducibility manifest"],
        "non_claim": "It does not establish the native biological conformation, protein function, or experimental activity.",
    },
    "membrane_orientation": {
        "use_when": "Use for a membrane or surface protein when you need a reproducible coordinate frame and residue-side labels.",
        "measures": "A context-declared membrane placement, orientation transform, membrane depth, and modeled residue accessibility.",
        "receives": ["StructQC evidence", "Orientation record", "Oriented coordinates", "Residue orientation layer"],
        "non_claim": "Modeled orientation is not direct evidence of intact-cell exposure or topology in the tested organism.",
    },
    "conformational_state": {
        "use_when": "Use when experimental references define two bounded conformational states and you want to compare a structure or ensemble.",
        "measures": "Sequence-mapped alignment, RMSD, RMSF, frame-to-reference distance, deterministic clusters, and interpretable-frame populations.",
        "receives": ["StructQC evidence", "State ensemble summary", "Frame metrics", "Clusters and residue layers"],
        "non_claim": "Active-like or inactive-like describes structural resemblance, not biochemical activation, inhibition, or efficacy.",
    },
    "functional_site_state": {
        "use_when": "Use when curated residues, ligands, metals, or cofactors must be located and checked in one structure.",
        "measures": "Exact residue mapping, role-specific site completeness, ligand/cofactor observations, geometry, and method-specific pocket evidence.",
        "receives": ["StructQC evidence", "Site summary", "Residue and pocket tables", "Catalytic-competence evidence panels"],
        "non_claim": "A complete or plausible site is not observed catalysis, binding, inhibition, or physiological function.",
    },
    "assembly_interface": {
        "use_when": "Use when monomer measurements must be compared with a declared biological assembly or homologous architecture.",
        "measures": "Assembly identity, stoichiometry evidence, heavy-atom contacts, interface residues, SASA, burial, and lower-bound status.",
        "receives": ["StructQC evidence", "Assembly summary", "Interface table", "Interface and occlusion layers"],
        "non_claim": "Assembly geometry does not establish native abundance, intact-cell accessibility, or measured binding.",
    },
    "sf_csa": {
        "use_when": "Use for a checksum-pinned comparative search where structural and sequence relationships must remain separate.",
        "measures": "Foldseek structural hits, DIAMOND sequence hits, declared mechanism-group context, divergence, and closed-vocabulary relationship evidence.",
        "receives": ["StructQC evidence", "Separate structural and sequence result legs", "Comparison matrices", "Checksums and run manifest"],
        "non_claim": "Similarity, homology, or a shared fold never becomes exact identity or automatic functional transfer.",
    },
}


_ROLE_GUIDANCE: Mapping[str, Mapping[str, Any]] = {
    "structure": {"description": "The coordinate model to analyze.", "why_needed": "All residue-level measurements must bind to exact coordinate bytes.", "absence_effect": "blocked", "accepted_artifact_types": ["pdb.coordinates", "alphafold.model"], "source_ids": ["pdb", "alphafold_db"], "format_guide": "PDB text or PDBx/mmCIF with at least one model and polymer chain.", "validator_id": "coordinate_structure", "sensitivity": "sensitive_by_default"},
    "query_structure": {"description": "The coordinate model used as the structural-search query.", "why_needed": "SF-CSA must checksum and search the exact query model.", "absence_effect": "blocked", "accepted_artifact_types": ["pdb.coordinates", "alphafold.model"], "source_ids": ["pdb", "alphafold_db"], "format_guide": "PDB text or PDBx/mmCIF; the selected chain must match the query FASTA.", "validator_id": "coordinate_structure", "sensitivity": "sensitive_by_default"},
    "assembly": {"description": "A deposited or explicitly expanded biological assembly.", "why_needed": "Contacts and burial cannot be inferred from an isolated chain.", "absence_effect": "blocked", "accepted_artifact_types": ["pdb.biological_assembly"], "source_ids": ["pdb"], "format_guide": "Prefer RCSB biological-assembly mmCIF and record its assembly number.", "validator_id": "coordinate_structure", "sensitivity": "sensitive_by_default"},
    "reference_fasta": {"description": "The identity reference sequence for exact coordinate mapping.", "why_needed": "It distinguishes missing coordinates from a genuinely shorter sequence.", "absence_effect": "completeness_unevaluated", "accepted_artifact_types": ["uniprot.sequence"], "source_ids": ["uniprot_proteomes"], "format_guide": "One FASTA record with an accession-bearing header and standard amino-acid sequence.", "validator_id": "protein_fasta", "sensitivity": "sensitive_by_default"},
    "query_fasta": {"description": "The exact amino-acid sequence for the SF-CSA query accession.", "why_needed": "The sequence leg and coordinate identity checks depend on the same declared protein.", "absence_effect": "blocked", "accepted_artifact_types": ["uniprot.sequence"], "source_ids": ["uniprot_proteomes"], "format_guide": "FASTA containing the accession entered below; sequence and coordinates must map exactly.", "validator_id": "protein_fasta", "sensitivity": "sensitive_by_default"},
    "source_proteome": {"description": "A local FASTA defining the sequence-comparison universe.", "why_needed": "DIAMOND results are meaningful only relative to a declared, checksum-pinned proteome.", "absence_effect": "blocked", "accepted_artifact_types": ["uniprot.proteome"], "source_ids": ["uniprot_proteomes"], "format_guide": "Multi-record protein FASTA from one declared organism/proteome release.", "validator_id": "protein_fasta", "sensitivity": "sensitive_by_default"},
    "provenance": {"description": "An explicit declaration of where the coordinate model came from and how confidence is encoded.", "why_needed": "B-factors cannot be interpreted as pLDDT and models cannot be called experimental without it.", "absence_effect": "provenance_unknown", "accepted_artifact_types": [], "source_ids": [], "format_guide": "Use the template and retain unknown values when evidence is unavailable.", "template_id": "explicit_structure_provenance", "validator_id": "explicit_structure_provenance", "sensitivity": "non_sensitive_metadata"},
    "pae": {"description": "Pairwise predicted aligned error for an AlphaFold model.", "why_needed": "PAE reveals uncertain relative domain placement that pLDDT alone cannot describe.", "absence_effect": "domain_confidence_unevaluated", "accepted_artifact_types": ["alphafold.pae"], "source_ids": ["alphafold_db"], "format_guide": "AlphaFold PAE JSON containing a square residue-by-residue error matrix.", "validator_id": "predicted_aligned_error", "sensitivity": "sensitive_by_default"},
    "validation_report": {"description": "A checksum-bound wwPDB, MolProbity, or Phenix validation result.", "why_needed": "Coordinate parsing alone cannot assess clashes, stereochemistry, rotamers, or experiment-specific fit.", "absence_effect": "scientifically_incomplete", "accepted_artifact_types": ["wwpdb.validation"], "source_ids": ["wwpdb_validation"], "format_guide": "wwPDB validation XML or a supported normalized local validation JSON.", "validator_id": "structure_validation_report", "sensitivity": "sensitive_by_default"},
    "trajectory": {"description": "Optional time-ordered MD trajectory paired with the submitted topology.", "why_needed": "It enables per-frame distances, RMSF, clustering, and population accounting.", "absence_effect": "single_structure_only", "accepted_artifact_types": [], "source_ids": [], "format_guide": "MDAnalysis-compatible XTC, DCD, or TRR generated locally; record PBC handling and stride.", "validator_id": "molecular_trajectory", "sensitivity": "sensitive_by_default"},
    "topology_evidence": {"description": "Exact transmembrane helix spans bound to the coordinate checksum.", "why_needed": "Alpha-helical placement must derive its axis from declared membrane-spanning residues rather than a whole-structure barrel heuristic.", "absence_effect": "alpha_helical_placement_incomplete", "accepted_artifact_types": [], "source_ids": ["opm_ppm", "uniprot_proteomes"], "format_guide": "JSON with coordinate_sha256, source id/citation, exact chain and author-residue spans, and optional extracellular marker residue.", "template_id": "membrane_topology_evidence", "validator_id": "membrane_topology_evidence", "sensitivity": "non_sensitive_public_or_sensitive_local"},
    "active_reference": {"description": "At least two independently supported experimental active-state coordinate references.", "why_needed": "State labels are disabled without multiple references on both sides of the comparison.", "absence_effect": "blocked", "accepted_artifact_types": ["pdb.coordinates"], "source_ids": ["pdb", "sifts"], "format_guide": "Two or more experimental PDB/mmCIF chains plus their citation and evidence basis.", "validator_id": "experimental_state_reference", "sensitivity": "non_sensitive_public_or_sensitive_local"},
    "inactive_reference": {"description": "At least two independently supported experimental inactive-state coordinate references.", "why_needed": "State labels are disabled without multiple references on both sides of the comparison.", "absence_effect": "blocked", "accepted_artifact_types": ["pdb.coordinates"], "source_ids": ["pdb", "sifts"], "format_guide": "Two or more experimental PDB/mmCIF chains plus their citation and evidence basis.", "validator_id": "experimental_state_reference", "sensitivity": "non_sensitive_public_or_sensitive_local"},
    "alignment_map": {"description": "An exact query/reference residue-equivalence map in a declared coordinate system.", "why_needed": "It prevents extra construct domains and numbering coincidences from silently driving state RMSD.", "absence_effect": "blocked", "accepted_artifact_types": [], "source_ids": ["sifts", "uniprot_proteomes"], "format_guide": "Reference Set v2 JSON map keyed by UniProt positions with exact chain, author residue, insertion code, source citation, and source checksum.", "template_id": "state_alignment_map_v2", "validator_id": "state_alignment_map_v2", "sensitivity": "non_sensitive_public_or_sensitive_local"},
    "site_annotations": {"description": "Declared functional residues, roles, and source identifiers.", "why_needed": "The workbench never guesses catalytic residues from amino-acid type alone.", "absence_effect": "blocked", "accepted_artifact_types": ["mcsa.sites", "uniprot.annotations"], "source_ids": ["mcsa", "uniprot_proteomes"], "format_guide": "JSON/TSV/CSV keyed by chain and residue number with role and evidence source.", "template_id": "functional_site_annotations", "validator_id": "functional_site_annotations", "sensitivity": "non_sensitive_public_or_sensitive_local"},
    "uniprot_annotations": {"description": "Optional UniProt feature export for independent ActState evidence.", "why_needed": "It preserves annotation evidence separately from observed coordinate chemistry.", "absence_effect": "activity_annotation_leg_missing", "accepted_artifact_types": ["uniprot.annotations"], "source_ids": ["uniprot_proteomes"], "format_guide": "UniProt TSV with accession and feature columns such as active site, binding, site, and cofactor.", "validator_id": "uniprot_annotation_export", "sensitivity": "non_sensitive_public_or_sensitive_local"},
    "component_map": {"description": "Exact mappings between observed residue/component IDs and CCD or ChEBI identifiers.", "why_needed": "Unknown ligand synonyms remain unresolved rather than being accepted by name similarity.", "absence_effect": "cofactor_identity_unresolved", "accepted_artifact_types": [], "source_ids": ["pdb_ccd", "chebi"], "format_guide": "Use the JSON template; consult CCD and ChEBI records, but do not upload raw CIF, OBO, or bulk mapping files into this normalized role.", "template_id": "component_identifier_map", "validator_id": "component_identifier_map", "sensitivity": "non_sensitive_metadata"},
    "pocket_result": {"description": "An optional result exported by one named pocket tool.", "why_needed": "It adds method-specific pocket predictions without merging scores across tools.", "absence_effect": "pocket_evidence_not_run", "accepted_artifact_types": [], "source_ids": [], "format_guide": "Normalized JSON naming fpocket or P2Rank, runtime version, parameters, and unmodified method score.", "validator_id": "pocket_tool_result", "sensitivity": "sensitive_by_default"},
    "interpretation_tables": {"description": "Organism-appropriate, closed-vocabulary mechanism and divergence tables.", "why_needed": "Campaign-specific defaults must never leak into an unrelated comparative analysis.", "absence_effect": "blocked", "accepted_artifact_types": [], "source_ids": ["mcsa"], "format_guide": "JSON containing mechanism_families, contested_groups, divergence_sets, and classification_vocabulary arrays.", "template_id": "sf_csa_interpretation_tables", "validator_id": "sf_csa_interpretation_tables", "sensitivity": "non_sensitive_metadata"},
}


def analysis_definitions() -> list[dict[str, Any]]:
    """Return human-readable, source-aware task definitions for CLI and UI."""
    definitions = _base_analysis_definitions()
    for definition in definitions:
        definition.update(_WORKFLOW_GUIDANCE[definition["analysis_type"]])
        definition["scientific_readiness"] = {
            "software_state": definition["readiness"],
            "external_benchmark": "qualification_v2_incomplete",
            "release_gate": "all_five_mark_1_release_blocking_scopes_must_pass",
        }
        for input_role in definition["inputs"]:
            input_role["accepted_extensions"] = list(input_role["extensions"])
            input_role.update(_ROLE_GUIDANCE[input_role["role"]])
            input_role.setdefault("template_id", None)
    return definitions


def _runtime_version(binary: str) -> str | None:
    path = shutil.which(binary)
    if not path:
        return None
    for flag in ("version", "--version"):
        try:
            result = subprocess.run([path, flag], capture_output=True, text=True, timeout=10, check=False)
        except (OSError, subprocess.SubprocessError):
            continue
        value = (result.stdout or result.stderr).strip().splitlines()
        if value:
            return value[0][:160]
    return "available_unknown_version"


def tool_readiness(workspace: str | Path) -> list[dict[str, Any]]:
    root = Path(workspace).resolve()
    source_root = _checkout_source_root(root)
    runtimes = {
        "gemmi": "available" if importlib.util.find_spec("gemmi") else "missing",
        "mdanalysis": "available" if importlib.util.find_spec("MDAnalysis") else "missing",
        "freesasa": _runtime_version("freesasa") or "missing",
        "mkdssp": _runtime_version("mkdssp") or "missing",
        "foldseek": _runtime_version("foldseek") or "missing",
        "diamond": _runtime_version("diamond") or "missing",
        "molprobity_or_phenix": "available" if shutil.which("phenix.molprobity") or shutil.which("clashscore") else "missing",
        "vina": _runtime_version("vina") or "missing",
        "meeko": "available" if importlib.util.find_spec("meeko") else "missing",
        "haddock3": _runtime_version("haddock3") or "missing",
    }
    packages = ({
        "structure_qc": source_root / "structqc",
        "membrane_orientation": source_root / "Membrane Orientor" / "memorient",
        "conformational_state": source_root / "state-atlas",
        "functional_site_state": source_root / "site-context",
        "assembly_interface": source_root / "assembly-context",
        "sf_csa": source_root / "sf-csa",
    } if source_root is not None else {
        "structure_qc": _installed_package_root("structqc"),
        "membrane_orientation": _installed_package_root("memorient"),
        "conformational_state": _installed_package_root("state_atlas"),
        "functional_site_state": _installed_package_root("site_context"),
        "assembly_interface": _installed_package_root("assembly_context"),
        "sf_csa": _installed_package_root("sf_csa"),
    })
    required = {
        "structure_qc": [], "membrane_orientation": [], "conformational_state": [],
        "functional_site_state": [], "assembly_interface": [], "sf_csa": ["foldseek", "diamond"],
    }
    optional = {
        "structure_qc": ["gemmi", "molprobity_or_phenix", "mkdssp"],
        "membrane_orientation": [], "conformational_state": ["mdanalysis"],
        "functional_site_state": ["gemmi"], "assembly_interface": ["gemmi", "freesasa"],
        "sf_csa": [],
    }
    scopes: dict[str, list[dict[str, Any]]] = {
        "structure_qc": [{
            "scope_id": "coordinate_provenance_and_validation",
            "scientific_state": "prototype",
            "benchmark_collection": "qualification-v2-structqc",
            "release_blocking": True,
            "supported_subject_class": "experimental and predicted protein coordinate models",
            "required_evidence": ["exact coordinates", "explicit provenance", "method-appropriate validation evidence"],
            "known_limitations": ["Qualification v2 public panel is not complete."],
        }],
        "membrane_orientation": [
            {
                "scope_id": "beta_barrel",
                "scientific_state": "conditionally_qualified",
                "benchmark_collection": "qualification-v2-membrane-beta-barrel",
                "release_blocking": True,
                "supported_subject_class": "transmembrane beta-barrel proteins",
                "required_evidence": ["exact coordinates", "declared membrane context"],
                "known_limitations": ["Independent second-machine reproduction remains required."],
            },
            {
                "scope_id": "alpha_helical",
                "scientific_state": "prototype",
                "benchmark_collection": "qualification-v2-membrane-alpha-helical",
                "release_blocking": False,
                "supported_subject_class": "single-pass and multipass alpha-helical membrane proteins",
                "required_evidence": ["exact coordinates", "checksum-bound transmembrane spans", "external topology for sidedness"],
                "known_limitations": ["Experimental method; not part of the Mark 1 qualified scope."],
            },
        ],
        "conformational_state": [
            {
                "scope_id": "abl_family",
                "scientific_state": "prototype",
                "benchmark_collection": "qualification-v2-abl-state-atlas",
                "release_blocking": True,
                "supported_subject_class": "ABL-family experimental coordinate structures and ensembles",
                "required_evidence": ["two-sided experimental references", "exact SIFTS or explicit residue map", "at least 90 percent mapped ABL1 242-495 coverage"],
                "known_limitations": ["Qualification v2 held-out ABL gate has not yet passed."],
            },
            {
                "scope_id": "other_proteins",
                "scientific_state": "prototype",
                "benchmark_collection": "none",
                "release_blocking": False,
                "supported_subject_class": "other proteins with user-curated two-sided references",
                "required_evidence": ["two-sided experimental references", "exact declared alignment map"],
                "known_limitations": ["Outside the Mark 1 qualified scope."],
            },
        ],
        "functional_site_state": [{
            "scope_id": "curated_functional_site_mapping",
            "scientific_state": "prototype",
            "benchmark_collection": "qualification-v2-site-context",
            "release_blocking": True,
            "supported_subject_class": "proteins with curated residue-role annotations",
            "required_evidence": ["exact residue mapping", "role-specific annotations", "exact component identifiers when applicable"],
            "known_limitations": ["Qualification v2 M-CSA panel is not complete."],
        }],
        "assembly_interface": [{
            "scope_id": "deposited_biological_assembly",
            "scientific_state": "prototype",
            "benchmark_collection": "qualification-v2-assembly-context",
            "release_blocking": True,
            "supported_subject_class": "deposited protein biological assemblies",
            "required_evidence": ["exact assembly operators", "entity and copy identity", "pinned FreeSASA runtime for qualification"],
            "known_limitations": ["Qualification v2 operator and FreeSASA panel is not complete."],
        }],
        "sf_csa": [{
            "scope_id": "curated_structure_sequence_comparison",
            "scientific_state": "prototype",
            "benchmark_collection": "qualification-v2-sf-csa",
            "release_blocking": True,
            "supported_subject_class": "curator-declared protein family comparison panels",
            "required_evidence": ["checksum-pinned Foldseek database", "declared proteome", "closed interpretation tables", "Foldseek and DIAMOND runtimes"],
            "known_limitations": ["Qualification v2 four-family panel is not complete."],
        }],
    }
    rows = []
    for definition in analysis_definitions():
        tool = definition["analysis_type"]
        missing = [name for name in required[tool] if runtimes[name] == "missing"]
        state = "blocked_missing_runtime" if missing else definition["readiness"]
        rows.append({
            "contract_id": "tool_readiness_record", "analysis_type": tool,
            "title": definition["title"], "state": state,
            "package_source_sha256": _tree_sha(packages[tool]) if packages[tool] is not None else None,
            "required_runtimes": {name: runtimes[name] for name in required[tool]},
            "optional_runtimes": {name: runtimes[name] for name in optional[tool]},
            "unit_tests": {"state": "not_run_in_this_process", "separate_from_scientific_validation": True},
            "external_benchmark": {"state": "qualification_v2_incomplete", "collection": "qualification-v2"},
            "scientific_scopes": scopes[tool],
            "version_control": "missing" if not (root / ".git").exists() else "available",
            "limitations": [definition["claim_ceiling"]],
        })
    rows.extend([
        {"contract_id": "tool_readiness_record", "analysis_type": "structural_conservation", "title": "Structural conservation", "state": "prototype", "labs": True, "limitations": ["Mapping utility; no alignment generation or public benchmark pack."]},
        {"contract_id": "tool_readiness_record", "analysis_type": "structural_relationships", "title": "Structural relationships", "state": "adapter_only", "labs": True, "limitations": ["Consumes supplied Foldseek/SIFTS evidence; it is not a structural search engine."]},
        {"contract_id": "tool_readiness_record", "analysis_type": "docking", "title": "Docking", "state": "blocked_missing_runtime", "labs": True, "required_runtimes": {"vina": runtimes["vina"], "meeko": runtimes["meeko"], "haddock3": runtimes["haddock3"]}, "limitations": ["Readiness only. Docking score is not affinity or activity."]},
        {"contract_id": "tool_readiness_record", "analysis_type": "design", "title": "Structural design", "state": "adapter_only", "labs": True, "limitations": ["Readiness only. Generated candidates are hypotheses and are never released automatically."]},
    ])
    return rows


class StructuralAnalysisStore:
    """Immutable analysis cases and deterministic output bundles."""

    def __init__(self, workspace: str | Path) -> None:
        self.workspace = Path(workspace).resolve()
        checkout = _checkout_source_root(self.workspace)
        self.installed_mode = checkout is None
        self.source_root = checkout or self.workspace
        self.root = self.workspace / "structural_analyses"
        self.cases_root = self.root / "cases"
        self.objects_root = self.root / "objects" / "sha256"
        self.ingests_root = self.root / "ingests"
        for directory in (self.cases_root, self.objects_root, self.ingests_root):
            directory.mkdir(parents=True, exist_ok=True)
        self._lock = threading.RLock()
        self._definitions = {d["analysis_type"]: d for d in analysis_definitions()}

    def _case_dir(self, analysis_id: str) -> Path:
        if not SAFE_ID.fullmatch(analysis_id):
            raise AnalysisError("analysis id must use lowercase letters, digits, hyphens, or underscores")
        return self.cases_root / analysis_id

    def _manifest_path(self, analysis_id: str) -> Path:
        return self._case_dir(analysis_id) / "ANALYSIS_CASE.json"

    def load(self, analysis_id: str) -> dict[str, Any]:
        try:
            value = json.loads(self._manifest_path(analysis_id).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise AnalysisError(f"analysis does not exist or is unreadable: {analysis_id}") from exc
        return value

    def list(self) -> list[dict[str, Any]]:
        rows = []
        for path in sorted(self.cases_root.glob("*/ANALYSIS_CASE.json")):
            try:
                value = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            rows.append({key: value.get(key) for key in ("analysis_id", "analysis_type", "question", "subject_id", "revision", "state", "latest_run_id")})
        return rows

    def create(self, analysis_id: str, *, analysis_type: str, question: str, subject_id: str = "") -> dict[str, Any]:
        if analysis_type not in self._definitions:
            raise AnalysisError(f"unknown structural analysis type: {analysis_type}")
        if not question.strip():
            raise AnalysisError("research question is required")
        directory = self._case_dir(analysis_id)
        with self._lock:
            if directory.exists():
                raise AnalysisError(f"analysis already exists: {analysis_id}")
            directory.mkdir(parents=True)
            manifest = {
                "schema_version": SCHEMA_VERSION, "contract_id": "structural_analysis_case",
                "analysis_id": analysis_id, "analysis_type": analysis_type,
                "question": question.strip(), "subject_id": subject_id.strip() or analysis_id,
                "revision": 1, "state": "draft", "inputs": [], "source_adoptions": [], "parameters": {}, "runs": [],
                "latest_run_id": None,
            }
            self._commit(directory, manifest)
            return manifest

    def _commit(self, directory: Path, manifest: dict[str, Any]) -> None:
        material = {k: v for k, v in manifest.items() if k != "revision_sha256"}
        manifest["revision_sha256"] = _sha_bytes(_canonical(material))
        _write_json(directory / "ANALYSIS_CASE.json", manifest)

    def update_parameters(self, analysis_id: str, parameters: Mapping[str, Any]) -> dict[str, Any]:
        manifest = self.load(analysis_id)
        allowed = {p["name"]: p for p in self._definitions[manifest["analysis_type"]].get("parameters", [])}
        cleaned = {}
        for key, value in parameters.items():
            if key not in allowed:
                raise AnalysisError(f"unknown parameter for {manifest['analysis_type']}: {key}")
            cleaned[key] = value
        manifest["parameters"] = cleaned
        manifest["revision"] = int(manifest["revision"]) + 1
        manifest["state"] = "draft"
        self._commit(self._case_dir(analysis_id), manifest)
        return manifest

    def begin_ingest(self, analysis_id: str, *, role: str, file_name: str, size: int, expected_sha256: str) -> dict[str, Any]:
        manifest = self.load(analysis_id)
        definition = self._definitions[manifest["analysis_type"]]
        role_def = next((item for item in definition["inputs"] if item["role"] == role), None)
        if role_def is None:
            raise AnalysisError(f"input role is not allowed for {manifest['analysis_type']}: {role}")
        suffix = Path(file_name).suffix.lower()
        if suffix not in role_def["extensions"]:
            raise AnalysisError(f"{role} requires one of: {', '.join(role_def['extensions'])}")
        if size < 1 or size > MAX_INPUT_BYTES:
            raise AnalysisError("input file is empty or exceeds the role-independent 20 GiB limit")
        if expected_sha256 and not SHA256.fullmatch(expected_sha256):
            raise AnalysisError("expected_sha256 must be empty or a lowercase SHA-256 digest")
        if not role_def.get("multiple") and any(item["role"] == role for item in manifest["inputs"]):
            raise AnalysisError(f"input role accepts only one file: {role}")
        upload_id = "upload_" + os.urandom(12).hex()
        directory = self.ingests_root / upload_id
        directory.mkdir()
        record = {
            "schema_version": SCHEMA_VERSION, "contract_id": "file_ingest_manifest",
            "upload_id": upload_id, "analysis_id": analysis_id, "role": role,
            "file_name": Path(file_name).name, "size": int(size),
            "expected_sha256": expected_sha256 or None,
            "checksum_policy": "client_and_server" if expected_sha256 else "server_finalized",
            "received": 0, "next_chunk": 0, "state": "receiving",
        }
        _write_json(directory / "INGEST.json", record)
        return record

    def ingest_chunk(self, upload_id: str, index: int, content: bytes) -> dict[str, Any]:
        if not SAFE_UPLOAD_ID.fullmatch(upload_id):
            raise AnalysisError("invalid upload id")
        directory = self.ingests_root / upload_id
        try:
            record = json.loads((directory / "INGEST.json").read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise AnalysisError("upload does not exist") from exc
        if record["state"] != "receiving" or index != int(record["next_chunk"]):
            raise AnalysisError("chunks must be uploaded once and in order")
        if not content or len(content) > MAX_CHUNK_BYTES:
            raise AnalysisError("chunk is empty or exceeds 8 MiB")
        if int(record["received"]) + len(content) > int(record["size"]):
            raise AnalysisError("chunk exceeds the declared file size")
        (directory / f"chunk-{index:08d}").write_bytes(content)
        record["received"] = int(record["received"]) + len(content)
        record["next_chunk"] = index + 1
        _write_json(directory / "INGEST.json", record)
        return record

    def finalize_ingest(self, upload_id: str) -> dict[str, Any]:
        if not SAFE_UPLOAD_ID.fullmatch(upload_id):
            raise AnalysisError("invalid upload id")
        directory = self.ingests_root / upload_id
        record = json.loads((directory / "INGEST.json").read_text(encoding="utf-8"))
        if int(record["received"]) != int(record["size"]):
            raise AnalysisError("upload is incomplete")
        digest = hashlib.sha256()
        chunks = sorted(directory.glob("chunk-*"))
        for chunk in chunks:
            digest.update(chunk.read_bytes())
        actual = digest.hexdigest()
        if record.get("expected_sha256") and actual != record["expected_sha256"]:
            raise AnalysisError("uploaded bytes do not match the declared checksum")
        record["actual_sha256"] = actual
        self._validate_ingested_content(record, chunks)
        object_path = self.objects_root / actual[:2] / actual
        object_path.parent.mkdir(parents=True, exist_ok=True)
        if not object_path.exists():
            temporary = object_path.with_suffix(".partial")
            with temporary.open("wb") as output:
                for chunk in chunks:
                    output.write(chunk.read_bytes())
            os.replace(temporary, object_path)
        manifest = self.load(record["analysis_id"])
        item = {
            "contract_id": "analysis_input_reference", "role": record["role"],
            "file_name": record["file_name"], "sha256": actual, "bytes": record["size"],
        }
        manifest["inputs"].append(item)
        manifest["inputs"].sort(key=lambda value: (value["role"], value["sha256"], value["file_name"]))
        manifest["revision"] = int(manifest["revision"]) + 1
        manifest["state"] = "draft"
        self._commit(self._case_dir(record["analysis_id"]), manifest)
        shutil.rmtree(directory)
        return item

    @staticmethod
    def _validate_ingested_content(record: Mapping[str, Any], chunks: list[Path]) -> None:
        """Reject obvious extension/content mismatches without claiming full parsing."""
        prefix = b"".join(path.read_bytes() for path in chunks[:1])[:131072]
        suffix = Path(str(record.get("file_name", ""))).suffix.lower()
        stripped = prefix.lstrip()
        text = prefix.decode("utf-8", errors="replace")
        mismatch = False
        if suffix == ".pdb":
            mismatch = not any(line.startswith(("ATOM", "HETATM", "HEADER", "MODEL", "REMARK", "TITLE")) for line in text.splitlines()[:200])
        elif suffix in {".cif", ".mmcif"}:
            mismatch = not stripped.lower().startswith(b"data_")
        elif suffix in {".fasta", ".fa", ".faa"}:
            mismatch = not stripped.startswith(b">")
        elif suffix == ".xml":
            mismatch = not stripped.startswith(b"<")
        elif suffix == ".json":
            try: json.loads(b"".join(path.read_bytes() for path in chunks))
            except (json.JSONDecodeError, UnicodeDecodeError): mismatch = True
        elif suffix == ".tsv":
            mismatch = "\t" not in text.splitlines()[0] if text.splitlines() else True
        elif suffix == ".csv":
            mismatch = "," not in text.splitlines()[0] if text.splitlines() else True
        if mismatch:
            raise AnalysisError(f"uploaded content does not match declared {suffix or 'file'} format")

    def add_file(self, analysis_id: str, *, role: str, path: str | Path) -> dict[str, Any]:
        """Ingest one local CLI-selected file through the same bounded contract as the UI."""
        source = Path(path).resolve()
        if not source.is_file():
            raise AnalysisError(f"input file is unavailable: {source}")
        size = source.stat().st_size
        digest = _sha_file(source)
        ingest = self.begin_ingest(
            analysis_id, role=role, file_name=source.name, size=size,
            expected_sha256=digest,
        )
        with source.open("rb") as handle:
            index = 0
            while chunk := handle.read(MAX_CHUNK_BYTES):
                self.ingest_chunk(ingest["upload_id"], index, chunk)
                index += 1
        return self.finalize_ingest(ingest["upload_id"])

    def adopt_source_artifact(
        self,
        analysis_id: str,
        *,
        role: str,
        path: str | Path,
        artifact_manifest: Mapping[str, Any],
    ) -> dict[str, Any]:
        """Adopt a verified cached public artifact into a new case revision.

        Acquisition and adoption are deliberately separate.  A refresh can add
        bytes to the shared cache, but only this explicit action binds them to a
        scientific analysis.
        """
        if artifact_manifest.get("contract_id") != "source_artifact_manifest":
            raise AnalysisError("source artifact manifest has the wrong contract")
        source_path = Path(path)
        expected = str(artifact_manifest.get("sha256", ""))
        if not SHA256.fullmatch(expected) or not source_path.is_file() or _sha_file(source_path) != expected:
            raise AnalysisError("source artifact bytes do not match their manifest")
        item = self.add_file(analysis_id, role=role, path=source_path)
        manifest = self.load(analysis_id)
        bound = next(
            (row for row in manifest["inputs"] if row["role"] == role and row["sha256"] == item["sha256"]),
            None,
        )
        if bound is None:
            raise AnalysisError("adopted artifact was not recorded in the analysis")
        bound["source_artifact"] = {
            key: artifact_manifest.get(key)
            for key in (
                "source_id", "provider", "artifact_type", "identifier", "release",
                "license", "origin", "manifest_sha256",
            )
        }
        manifest["revision"] = int(manifest["revision"]) + 1
        adoption = {
            "schema_version": SCHEMA_VERSION, "contract_id": "source_adoption_record",
            "analysis_id": analysis_id, "role": role, "analysis_revision": manifest["revision"],
            "artifact_manifest_sha256": artifact_manifest.get("manifest_sha256"),
            "artifact_sha256": expected,
        }
        manifest.setdefault("source_adoptions", []).append(adoption)
        manifest["source_adoptions"].sort(key=lambda row: (row["analysis_revision"], row["role"], row["artifact_sha256"]))
        manifest["state"] = "draft"
        self._commit(self._case_dir(analysis_id), manifest)
        return bound

    def object_path(self, digest: str) -> Path:
        if not SHA256.fullmatch(digest):
            raise AnalysisError("invalid object checksum")
        path = self.objects_root / digest[:2] / digest
        if not path.is_file():
            raise AnalysisError("content-addressed object is unavailable")
        return path

    def _inputs_by_role(self, manifest: Mapping[str, Any]) -> dict[str, list[dict[str, Any]]]:
        found: dict[str, list[dict[str, Any]]] = {}
        for item in manifest.get("inputs", []):
            found.setdefault(str(item["role"]), []).append(dict(item))
        return found

    def _materialize_inputs(self, inputs: Mapping[str, list[dict[str, Any]]], run_dir: Path) -> dict[str, list[dict[str, Any]]]:
        """Give CLI tools safe case-local names while objects remain addressed by digest."""
        root = run_dir / "generated" / "inputs"
        root.mkdir(parents=True, exist_ok=True)
        materialized: dict[str, list[dict[str, Any]]] = {}
        for role, items in sorted(inputs.items()):
            for index, item in enumerate(items, 1):
                suffix = Path(str(item.get("file_name", ""))).suffix.lower()
                target = root / f"{role}-{index:03d}-{item['sha256'][:12]}{suffix}"
                if not target.is_file():
                    shutil.copyfile(self.object_path(item["sha256"]), target)
                copied = dict(item); copied["materialized_path"] = str(target)
                materialized.setdefault(role, []).append(copied)
        return materialized

    def _sf_pack(self, pack_id: str) -> tuple[dict[str, Any], Path]:
        if not SAFE_ID.fullmatch(pack_id):
            raise AnalysisError("SF-CSA database pack ID is invalid")
        pack_root = self.workspace / ".yauvi-cache" / "sources" / "packs" / pack_id
        pack_path = pack_root / "PACK_MANIFEST.json"
        try:
            pack = json.loads(pack_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise AnalysisError("pinned SF-CSA pack manifest is unavailable or unreadable") from exc
        relative = str(pack.get("sf_csa_database_manifest") or pack.get("database_manifest") or "")
        if not relative or Path(relative).is_absolute() or ".." in Path(relative).parts:
            raise AnalysisError("SF-CSA pack does not declare a safe database manifest")
        database = (pack_root / relative).resolve()
        try: database.relative_to(pack_root.resolve())
        except ValueError as exc: raise AnalysisError("SF-CSA database manifest escapes its pack") from exc
        expected = str(pack.get("database_manifest_sha256") or "")
        files = pack.get("files", {})
        if not expected and isinstance(files, Mapping): expected = str(files.get(relative, ""))
        if not expected and isinstance(files, list):
            expected = str(next((item.get("sha256") for item in files if isinstance(item, Mapping) and item.get("path") == relative), ""))
        if not SHA256.fullmatch(expected) or not database.is_file() or _sha_file(database) != expected:
            raise AnalysisError("SF-CSA database manifest checksum is missing or mismatched")
        return pack, database

    def _sf_interpretation(self, inputs: Mapping[str, list[dict[str, Any]]]) -> dict[str, Any]:
        source = self._one(inputs, "interpretation_tables")
        try: value = json.loads(source.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc: raise AnalysisError("SF-CSA interpretation tables are not valid JSON") from exc
        required = ("mechanism_families", "contested_groups", "divergence_sets", "classification_vocabulary")
        missing = [key for key in required if not isinstance(value.get(key), list)]
        if missing: raise AnalysisError("organism-neutral SF-CSA interpretation tables are incomplete: " + ", ".join(missing))
        return value

    def _membrane_topology(self, inputs: Mapping[str, list[dict[str, Any]]]) -> tuple[Path | None, dict[str, Any] | None]:
        source = self._one(inputs, "topology_evidence", required=False)
        if source is None:
            return None, None
        try:
            document = json.loads(source.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise AnalysisError("membrane topology evidence is not valid JSON") from exc
        structure_rows = inputs.get("structure", [])
        coordinate_sha = str(document.get("coordinate_sha256", ""))
        if not structure_rows or coordinate_sha != str(structure_rows[0].get("sha256", "")):
            raise AnalysisError("membrane topology evidence is not bound to the submitted coordinate checksum")
        source_record = document.get("source", {})
        if not isinstance(source_record, Mapping) or not str(source_record.get("id", "")).strip() or not str(source_record.get("citation", "")).strip():
            raise AnalysisError("membrane topology evidence requires a source id and citation")
        spans = document.get("spans", [])
        if not isinstance(spans, list) or not spans:
            raise AnalysisError("membrane topology evidence has no declared transmembrane spans")
        declared: set[tuple[str, int, str]] = set()
        for span_number, span in enumerate(spans, 1):
            if not isinstance(span, Mapping):
                raise AnalysisError(f"membrane topology span {span_number} is not an object")
            keys: list[tuple[str, int, str]] = []
            residues = span.get("residues")
            if isinstance(residues, list) and residues:
                for residue in residues:
                    if not isinstance(residue, Mapping):
                        raise AnalysisError(f"membrane topology span {span_number} has a malformed residue")
                    try:
                        key = (
                            str(residue.get("chain_id", "")), int(residue.get("auth_seq_id")),
                            str(residue.get("insertion_code", "")),
                        )
                    except (TypeError, ValueError) as exc:
                        raise AnalysisError(f"membrane topology span {span_number} has a malformed residue") from exc
                    if not key[0]:
                        raise AnalysisError(f"membrane topology span {span_number} has a residue without a chain")
                    keys.append(key)
            else:
                chain = str(span.get("chain_id", ""))
                try:
                    start = int(span.get("start_auth_seq_id")); end = int(span.get("end_auth_seq_id"))
                except (TypeError, ValueError) as exc:
                    raise AnalysisError(f"membrane topology span {span_number} has an invalid range") from exc
                if not chain or start > end:
                    raise AnalysisError(f"membrane topology span {span_number} has an invalid chain or range")
                keys = [(chain, position, "") for position in range(start, end + 1)]
            if len(keys) < 6:
                raise AnalysisError(f"membrane topology span {span_number} maps fewer than six residues")
            if len(set(keys)) != len(keys) or declared.intersection(keys):
                raise AnalysisError(f"membrane topology span {span_number} contains duplicate or overlapping residues")
            declared.update(keys)
        structure_item = structure_rows[0]
        structure_path = (
            Path(structure_item["materialized_path"])
            if structure_item.get("materialized_path") else self.object_path(structure_item["sha256"])
        )
        try:
            import gemmi  # type: ignore
            structure_text = structure_path.read_text(encoding="utf-8")
            suffix = Path(str(structure_item.get("file_name", ""))).suffix.lower()
            if suffix == ".pdb":
                coordinate = gemmi.read_pdb_string(structure_text)
            elif suffix in {".cif", ".mmcif"}:
                coordinate = gemmi.make_structure_from_block(gemmi.cif.read_string(structure_text).sole_block())
            else:
                raise AnalysisError("membrane topology coordinate format is unsupported")
            model = coordinate[0]
            coordinate_keys = [
                (str(chain.name), int(residue.seqid.num), str(residue.seqid.icode).strip())
                for chain in model for residue in chain
                if any(str(atom.name).strip() == "CA" for atom in residue)
            ]
        except (ImportError, IndexError, RuntimeError, OSError, UnicodeError) as exc:
            raise AnalysisError(f"membrane topology coordinate mapping cannot be validated: {exc}") from exc
        if len(set(coordinate_keys)) != len(coordinate_keys):
            raise AnalysisError("membrane topology coordinate model has duplicate residue identities")
        missing = sorted(declared - set(coordinate_keys))
        if missing:
            preview = ", ".join(f"{chain}:{position}{icode}" for chain, position, icode in missing[:8])
            raise AnalysisError(f"membrane topology residues are absent from the coordinate model: {preview}")
        return source, document

    def _state_alignment(self, inputs: Mapping[str, list[dict[str, Any]]]) -> tuple[Path, dict[str, Any], str]:
        source = self._one(inputs, "alignment_map")
        try:
            document = json.loads(source.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise AnalysisError("state alignment map is not valid JSON") from exc
        if str(document.get("coordinate_system", "")) != "uniprot":
            raise AnalysisError("state alignment map coordinate_system must be uniprot")
        domain = document.get("domain", {})
        try:
            start, end = int(domain.get("uniprot_start")), int(domain.get("uniprot_end"))
        except (AttributeError, TypeError, ValueError) as exc:
            raise AnalysisError("state alignment map requires integer UniProt domain boundaries") from exc
        if (start, end) != (242, 495):
            raise AnalysisError("Mark 1 ABL analysis requires the curated UniProt ABL1 242-495 alignment mask")
        source_record = document.get("source", {})
        if not isinstance(source_record, Mapping) or not str(source_record.get("id", "")).strip() or not str(source_record.get("citation", "")).strip():
            raise AnalysisError("state alignment map requires a source id and citation")
        if not SHA256.fullmatch(str(source_record.get("sha256", ""))):
            raise AnalysisError("state alignment map source.sha256 must be a SHA-256 digest")
        expected_positions = end - start + 1
        required_exact = int(expected_positions * 0.9 + 0.999999)

        def exact_rows(rows: Any, label: str) -> tuple[set[int], set[str]]:
            if not isinstance(rows, list):
                raise AnalysisError(f"state alignment map {label} rows must be a list")
            positions: set[int] = set(); chains: set[str] = set()
            for row in rows:
                if not isinstance(row, Mapping) or str(row.get("mapping_state", "exact")) != "exact":
                    continue
                try:
                    position = int(row.get("uniprot_position")); int(row.get("auth_seq_id"))
                except (TypeError, ValueError) as exc:
                    raise AnalysisError(f"state alignment map {label} has a malformed exact row") from exc
                chain = str(row.get("chain_id", ""))
                if start <= position <= end and chain:
                    if position in positions:
                        raise AnalysisError(f"state alignment map {label} repeats UniProt position {position}")
                    positions.add(position); chains.add(chain)
            if len(positions) < required_exact:
                raise AnalysisError(
                    f"state alignment map {label} exact coverage is below 90 percent of ABL1 242-495"
                )
            return positions, chains

        _query_positions, query_chains = exact_rows(document.get("query"), "query")
        if len(query_chains) != 1:
            raise AnalysisError("state alignment map query must resolve to exactly one chain")
        reference_rows = document.get("references", {})
        metadata = document.get("reference_metadata", {})
        if not isinstance(reference_rows, Mapping) or not isinstance(metadata, Mapping):
            raise AnalysisError("state alignment map requires reference rows and reference_metadata keyed by reference ID")
        expected_ids = [
            *(f"ACTIVE_{index}" for index in range(1, len(inputs.get("active_reference", [])) + 1)),
            *(f"INACTIVE_{index}" for index in range(1, len(inputs.get("inactive_reference", [])) + 1)),
        ]
        pdb_entries: set[str] = set()
        for reference_id in expected_ids:
            _positions, chains = exact_rows(reference_rows.get(reference_id), f"reference {reference_id}")
            meta = metadata.get(reference_id, {})
            pdb_entry = str(meta.get("pdb_entry_id", "")).strip().upper() if isinstance(meta, Mapping) else ""
            chain = str(meta.get("chain_id", "")).strip() if isinstance(meta, Mapping) else ""
            if not pdb_entry or not chain:
                raise AnalysisError(f"state alignment map reference_metadata is incomplete for {reference_id}")
            if chain not in chains:
                raise AnalysisError(f"state alignment map chain metadata conflicts with rows for {reference_id}")
            if pdb_entry in pdb_entries:
                raise AnalysisError(f"state alignment map repeats PDB entry {pdb_entry}")
            pdb_entries.add(pdb_entry)
        return source, document, next(iter(query_chains))

    def preflight(self, analysis_id: str) -> dict[str, Any]:
        manifest = self.load(analysis_id)
        definition = self._definitions[manifest["analysis_type"]]
        inputs = self._inputs_by_role(manifest)
        checks = []
        for role in definition["inputs"]:
            present = inputs.get(role["role"], [])
            minimum_files = int(role.get("minimum_files", 1 if role["required"] else 0))
            ok = len(present) >= minimum_files
            detail = f"{len(present)} checksum-verified file(s)" if present else (
                f"optional; missing effect: {role.get('absence_effect', 'reduced evidence')}" if ok
                else f"required input is missing; analysis cannot run ({role.get('absence_effect', 'blocked')})"
            )
            if present and not ok:
                detail = f"{len(present)} checksum-verified file(s); at least {minimum_files} required"
            checks.append({"name": f"input:{role['role']}", "category": "files", "ok": ok, "detail": detail})
        parameters = manifest.get("parameters", {})
        for spec in definition.get("parameters", []):
            if spec.get("required"):
                value = parameters.get(spec["name"], spec.get("default"))
                category = "identity_mapping" if spec["name"] in {"chain", "subject_chain", "accession", "organism"} else "scientific_evidence"
                checks.append({"name": f"parameter:{spec['name']}", "category": category, "ok": value not in (None, ""),
                               "detail": "recorded" if value not in (None, "") else "required parameter is missing"})
        readiness = next(item for item in tool_readiness(self.workspace) if item["analysis_type"] == manifest["analysis_type"])
        checks.append({"name": "runtime:registered_tool", "category": "runtimes", "ok": readiness["state"] != "blocked_missing_runtime",
                       "detail": readiness["state"]})
        if manifest["analysis_type"] == "conformational_state" and inputs.get("trajectory"):
            checks.append({"name": "runtime:mdanalysis", "category": "runtimes", "ok": importlib.util.find_spec("MDAnalysis") is not None,
                           "detail": "required for trajectory input"})
        if manifest["analysis_type"] == "membrane_orientation":
            alpha_scope = str(parameters.get("context", "")) in {"eukaryotic_pm", "tm_receptor"}
            topology = None
            try:
                _topology_path, topology = self._membrane_topology(inputs)
                topology_ok, topology_detail = True, "coordinate-bound transmembrane spans verified" if topology else "not supplied"
            except AnalysisError as exc:
                topology_ok, topology_detail = False, str(exc)
            if alpha_scope and topology is None and topology_ok:
                topology_ok = False
                topology_detail = "alpha-helical analysis requires checksum-bound transmembrane spans"
            checks.append({"name": "evidence:membrane_topology", "category": "scientific_evidence",
                           "ok": topology_ok, "detail": topology_detail})
        if manifest["analysis_type"] == "conformational_state":
            family_ok = str(parameters.get("subject_family", "")) == "ABL1"
            checks.append({"name": "scope:abl_family", "category": "scientific_evidence", "ok": family_ok,
                           "detail": "Mark 1 qualified scope" if family_ok else "other protein families are prototype-only"})
            try:
                thresholds_ok = (
                    float(parameters.get("max_rmsd_A", 2.5)) == 2.5
                    and float(parameters.get("min_margin_A", 0.25)) == 0.25
                )
            except (TypeError, ValueError):
                thresholds_ok = False
            checks.append({"name": "evidence:abl_frozen_decision_rules", "category": "scientific_evidence",
                           "ok": thresholds_ok,
                           "detail": "RMSD 2.5 A and opposite-state margin 0.25 A"
                           if thresholds_ok else "Mark 1 ABL thresholds cannot be changed"})
            try:
                _map_path, _map_document, map_chain = self._state_alignment(inputs)
                map_ok, map_detail = True, f"exact ABL1 242-495 map verified for query chain {map_chain}"
            except AnalysisError as exc:
                map_ok, map_detail = False, str(exc)
            checks.append({"name": "evidence:state_alignment_map_v2", "category": "identity_mapping",
                           "ok": map_ok, "detail": map_detail})
        if manifest["analysis_type"] == "sf_csa":
            pack = str(parameters.get("database_pack", ""))
            try:
                self._sf_pack(pack)
                pack_ok, pack_detail = True, f"{pack}: checksum verified"
            except AnalysisError as exc:
                pack_ok, pack_detail = False, str(exc)
            checks.append({"name": "reference:sf_csa_database_pack", "category": "scientific_evidence", "ok": pack_ok, "detail": pack_detail})
            try:
                interpretation = self._sf_interpretation(inputs)
                declared_group = str(parameters.get("mechanism_group", ""))
                groups = {str(item.get("group", "")) for item in interpretation["mechanism_families"] if isinstance(item, Mapping)}
                if declared_group not in groups:
                    raise AnalysisError("declared mechanism group is absent from the supplied interpretation table")
                interpretation_ok, interpretation_detail = True, "closed organism-neutral tables supplied"
            except AnalysisError as exc:
                interpretation_ok, interpretation_detail = False, str(exc)
            checks.append({"name": "evidence:sf_csa_interpretation_tables", "category": "scientific_evidence", "ok": interpretation_ok,
                           "detail": interpretation_detail})
        valid = all(item["ok"] for item in checks)
        record = {
            "schema_version": SCHEMA_VERSION, "contract_id": "analysis_preflight",
            "analysis_id": analysis_id, "analysis_type": manifest["analysis_type"],
            "revision_sha256": manifest["revision_sha256"], "valid": valid, "checks": checks,
            "readiness": readiness,
        }
        _write_json(self._case_dir(analysis_id) / "PREFLIGHT.json", record)
        manifest["state"] = "validated" if valid else "blocked"
        self._commit(self._case_dir(analysis_id), manifest)
        return record

    def _one(self, inputs: Mapping[str, list[dict[str, Any]]], role: str, required: bool = True) -> Path | None:
        values = inputs.get(role, [])
        if not values:
            if required:
                raise AnalysisError(f"required input is missing: {role}")
            return None
        return Path(values[0]["materialized_path"]) if values[0].get("materialized_path") else self.object_path(values[0]["sha256"])

    def _package_env(self) -> dict[str, str]:
        if self.installed_mode:
            return os.environ.copy()
        paths = [
            self.source_root / "structqc" / "src", self.source_root / "assembly-context" / "src",
            self.source_root / "site-context" / "src", self.source_root / "state-atlas" / "src",
            self.source_root / "activity-state" / "src", self.source_root / "sf-csa" / "src",
            self.source_root / "Membrane Orientor" / "memorient" / "src",
        ]
        env = os.environ.copy()
        existing = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = os.pathsep.join([*(str(p) for p in paths), *([existing] if existing else [])])
        return env

    def _execute(self, command: list[str], *, cwd: Path, log_path: Path,
                 cancel_event: Any | None = None, on_process: Any | None = None) -> int:
        process = subprocess.Popen(command, cwd=cwd, env=self._package_env(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if on_process is not None:
            on_process(process)
        while True:
            try:
                stdout, stderr = process.communicate(timeout=0.1)
                break
            except subprocess.TimeoutExpired:
                pass
            if cancel_event is not None and cancel_event.is_set():
                process.terminate()
                try: process.wait(timeout=5)
                except subprocess.TimeoutExpired: process.kill()
                stdout, stderr = process.communicate()
                log_path.parent.mkdir(parents=True, exist_ok=True)
                log_path.write_text((stdout + "\n" + stderr).replace(str(self.workspace), "<workspace>")
                                    .replace(str(self.source_root), "<source-tree>"), encoding="utf-8")
                raise InterruptedError("analysis cancellation was requested")
        scrubbed = (stdout + ("\n" if stdout and stderr else "") + stderr)
        scrubbed = scrubbed.replace(str(self.workspace), "<workspace>").replace(str(self.root), "<analysis-store>")
        scrubbed = scrubbed.replace(str(self.source_root), "<source-tree>")
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(scrubbed, encoding="utf-8")
        return int(process.returncode)

    def _run_structqc(self, manifest: Mapping[str, Any], inputs: Mapping[str, list[dict[str, Any]]], run_dir: Path,
                      cancel_event: Any | None = None, on_process: Any | None = None) -> tuple[Path, int]:
        output = run_dir / "outputs" / "structqc"
        structure = self._one(inputs, "structure")
        command = [sys.executable, "-m", "structqc.cli", "run", "--structure", str(structure),
                   "--subject-id", str(manifest["subject_id"]), "--out", str(output)]
        for role, flag in (("reference_fasta", "--reference-fasta"), ("provenance", "--provenance"),
                           ("pae", "--pae"), ("validation_report", "--validation-report")):
            value = self._one(inputs, role, required=False)
            if value:
                command.extend((flag, str(value)))
        command.append("--require-external-validation")
        params = manifest.get("parameters", {})
        if params.get("model") not in (None, ""):
            command.extend(("--model", str(params["model"])))
        if params.get("chain"):
            command.extend(("--chain", str(params["chain"])))
        return output / "STRUCTURE_EVIDENCE.json", self._execute(command, cwd=self.source_root, log_path=run_dir / "logs" / "structqc.log",
                                                                  cancel_event=cancel_event, on_process=on_process)

    def _write_memorient_manifest(self, run_dir: Path, inputs: Mapping[str, list[dict[str, Any]]], parameters: Mapping[str, Any]) -> None:
        output = run_dir / "outputs" / "memorient"
        alpha_scope = str(parameters.get("context", "")) in {"eukaryotic_pm", "tm_receptor"}
        document = {
            "schema_version": SCHEMA_VERSION, "module_id": "membrane_orientation", "version": "0.3.0",
            "input_sha256": {
                item["role"]: item["sha256"]
                for values in inputs.values() for item in values
                if item["role"] in {"structure", "topology_evidence"}
            },
            "parameters": {"context": parameters.get("context", "gram_negative_om"), "chain": parameters.get("chain", "")},
            "scientific_scope": {
                "scope_id": "alpha_helical" if alpha_scope else "beta_barrel",
                "scientific_state": "prototype" if alpha_scope else "conditionally_qualified",
                "release_blocking": not alpha_scope,
            },
            "runtime_versions": {"python": platform.python_version()},
            "outputs": sorted(p.name for p in output.iterdir() if p.is_file() and p.name != "RUN_MANIFEST.json"),
            "missing_evidence": [
                "independent_second_machine_reproduction",
                *(["external_topology_sidedness"] if alpha_scope else []),
            ],
            "limitations": [
                "Modeled orientation does not prove native intact-cell exposure.",
                *(["Alpha-helical orientation is experimental and is not part of the Mark 1 qualified scope."] if alpha_scope else []),
            ],
        }
        _write_json(output / "RUN_MANIFEST.json", document)

    def _generated_references(
        self, inputs: Mapping[str, list[dict[str, Any]]], params: Mapping[str, Any], run_dir: Path,
        alignment_map: Mapping[str, Any], alignment_map_sha256: str,
    ) -> Path:
        generated = run_dir / "generated"; generated.mkdir(parents=True, exist_ok=True)
        metadata = alignment_map["reference_metadata"]
        rows = []
        for state, role, evidence_key in (("active", "active_reference", "active_state_evidence"), ("inactive", "inactive_reference", "inactive_state_evidence")):
            for index, item in enumerate(inputs.get(role, []), 1):
                reference_id = f"{state.upper()}_{index}"
                reference_meta = metadata[reference_id]
                source = Path(item["materialized_path"]) if item.get("materialized_path") else self.object_path(item["sha256"])
                suffix = Path(item["file_name"]).suffix.lower()
                target = generated / f"{state}-{index}{suffix}"
                shutil.copyfile(source, target)
                rows.append({
                    "reference_id": reference_id, "state": state, "structure": target.name,
                    "structure_sha256": _sha_file(target),
                    "chain": str(reference_meta["chain_id"]),
                    "pdb_entry_id": str(reference_meta["pdb_entry_id"]).upper(),
                    "provenance": {"class": "experimental", "method": str(params.get("reference_method", ""))},
                    "state_evidence": {
                        "basis": str(params.get(evidence_key, "")),
                        "citation": str(params.get(f"{state}_reference_citation", "")),
                    },
                })
        domain = alignment_map["domain"]
        document = {
            "schema_version": "2.0",
            "reference_set_id": "ui-generated-abl-reference-set-v2",
            "qualification_scope": "abl_family",
            "subject_family": "ABL1",
            "alignment_map_sha256": alignment_map_sha256,
            "alignment_mask": {
                "coordinate_system": "uniprot",
                "uniprot_start": int(domain["uniprot_start"]),
                "uniprot_end": int(domain["uniprot_end"]),
                "minimum_coverage": 0.9,
            },
            "decision_rules": {"max_rmsd_A": float(params.get("max_rmsd_A", 2.5)),
                               "min_margin_A": float(params.get("min_margin_A", 0.25))},
            "references": rows,
        }
        path = generated / "REFERENCE_SET.json"; _write_json(path, document); return path

    @staticmethod
    def _fasta_sequence(path: Path, accession: str) -> str:
        records: list[tuple[str, str]] = []
        header = ""; chunks: list[str] = []
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                if header: records.append((header, "".join(chunks).upper()))
                header, chunks = line[1:].strip(), []
            elif header:
                chunks.append(line.strip())
        if header: records.append((header, "".join(chunks).upper()))
        for record_header, sequence in records:
            first = record_header.split()[0]
            identifiers = {first, *first.split("|")}
            if accession in identifiers:
                if not sequence or re.search(r"[^A-Z*]", sequence):
                    raise AnalysisError("SF-CSA query FASTA sequence is empty or invalid")
                return sequence.replace("*", "")
        raise AnalysisError(f"SF-CSA accession is absent from query FASTA: {accession}")

    def _sf_manifests(self, manifest: Mapping[str, Any], inputs: Mapping[str, list[dict[str, Any]]], run_dir: Path) -> tuple[Path, Path]:
        params = manifest.get("parameters", {})
        accession = str(params.get("accession", ""))
        _pack, database_source = self._sf_pack(str(params.get("database_pack", "")))
        interpretation = self._sf_interpretation(inputs)
        structure = self._one(inputs, "query_structure")
        fasta = self._one(inputs, "query_fasta")
        source_proteome = self._one(inputs, "source_proteome")
        sequence = self._fasta_sequence(fasta, accession)
        config = run_dir / "generated" / "sf_csa"; config.mkdir(parents=True, exist_ok=True)
        input_root = structure.parent
        query_manifest = {
            "schema_version": 1,
            "path_base": Path(os.path.relpath(input_root, config)).as_posix(),
            "release_scope": manifest["question"],
            "queries": [{
                "accession": accession, "common_name": accession,
                "organism": str(params.get("organism", "")), "strain": "not_declared",
                "uniprot_accession": accession, "decision_status": "REVIEW_SELECTED",
                "mechanism_group": str(params.get("mechanism_group", "")),
                "protein_specific_boundary": str(params.get("protein_specific_boundary", "")),
                "fasta_path": fasta.name, "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "structure_path": structure.name, "structure_sha256": _sha_file(structure),
                "source_proteome_path": source_proteome.name,
                "structure_class": "user_supplied_coordinate_model", "chain": str(params.get("chain", "A")),
                "residue_mapping": "exact sequence-to-coordinate mapping required and verified by SF-CSA",
                "orientation_artifact": "",
            }],
        }
        try: database_manifest = json.loads(database_source.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc: raise AnalysisError("SF-CSA database manifest is unreadable") from exc
        original_base = (database_source.parent / str(database_manifest.get("path_base", "."))).resolve()
        database_manifest["path_base"] = Path(os.path.relpath(original_base, config)).as_posix()
        for key in ("mechanism_families", "contested_groups", "divergence_sets", "classification_vocabulary"):
            database_manifest[key] = interpretation[key]
        database_manifest["source_pack_database_manifest_sha256"] = _sha_file(database_source)
        query_path = config / "target_manifest.json"; database_path = config / "database_manifest.json"
        _write_json(query_path, query_manifest); _write_json(database_path, database_manifest)
        return query_path, database_path

    def _run_registered(self, manifest: Mapping[str, Any], run_dir: Path,
                        cancel_event: Any | None = None, on_process: Any | None = None) -> tuple[list[dict[str, Any]], int]:
        inputs = self._materialize_inputs(self._inputs_by_role(manifest), run_dir); params = manifest.get("parameters", {})
        tool = manifest["analysis_type"]; steps: list[dict[str, Any]] = []; codes: list[int] = []
        if tool == "sf_csa":
            qc_inputs = dict(inputs)
            qc_inputs["structure"] = inputs.get("query_structure", [])
            qc_inputs["reference_fasta"] = inputs.get("query_fasta", [])
            struct_manifest, qc_code = self._run_structqc(
                manifest, qc_inputs, run_dir, cancel_event=cancel_event, on_process=on_process,
            )
            steps.append({"module_id": "structure_quality", "exit_code": qc_code}); codes.append(qc_code)
            if qc_code == 2 or not struct_manifest.is_file():
                return steps, 2
            query_manifest, database_manifest = self._sf_manifests(manifest, inputs, run_dir)
            output = run_dir / "outputs" / "sf_csa"
            command = [sys.executable, "-m", "sf_csa.cli", "run", "--queries", str(query_manifest),
                       "--databases", str(database_manifest), "--output", str(output)]
            code = self._execute(command, cwd=self.source_root, log_path=run_dir / "logs" / "sf-csa.log",
                                 cancel_event=cancel_event, on_process=on_process)
            steps.append({"module_id": "sf_csa", "exit_code": code}); codes.append(code)
            return steps, 2 if 2 in codes else 1 if 1 in codes else 0
        struct_manifest, code = self._run_structqc(manifest, inputs, run_dir, cancel_event=cancel_event, on_process=on_process)
        steps.append({"module_id": "structure_quality", "exit_code": code}); codes.append(code)
        if code == 2 or not struct_manifest.is_file():
            return steps, 2
        structure = self._one(inputs, "structure")
        if tool == "structure_qc":
            return steps, code
        if tool == "membrane_orientation":
            output = run_dir / "outputs" / "memorient"; output.mkdir(parents=True, exist_ok=True)
            command = [sys.executable, "-m", "memorient.cli", "orient", str(structure), "--context", str(params.get("context", "gram_negative_om")),
                       "--out-json", str(output / "MEMBRANE_ORIENTATION.json"), "--out-pdb", str(output / "ORIENTED_STRUCTURE.pdb"),
                       "--out-viz", str(output / "MEMBRANE_LAYER.json"), "--max-rows", "0"]
            if params.get("chain"): command.extend(("--chain", str(params["chain"])))
            topology_path = self._one(inputs, "topology_evidence", required=False)
            if topology_path: command.extend(("--topology-evidence", str(topology_path)))
            code = self._execute(command, cwd=self.source_root, log_path=run_dir / "logs" / "memorient.log", cancel_event=cancel_event, on_process=on_process)
            log = (run_dir / "logs" / "memorient.log").read_text(encoding="utf-8")
            table = log.split("\n# summary", 1)[0].strip()
            (output / "RESIDUE_ORIENTATION.tsv").write_text(table + ("\n" if table else ""), encoding="utf-8")
            self._write_memorient_manifest(run_dir, inputs, params)
            steps.append({"module_id": "membrane_orientation", "exit_code": code}); codes.append(code)
        elif tool == "conformational_state":
            alignment_path, alignment_document, mapped_query_chain = self._state_alignment(inputs)
            alignment_digest = _sha_file(alignment_path)
            refs = self._generated_references(inputs, params, run_dir, alignment_document, alignment_digest)
            output = run_dir / "outputs" / "state_atlas"
            command = [sys.executable, "-m", "state_atlas.cli", "run", "--manifest", str(struct_manifest),
                       "--reference-set", str(refs), "--alignment-map", str(alignment_path),
                       "--cluster-cutoff-A", str(params.get("cluster_cutoff_A", 2.0)),
                       "--stride", str(params.get("stride", 1)), "--out", str(output)]
            trajectory = self._one(inputs, "trajectory", required=False)
            if trajectory:
                command.extend(("--topology", str(structure), "--trajectory", str(trajectory), "--pbc", str(params.get("pbc", "none"))))
            else:
                command.extend(("--structure", str(structure)))
            query_chain = str(params.get("chain", "") or mapped_query_chain)
            if query_chain != mapped_query_chain:
                raise AnalysisError("declared query chain conflicts with the exact state alignment map")
            command.extend(("--chain", query_chain))
            code = self._execute(command, cwd=self.source_root, log_path=run_dir / "logs" / "state-atlas.log", cancel_event=cancel_event, on_process=on_process)
            steps.append({"module_id": "conformational_state", "exit_code": code}); codes.append(code)
        elif tool == "functional_site_state":
            output = run_dir / "outputs" / "site_context"
            command = [sys.executable, "-m", "site_context.cli", "run", "--manifest", str(struct_manifest),
                       "--structure", str(structure), "--annotations", str(self._one(inputs, "site_annotations")), "--out", str(output)]
            component = self._one(inputs, "component_map", required=False)
            if component: command.extend(("--component-map", str(component)))
            for item in inputs.get("pocket_result", []): command.extend(("--pocket-result", str(Path(item["materialized_path"]))))
            code = self._execute(command, cwd=self.source_root, log_path=run_dir / "logs" / "site-context.log", cancel_event=cancel_event, on_process=on_process)
            steps.append({"module_id": "site_context", "exit_code": code}); codes.append(code)
            annotation = self._one(inputs, "uniprot_annotations", required=False)
            if annotation:
                act_output = run_dir / "outputs" / "activity_state"
                act_input = run_dir / "generated" / "actstate-input"; structures = act_input / "structures"
                structures.mkdir(parents=True, exist_ok=True); act_output.mkdir(parents=True, exist_ok=True)
                if inputs.get("uniprot_annotations", [{}])[0].get("file_name", "").lower().endswith(".csv"):
                    with annotation.open(encoding="utf-8", newline="") as source_handle:
                        rows = list(csv.reader(source_handle))
                    with (act_input / "annotations.tsv").open("w", encoding="utf-8", newline="") as target_handle:
                        csv.writer(target_handle, delimiter="\t", lineterminator="\n").writerows(rows)
                else:
                    shutil.copyfile(annotation, act_input / "annotations.tsv")
                fasta = self._one(inputs, "reference_fasta"); shutil.copyfile(fasta, act_input / "sequences.fasta")
                shutil.copyfile(structure, structures / f"{manifest['subject_id']}{structure.suffix or '.pdb'}")
                command = [sys.executable, "-m", "actstate.cli", "run", "--in", str(act_input), "--out", str(act_output)]
                acode = self._execute(command, cwd=self.source_root, log_path=run_dir / "logs" / "actstate.log", cancel_event=cancel_event, on_process=on_process)
                steps.append({"module_id": "activity_state", "exit_code": acode}); codes.append(acode)
        elif tool == "assembly_interface":
            output = run_dir / "outputs" / "assembly_context"
            command = [sys.executable, "-m", "assembly_context.cli", "run", "--manifest", str(struct_manifest),
                       "--isolated", str(structure), "--assembly", str(self._one(inputs, "assembly")),
                       "--subject-chain", str(params.get("subject_chain", "")), "--relationship", str(params.get("relationship", "exact_protein")), "--out", str(output)]
            if params.get("assembly_id"): command.extend(("--assembly-id", str(params["assembly_id"])))
            if params.get("expected_chains"): command.extend(("--expected-chains", str(params["expected_chains"])))
            code = self._execute(command, cwd=self.source_root, log_path=run_dir / "logs" / "assembly-context.log", cancel_event=cancel_event, on_process=on_process)
            steps.append({"module_id": "assembly_context", "exit_code": code}); codes.append(code)
        else:
            raise AnalysisError(f"no registered execution adapter for {tool}")
        return steps, 2 if 2 in codes else 1 if 1 in codes else 0

    def _source_digests(self, analysis_type: str) -> dict[str, str | None]:
        """Bind a run to the exact scientific source trees that can affect it."""
        readiness = {row["analysis_type"]: row for row in tool_readiness(self.workspace) if not row.get("labs")}
        tools = [analysis_type]
        if analysis_type != "structure_qc":
            tools.insert(0, "structure_qc")
        digests = {
            tool: readiness[tool].get("package_source_sha256")
            for tool in tools
        }
        digests["workbench_orchestrator"] = _tree_sha(Path(__file__).resolve().parent)
        return digests

    def run(self, analysis_id: str, *, cancel_event: Any | None = None, on_process: Any | None = None) -> dict[str, Any]:
        preflight = self.preflight(analysis_id)
        if not preflight["valid"]:
            raise AnalysisError("analysis preflight is blocked")
        manifest = self.load(analysis_id)
        source_digests = self._source_digests(manifest["analysis_type"])
        identity = _sha_bytes(_canonical({
            "analysis_type": manifest["analysis_type"],
            "question": manifest["question"],
            "subject_id": manifest["subject_id"],
            "inputs": manifest["inputs"],
            "parameters": manifest["parameters"],
            "scientific_source_sha256": source_digests,
        }))
        run_id = f"run-{identity[:16]}"; run_dir = self._case_dir(analysis_id) / "runs" / run_id
        if (run_dir / "ANALYSIS_RUN.json").is_file():
            return json.loads((run_dir / "ANALYSIS_RUN.json").read_text(encoding="utf-8"))
        run_dir.mkdir(parents=True, exist_ok=True)
        try:
            steps, exit_code = self._run_registered(manifest, run_dir, cancel_event=cancel_event, on_process=on_process)
            status = "completed" if exit_code == 0 else "scientifically_incomplete" if exit_code == 1 else "failed"
            error = ""
        except InterruptedError as exc:
            steps, exit_code, status, error = [], 130, "cancelled", str(exc)
        except AnalysisError as exc:
            steps, exit_code, status, error = [], 1, "blocked", str(exc)
        record = {
            "schema_version": SCHEMA_VERSION, "contract_id": "analysis_run_record",
            "analysis_id": analysis_id, "run_id": run_id, "analysis_type": manifest["analysis_type"],
            "input_revision_sha256": manifest["revision_sha256"], "status": status,
            "scientific_source_sha256": source_digests,
            "scientifically_incomplete": status in {"scientifically_incomplete", "blocked"},
            "exit_code": exit_code, "steps": steps, "error": error,
            "limitations": [self._definitions[manifest["analysis_type"]]["claim_ceiling"]],
        }
        _write_json(run_dir / "ANALYSIS_RUN.json", record)
        self._render_report(manifest, record, run_dir)
        manifest = self.load(analysis_id); manifest["runs"].append({"run_id": run_id, "path": f"runs/{run_id}/ANALYSIS_RUN.json", "status": status})
        manifest["latest_run_id"] = run_id; manifest["state"] = status; manifest["revision"] = int(manifest["revision"]) + 1
        self._commit(self._case_dir(analysis_id), manifest)
        return record

    def _read_outputs(self, run_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
        documents, tables = [], []
        for path in sorted((run_dir / "outputs").rglob("*")) if (run_dir / "outputs").is_dir() else []:
            if not path.is_file(): continue
            relative = path.relative_to(run_dir).as_posix()
            if path.suffix.lower() == ".json" and path.stat().st_size <= 10 * 1024 * 1024:
                try: value = json.loads(path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError): continue
                documents.append({"path": relative, "document": value})
            elif path.suffix.lower() in {".tsv", ".csv"} and path.stat().st_size <= 20 * 1024 * 1024:
                delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
                try:
                    with path.open(encoding="utf-8", newline="") as handle: rows = list(csv.DictReader(handle, delimiter=delimiter))[:1000]
                except (OSError, csv.Error): continue
                tables.append({"path": relative, "columns": list(rows[0]) if rows else [], "rows": rows, "truncated": len(rows) == 1000})
        return documents, tables

    @staticmethod
    def _summary_value(documents: Iterable[Mapping[str, Any]], key: str) -> Any:
        for item in documents:
            doc = item.get("document", {})
            if key in doc: return doc[key]
        return None

    def _render_report(self, manifest: Mapping[str, Any], run: Mapping[str, Any], run_dir: Path) -> None:
        documents, tables = self._read_outputs(run_dir)
        definition = self._definitions[manifest["analysis_type"]]
        inputs = [{k: value[k] for k in ("role", "file_name", "sha256", "bytes")} for value in manifest["inputs"]]
        missing = []
        for item in documents:
            value = item["document"]
            if isinstance(value.get("missing_evidence"), list): missing.extend(map(str, value["missing_evidence"]))
            if isinstance(value.get("warnings"), list): missing.extend(map(str, value["warnings"]))
        report = {
            "schema_version": SCHEMA_VERSION, "contract_id": "scientific_report_manifest",
            "platform_identity": {"platform_id": PLATFORM_ID, "display_name": PLATFORM_DISPLAY_NAME,
                                  "edition": "Mark 1", "scientific_suite_name": PLATFORM_SCIENTIFIC_SUITE},
            "analysis_id": manifest["analysis_id"], "run_id": run["run_id"], "analysis_type": manifest["analysis_type"],
            "title": definition["title"], "research_question": manifest["question"], "subject_id": manifest["subject_id"],
            "status": run["status"], "scientifically_incomplete": run["scientifically_incomplete"],
            "claim_ceiling": definition["claim_ceiling"], "parameters": manifest["parameters"],
            "inputs": inputs, "steps": run["steps"], "missing_evidence": sorted(set(missing)),
            "limitations": run["limitations"], "documents": documents, "tables": tables,
            "metric_definitions": metric_definitions(),
        }
        _write_json(run_dir / "REPORT_DATA.json", report)
        _write_json(run_dir / "RUN_MANIFEST.json", {
            "schema_version": SCHEMA_VERSION, "analysis_id": manifest["analysis_id"], "run_id": run["run_id"],
            "platform_id": PLATFORM_ID,
            "analysis_type": manifest["analysis_type"], "input_revision_sha256": run["input_revision_sha256"],
            "scientific_source_sha256": run["scientific_source_sha256"],
            "input_sha256": {f"{item['role']}:{index + 1}": item["sha256"] for index, item in enumerate(inputs)},
            "parameters": manifest["parameters"], "runtime_versions": {"python": platform.python_version()},
            "outputs": ["REPORT_DATA.json", "REPORT.html", "RAW_EVIDENCE.zip", "CHECKSUMS.json"],
            "missing_evidence": report["missing_evidence"], "version_control": "missing" if not (self.source_root / ".git").exists() else "available",
        })
        report_html = self._report_html(report)
        (run_dir / "REPORT.html").write_text(report_html, encoding="utf-8")
        checksums = {}
        for path in sorted(p for p in run_dir.rglob("*") if p.is_file() and p.name not in {"CHECKSUMS.json", "RAW_EVIDENCE.zip"}):
            checksums[path.relative_to(run_dir).as_posix()] = _sha_file(path)
        _write_json(run_dir / "CHECKSUMS.json", {"schema_version": SCHEMA_VERSION, "files": checksums})
        self._write_bundle(manifest, run_dir)

    def _report_html(self, report: Mapping[str, Any]) -> str:
        esc = lambda value: html.escape(str(value))
        input_rows = "".join(f"<tr><td>{esc(i['role'])}</td><td>{esc(i['file_name'])}</td><td><code>{esc(i['sha256'])}</code></td><td>{i['bytes']:,}</td></tr>" for i in report["inputs"])
        step_rows = "".join(f"<tr><td>{esc(i['module_id'])}</td><td>{esc(i['exit_code'])}</td></tr>" for i in report["steps"])
        gaps = "".join(f"<li>{esc(value)}</li>" for value in report["missing_evidence"]) or "<li>None recorded.</li>"
        table_blocks = []
        for table in report["tables"][:8]:
            head = "".join(f"<th>{esc(c)}</th>" for c in table["columns"])
            rows = "".join("<tr>" + "".join(f"<td>{esc(row.get(c, ''))}</td>" for c in table["columns"]) + "</tr>" for row in table["rows"][:50])
            table_blocks.append(f"<section><h2>{esc(table['path'])}</h2><div class=table-wrap><table><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table></div></section>")
        style = """@page{size:auto;margin:16mm}*{box-sizing:border-box}body{font:14px/1.5 Inter,Arial,sans-serif;color:#17231f;max-width:1080px;margin:0 auto;padding:32px}h1{font-size:30px;margin:0 0 4px}h2{font-size:18px;border-bottom:1px solid #ccd7d2;padding-bottom:6px;margin-top:28px}.eyebrow{color:#087965;font-weight:700;text-transform:uppercase;letter-spacing:.1em}.status{display:inline-block;padding:4px 10px;border-radius:99px;background:#e8f4f0;font-weight:700}.notice{border-left:4px solid #d18324;background:#fff8eb;padding:12px 16px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:16px}.card{border:1px solid #ccd7d2;border-radius:10px;padding:14px}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-size:12px}th,td{text-align:left;vertical-align:top;padding:7px;border-bottom:1px solid #dde5e1}code{font:11px ui-monospace,monospace;word-break:break-all}.print{position:fixed;right:20px;top:20px}@media print{body{padding:0}.print{display:none}.grid{display:block}.card{break-inside:avoid}.table-wrap{overflow:visible}table{font-size:9px}tr{break-inside:avoid}}"""
        return "<!doctype html><html lang=en><head><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'><title>" + esc(PLATFORM_DISPLAY_NAME) + " — " + esc(report["title"]) + " report</title><style>" + style + "</style></head><body><button class=print onclick='window.print()'>Print / Save as PDF</button><p class=eyebrow>" + esc(PLATFORM_DISPLAY_NAME) + " · structural evidence report</p><h1>" + esc(report["title"]) + "</h1><p>" + esc(report["research_question"]) + "</p><p class=status>" + esc(report["status"]) + "</p><div class=notice><strong>Claim ceiling:</strong> " + esc(report["claim_ceiling"]) + "</div><div class=grid><section class=card><h2>Analysis</h2><p><strong>Subject:</strong> " + esc(report["subject_id"]) + "</p><p><strong>Run:</strong> <code>" + esc(report["run_id"]) + "</code></p></section><section class=card><h2>Missing or limited evidence</h2><ul>" + gaps + "</ul></section></div><section><h2>Methods executed</h2><table><thead><tr><th>Module</th><th>Exit code</th></tr></thead><tbody>" + step_rows + "</tbody></table></section><section><h2>Exact inputs</h2><table><thead><tr><th>Role</th><th>File</th><th>SHA-256</th><th>Bytes</th></tr></thead><tbody>" + input_rows + "</tbody></table></section>" + "".join(table_blocks) + "<section><h2>Interpretation boundaries</h2><ul>" + "".join(f"<li>{esc(v)}</li>" for v in report["limitations"]) + "</ul><p>No combined protein, activity, docking, druggability, or design score was calculated.</p></section></body></html>"

    def _write_bundle(self, manifest: Mapping[str, Any], run_dir: Path) -> None:
        output = run_dir / "RAW_EVIDENCE.zip"
        entries: list[tuple[str, bytes]] = []
        for index, item in enumerate(manifest["inputs"], 1):
            safe_name = Path(item["file_name"]).name
            entries.append((f"inputs/{index:03d}-{item['role']}-{safe_name}", self.object_path(item["sha256"]).read_bytes()))
        for path in sorted(p for p in run_dir.rglob("*") if p.is_file() and p != output):
            entries.append((path.relative_to(run_dir).as_posix(), path.read_bytes()))
        with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            for name, content in sorted(entries):
                info = zipfile.ZipInfo(name, FIXED_ZIP_TIME); info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                archive.writestr(info, content)

    def snapshot(self, analysis_id: str) -> dict[str, Any]:
        manifest = self.load(analysis_id); preflight = None; run = None; report = None
        preflight_path = self._case_dir(analysis_id) / "PREFLIGHT.json"
        if preflight_path.is_file(): preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
        if manifest.get("latest_run_id"):
            run_dir = self._case_dir(analysis_id) / "runs" / manifest["latest_run_id"]
            if (run_dir / "ANALYSIS_RUN.json").is_file(): run = json.loads((run_dir / "ANALYSIS_RUN.json").read_text(encoding="utf-8"))
            if (run_dir / "REPORT_DATA.json").is_file(): report = json.loads((run_dir / "REPORT_DATA.json").read_text(encoding="utf-8"))
        return {"analysis": manifest, "definition": self._definitions[manifest["analysis_type"]], "preflight": preflight, "run": run, "report": report}

    def artifact_path(self, analysis_id: str, run_id: str, relative: str) -> Path:
        if not re.fullmatch(r"run-[0-9a-f]{16}", run_id): raise AnalysisError("invalid run id")
        if relative.startswith("/") or ".." in Path(relative).parts: raise AnalysisError("invalid artifact path")
        root = (self._case_dir(analysis_id) / "runs" / run_id).resolve(); path = (root / relative).resolve()
        try: path.relative_to(root)
        except ValueError as exc: raise AnalysisError("artifact path escapes its run") from exc
        if not path.is_file(): raise AnalysisError("artifact is unavailable")
        return path

    def input_path(self, analysis_id: str, digest: str) -> Path:
        """Return a content object only when it belongs to the requested analysis."""
        manifest = self.load(analysis_id)
        if digest not in {item.get("sha256") for item in manifest.get("inputs", [])}:
            raise AnalysisError("input checksum is not attached to this analysis")
        return self.object_path(digest)

    def export(self, analysis_id: str, out_dir: str | Path) -> dict[str, Any]:
        manifest = self.load(analysis_id); run_id = manifest.get("latest_run_id")
        if not run_id: raise AnalysisError("analysis has no completed or incomplete run to export")
        source = self._case_dir(analysis_id) / "runs" / run_id; target = Path(out_dir).resolve()
        if target.exists() and any(target.iterdir()): raise AnalysisError("export directory must be empty")
        target.mkdir(parents=True, exist_ok=True)
        for name in ("REPORT_DATA.json", "REPORT.html", "RAW_EVIDENCE.zip", "CHECKSUMS.json", "RUN_MANIFEST.json"):
            shutil.copyfile(source / name, target / name)
        return {"analysis_id": analysis_id, "run_id": run_id, "out": str(target), "files": sorted(p.name for p in target.iterdir())}
