"""Build checksum-pinned SF-CSA input manifests from a campaign specification.

The previous version of this file *was* the campaign: eight accessions, their
structure paths, their organisms, strains and mechanism groups, and the exact
result directories of one pipeline run, all as module-level dictionaries. The
descriptor said so plainly — *"organism specificity still lives in
build_manifests.py and must be rewritten per organism"* — which meant running
SF-CSA on anything else was a source edit.

Now the campaign is data. A spec file names the targets and where their inputs
live; this module resolves, checksums and validates them. Nothing here knows
what organism it is looking at.

The split also matters for what may be published: the target dictionaries are
private campaign material, and keeping them out of source is the same severance
that `redvax.construct.pg_tf_panel` -> `panel` performed for construct assembly.

Spec shape (JSON or YAML-free JSON; see `examples/campaign_spec.json`):

    {
      "schema_version": 1,
      "release_scope": "free text, recorded in the manifest",
      "path_base": "../..",
      "root": ".",                      # resolved relative to the spec file
      "sequence_manifest": "results/.../SEQUENCE_MANIFEST.tsv",
      "decision_ledger":   "results/.../SELECTION_LEDGER.tsv",
      "default_orientation_artifact": "results/.../orient_manifest.json",
      "targets": [
        {
          "accession": "...",
          "common_name": "...", "organism": "...", "strain": "...",
          "uniprot_accession": "...", "mechanism_group": "...",
          "protein_specific_boundary": "...",
          "structure_path": "results/.../x_oriented.pdb",
          "source_proteome_path": "proteomes/.../y.faa",
          "orientation_artifact": "optional override"
        }
      ],
      "database": { ... database-manifest fields ... }
    }
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Mapping

from . import core
from .core import select_fasta, sequence_sha, sha256

REQUIRED_TARGET_FIELDS = (
    "accession",
    "common_name",
    "organism",
    "strain",
    "uniprot_accession",
    "mechanism_group",
    "protein_specific_boundary",
    "structure_path",
    "source_proteome_path",
)


class SpecError(RuntimeError):
    """The campaign spec is malformed, or refers to something that is not there."""


def load_spec(path: str | Path) -> tuple[dict[str, Any], Path]:
    """Read a spec and resolve its root. Returns (spec, root)."""
    spec_path = Path(path).resolve()
    if not spec_path.is_file():
        raise SpecError(f"campaign spec not found: {spec_path}")
    try:
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
    except ValueError as exc:
        raise SpecError(f"campaign spec is not valid JSON: {spec_path}: {exc}") from exc
    if not isinstance(spec, Mapping):
        raise SpecError(f"campaign spec must be a JSON object: {spec_path}")

    targets = spec.get("targets")
    if not isinstance(targets, list) or not targets:
        raise SpecError(f"campaign spec declares no targets: {spec_path}")

    seen: set[str] = set()
    for index, target in enumerate(targets):
        if not isinstance(target, Mapping):
            raise SpecError(f"target #{index} is not an object")
        missing = [f for f in REQUIRED_TARGET_FIELDS if not target.get(f)]
        if missing:
            label = target.get("accession") or f"#{index}"
            raise SpecError(f"target {label} is missing required field(s): {', '.join(missing)}")
        accession = str(target["accession"])
        if accession in seen:
            raise SpecError(f"target {accession} is declared twice")
        seen.add(accession)

    root = (spec_path.parent / str(spec.get("root", "."))).resolve()
    if not root.is_dir():
        raise SpecError(f"campaign root does not exist: {root}")
    return dict(spec), root


def _read_tsv_index(path: Path, key: str) -> dict[str, dict[str, str]]:
    if not path.is_file():
        raise SpecError(f"table not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or key not in reader.fieldnames:
            raise SpecError(f"{path} has no {key!r} column")
        return {row[key]: row for row in reader if row.get(key)}


def build_target_manifest(spec: Mapping[str, Any], root: Path) -> dict[str, Any]:
    """Resolve every target to checksummed paths, verifying the sequence ledger."""
    sequences = _read_tsv_index(root / spec["sequence_manifest"], "accession")
    decisions = _read_tsv_index(root / spec["decision_ledger"], "primary_accession")
    default_orientation = spec.get("default_orientation_artifact", "")

    queries: list[dict[str, Any]] = []
    for target in spec["targets"]:
        accession = str(target["accession"])

        ledger_row = sequences.get(accession)
        if ledger_row is None:
            raise SpecError(f"{accession} is not in the sequence manifest")
        if "primary_accession" in next(iter(decisions.values()), {}) or decisions:
            if accession not in decisions:
                raise SpecError(f"{accession} is not in the decision ledger")

        fasta = root / ledger_row["source_file"]
        record = select_fasta(fasta, accession)
        # The ledger is the authority on what sequence this accession is. A
        # mismatch means the FASTA moved underneath the campaign, and continuing
        # would checksum the wrong protein.
        if sequence_sha(record["sequence"]) != ledger_row["sequence_sha256"]:
            raise SpecError(f"sequence ledger drift: {accession}")

        structure = root / str(target["structure_path"])
        if not structure.is_file():
            raise SpecError(f"{accession}: structure not found: {structure}")

        queries.append(
            {
                "accession": accession,
                "common_name": target["common_name"],
                "organism": target["organism"],
                "strain": target["strain"],
                "uniprot_accession": target["uniprot_accession"],
                "decision_status": decisions[accession]["decision_status"],
                "mechanism_group": target["mechanism_group"],
                "protein_specific_boundary": target["protein_specific_boundary"],
                "fasta_path": fasta.relative_to(root).as_posix(),
                "sequence_sha256": ledger_row["sequence_sha256"],
                "structure_path": structure.relative_to(root).as_posix(),
                "source_proteome_path": target["source_proteome_path"],
                "structure_sha256": sha256(structure),
                "structure_class": target.get(
                    "structure_class", "exact_predicted_monomer_membrane_oriented"
                ),
                "chain": target.get("chain", "A"),
                "residue_mapping": target.get(
                    "residue_mapping",
                    "exact full-length 1-based sequence mapping; verified at run time",
                ),
                "orientation_artifact": target.get(
                    "orientation_artifact", default_orientation
                ),
            }
        )

    return {
        "schema_version": 1,
        "path_base": spec.get("path_base", "../.."),
        "release_scope": spec.get("release_scope", ""),
        "queries": queries,
    }


def build_database_manifest(
    spec: Mapping[str, Any], root: Path, *, query_count: int
) -> dict[str, Any]:
    """Assemble the database manifest, checksumming the structure database."""
    database = dict(spec.get("database") or {})
    if not database:
        raise SpecError("campaign spec declares no `database` block")

    pdb_relative = database.get("pdb_database")
    if not pdb_relative:
        raise SpecError("`database.pdb_database` is required")
    pdb = (root / str(pdb_relative)).resolve()
    if not pdb.exists():
        raise SpecError(f"structure database not found: {pdb}")

    version_file = pdb.parent / "pdb.version"
    version = (
        version_file.read_text(encoding="utf-8").strip().replace("\n", "; ")
        if version_file.is_file()
        else database.get("pdb_database_version", "")
    )
    if not version:
        raise SpecError(
            f"no version recorded for the structure database. Expected {version_file}, "
            f"or `database.pdb_database_version` in the spec. An unversioned database "
            f"cannot be cited."
        )

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "path_base": spec.get("path_base", "../.."),
        **database,
        "pdb_database_version": version,
        "pdb_database_checksum": core.database_bundle_checksum(pdb),
        "pdb_database_file_checksums": core.database_bundle_file_checksums(pdb),
    }
    # `core.run_pipeline` requires this key. Defaulting it to the closed
    # vocabulary means a spec need not restate it, while a spec that does state
    # one still wins — a release built against a different vocabulary must say so.
    manifest.setdefault("classification_vocabulary", list(core.CLASSIFICATION_VOCABULARY))

    # These tables were literals in core.py. They live in the manifest so that a
    # campaign against other organisms is a manifest change, and so `verify`
    # audits a release against something the release did not write itself.
    manifest.setdefault("mechanism_families", core.DEFAULT_MECHANISM_FAMILIES)
    manifest.setdefault("contested_groups", core.DEFAULT_CONTESTED_GROUPS)
    manifest.setdefault("divergence_sets", core.DEFAULT_DIVERGENCE_SETS)

    expectations = dict(manifest.get("release_expectations") or {})
    expectations["query_count"] = query_count
    expectations.setdefault("title_traps", core.DEFAULT_TITLE_TRAPS)
    manifest["release_expectations"] = expectations
    return manifest


def build(spec_path: str | Path, out_dir: str | Path) -> tuple[Path, Path]:
    """Build both manifests from a spec. Returns (target_manifest, database_manifest)."""
    spec, root = load_spec(spec_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    targets = build_target_manifest(spec, root)
    databases = build_database_manifest(spec, root, query_count=len(targets["queries"]))

    target_path = out / "target_manifest.json"
    database_path = out / "database_manifest.json"
    target_path.write_text(json.dumps(targets, indent=2) + "\n", encoding="utf-8")
    database_path.write_text(json.dumps(databases, indent=2) + "\n", encoding="utf-8")
    return target_path, database_path
