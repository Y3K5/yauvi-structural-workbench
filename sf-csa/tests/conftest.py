"""Fixtures for the SF-CSA tests."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


SEQUENCE = "MKVLAAGIVGLTTHAADQPRSTWY"


@pytest.fixture
def campaign(tmp_path: Path) -> Path:
    """A complete, minimal campaign on disk, with no real organism in it.

    Everything the manifest builder touches is present: a sequence ledger, a
    decision ledger, a FASTA whose digest matches the ledger, a structure, and a
    stand-in structure database with a version file.
    """
    root = tmp_path / "campaign"
    (root / "results").mkdir(parents=True)
    (root / "proteomes").mkdir()
    (root / "db").mkdir()

    (root / "proteomes" / "example.faa").write_text(
        f">sp|EX0001|EXAMPLE_ORG\n{SEQUENCE}\n", encoding="utf-8"
    )
    (root / "results" / "SEQUENCE_MANIFEST.tsv").write_text(
        "accession\tsource_file\tsequence_sha256\n"
        f"EX0001\tproteomes/example.faa\t{_sha(SEQUENCE)}\n",
        encoding="utf-8",
    )
    (root / "results" / "SELECTION_LEDGER.tsv").write_text(
        "primary_accession\tdecision_status\nEX0001\tSELECTED\n", encoding="utf-8"
    )
    (root / "results" / "EX0001.pdb").write_text(
        "ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00 50.00           C\nEND\n",
        encoding="utf-8",
    )
    (root / "results" / "orient_manifest.json").write_text("{}\n", encoding="utf-8")

    database = root / "db" / "structdb"
    database.write_text("stand-in structure database\n", encoding="utf-8")
    (root / "db" / "pdb.version").write_text("2026-01-01\nfrozen snapshot\n", encoding="utf-8")

    spec = {
        "schema_version": 1,
        "release_scope": "one synthetic target",
        "root": ".",
        "sequence_manifest": "results/SEQUENCE_MANIFEST.tsv",
        "decision_ledger": "results/SELECTION_LEDGER.tsv",
        "default_orientation_artifact": "results/orient_manifest.json",
        "targets": [
            {
                "accession": "EX0001",
                "common_name": "ExampleProtein",
                "organism": "Examplus fictus",
                "strain": "T1",
                "uniprot_accession": "EX0001",
                "mechanism_group": "example_family",
                "protein_specific_boundary": "architecture does not establish substrate.",
                "structure_path": "results/EX0001.pdb",
                "source_proteome_path": "proteomes/example.faa",
            }
        ],
        "database": {
            "pdb_database": "db/structdb",
            "proteome_globs": ["proteomes/*.faa"],
            "thresholds": {"structure_evalue": "0.01", "same_fold_tm": 0.5},
        },
    }
    spec_path = root / "campaign.json"
    spec_path.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
    return spec_path


@pytest.fixture
def spec_document(campaign: Path) -> dict:
    return json.loads(campaign.read_text(encoding="utf-8"))


def write_spec(path: Path, document: dict) -> Path:
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    return path
