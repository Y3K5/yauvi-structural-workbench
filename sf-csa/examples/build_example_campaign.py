#!/usr/bin/env python3
"""Build a complete, runnable example campaign in a directory of your choosing.

Everything it writes is synthetic: a fictional organism, one target, a stand-in
structure database. The point is to show the shape a campaign spec must have and
to give `build-manifests` and `validate` something real to work on, with no
download and no licence attached.

    python3 examples/build_example_campaign.py /tmp/demo
    sf-csa build-manifests --spec /tmp/demo/campaign.json --out /tmp/demo/config
    sf-csa validate --queries /tmp/demo/config/target_manifest.json \
                    --databases /tmp/demo/config/database_manifest.json

`sf-csa run` is deliberately NOT part of this demo: it needs Foldseek and DIAMOND
on PATH and a real structure database, and a demo that pretended otherwise would
be demonstrating the wrong thing.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

SEQUENCE = "MKVLAAGIVGLTTHAADQPRSTWYCNDEQKRHILMFPSTWYV"


def build(root: Path) -> Path:
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "proteomes").mkdir(exist_ok=True)
    (root / "db").mkdir(exist_ok=True)

    (root / "proteomes" / "examplus.faa").write_text(
        f">sp|EX0001|EXPR_EXAMP Example protein\n{SEQUENCE}\n", encoding="utf-8"
    )

    # The sequence ledger is the authority on what an accession is. build-manifests
    # re-reads the FASTA and refuses if the digest has moved.
    (root / "results" / "SEQUENCE_MANIFEST.tsv").write_text(
        "accession\tsource_file\tsequence_sha256\n"
        f"EX0001\tproteomes/examplus.faa\t"
        f"{hashlib.sha256(SEQUENCE.encode('utf-8')).hexdigest()}\n",
        encoding="utf-8",
    )
    (root / "results" / "SELECTION_LEDGER.tsv").write_text(
        "primary_accession\tdecision_status\nEX0001\tSELECTED\n", encoding="utf-8"
    )
    (root / "results" / "EX0001.pdb").write_text(
        "HEADER    EXAMPLE                                 01-JAN-26   0EXA\n"
        "ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00 50.00           C\n"
        "ATOM      2  CA  LYS A   2       3.800   0.000   0.000  1.00 50.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    (root / "results" / "orient_manifest.json").write_text("{}\n", encoding="utf-8")

    # A stand-in for the frozen structure database. Its version file is not
    # optional: an unversioned database cannot be cited, so the builder refuses one.
    (root / "db" / "structdb").write_text("stand-in structure database\n", encoding="utf-8")
    (root / "db" / "pdb.version").write_text(
        "2026-01-01\nsynthetic example snapshot\n", encoding="utf-8"
    )

    spec = {
        "schema_version": 1,
        "release_scope": "one synthetic target, for demonstrating the spec shape",
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
                "protein_specific_boundary":
                    "Shared architecture does not establish a shared substrate.",
                "structure_path": "results/EX0001.pdb",
                "source_proteome_path": "proteomes/examplus.faa",
            }
        ],
        "database": {
            "pdb_database": "db/structdb",
            "proteome_globs": ["proteomes/*.faa"],
            "thresholds": {
                "structure_evalue": "0.01",
                "same_fold_tm": 0.5,
                "whole_architecture_coverage": 0.7,
                "max_structure_hits": 50,
                "sequence_evalue": "1e-5",
                "sequence_min_identity": 20,
                "sequence_min_query_coverage": 40,
                "sequence_max_hits": 5000,
                "sequence_hits_per_proteome": 10,
            },
            # A campaign that is not about periodontal pathogens should say so
            # here rather than inherit the built-in defaults.
            "mechanism_families": [
                {"group": "example_family", "pattern": r"example|demonstration"}
            ],
            "contested_groups": [],
            "divergence_sets": [],
            "coverage_boundary":
                "search is exhaustive only within the supplied database; missing models "
                "are reported, not treated as negatives",
        },
    }
    spec_path = root / "campaign.json"
    spec_path.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
    return spec_path


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "sf_csa_example").resolve()
    spec = build(root)
    print(f"wrote an example campaign under {root}")
    print(f"  sf-csa build-manifests --spec {spec} --out {root / 'config'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
