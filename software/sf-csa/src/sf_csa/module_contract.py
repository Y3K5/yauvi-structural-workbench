"""The module's declared interface, and the checks behind it.

`describe()` output is what `catalogs/modules/sf_csa.yaml` is checked against by
`tools/tests/test_module_contract.py`. A module that quietly stops emitting what
it declares should fail a test, not surprise someone reading the results.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

MODULE_ID = "sf_csa"
SOURCES_MANIFEST = Path(__file__).with_name("sources.yaml")


def describe() -> dict[str, Any]:
    from . import __version__
    from .core import CLASSIFICATION_VOCABULARY

    return {
        # Two names, deliberately both reported. `module_id` is how the platform
        # addresses this module; `package` is how pip and `yauvi-fetch --for`
        # address it.
        "module_id": MODULE_ID,
        "package": "sf_csa",
        "display_name": "Structure-Function Comparative Species Analysis",
        "version": __version__,
        "one_line": (
            "How do this protein's structure and function compare across the species "
            "we target?"
        ),
        "commands": ["run", "verify", "build-manifests", "describe", "validate", "fetch"],
        "inputs": [
            {"name": "--queries", "format": "json", "required": True,
             "note": "target manifest; built from a campaign spec by `build-manifests`"},
            {"name": "--databases", "format": "json", "required": True,
             "note": "database manifest, carrying thresholds and the pinned database checksum"},
        ],
        "outputs": [
            {"name": "SF_CSA_RELEASE_MANIFEST.json", "contract": "sf_csa_release_manifest",
             "format": "json"},
            {"name": "RELEASE_COMPARISON_MATRIX.tsv", "contract": "sf_csa_comparison_matrix",
             "format": "tsv"},
            {"name": "CHECKSUMS.json", "contract": "sf_csa_checksums", "format": "json"},
        ],
        "classification_vocabulary": list(CLASSIFICATION_VOCABULARY),
        "runtimes": ["yauvi-python", "foldseek", "diamond"],
        "optional_runtimes": [],
        "reproducible": True,
        "limitations": [
            "Structural and sequence similarity are reported separately and must not be "
            "merged into a single similarity claim.",
            "Query structures are exact predicted monomers, not experimental assemblies "
            "or active poses. Proteins without a local structure are emitted as "
            "candidate_missing_structure, never as structural negatives.",
            "Structural-category thresholds are set in the database manifest, so two "
            "releases built against different manifests are not comparable on category "
            "alone.",
            "The interpretation vocabulary is closed. A hit that does not fit one of the "
            "six labels is unresolved_or_conflicted, not a new category.",
            "Mechanism families, contested groups and divergence sets are manifest "
            "entries. Their default values are periodontal-pathogen biology; a campaign "
            "against other organisms that does not override them is using the wrong "
            "table, and the release records which table it used.",
            "The module verifies and records an existing release. It does not launch a "
            "structure sweep implicitly; with no release present it blocks.",
        ],
    }


def validate_manifests(queries: str | Path, databases: str | Path) -> list[tuple[str, bool, str]]:
    """Check both manifests are present, parseable, and internally consistent.

    Returns (check_name, ok, detail) triples; the caller decides how to report.
    """
    checks: list[tuple[str, bool, str]] = []

    documents: dict[str, Any] = {}
    for name, path in (("queries", Path(queries)), ("databases", Path(databases))):
        if not path.is_file():
            checks.append((f"manifest:{name}", False, f"no such file: {path}"))
            continue
        try:
            documents[name] = json.loads(path.read_text(encoding="utf-8"))
        except ValueError as exc:
            checks.append((f"manifest:{name}", False, f"{path.name}: {exc}"))
            continue
        checks.append((f"manifest:{name}", True, str(path)))

    target_manifest = documents.get("queries")
    if isinstance(target_manifest, dict):
        queries_list = target_manifest.get("queries") or []
        checks.append(
            ("queries:present", bool(queries_list), f"{len(queries_list)} query/queries")
            if queries_list
            else ("queries:present", False, "target manifest declares no queries")
        )
        missing_groups = [
            q.get("accession", "?") for q in queries_list if not q.get("mechanism_group")
        ]
        checks.append(
            ("queries:mechanism_group", not missing_groups,
             "every query names a mechanism group" if not missing_groups
             else f"no mechanism_group for: {missing_groups}")
        )

    database_manifest = documents.get("databases")
    if isinstance(database_manifest, dict):
        for key in ("pdb_database", "pdb_database_checksum", "thresholds",
                    "classification_vocabulary"):
            checks.append(
                (f"databases:{key}", key in database_manifest,
                 "present" if key in database_manifest else "missing from the database manifest")
            )
        # A campaign that has not overridden the default tables is using
        # periodontal-pathogen biology. That is legitimate for that campaign and
        # wrong for any other, so it is surfaced rather than assumed.
        for key in ("mechanism_families", "contested_groups", "divergence_sets"):
            checks.append(
                (f"databases:{key}", True,
                 "declared in the manifest" if key in database_manifest
                 else "not declared — the built-in periodontal-pathogen default will be used")
            )

    return checks
