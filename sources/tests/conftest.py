"""Shared fixtures.

Two kinds of test live here. Most run against a small synthetic registry built
in-process, so they are fast and state exactly what they depend on. A few run
against the real `catalogs/sources.yaml`, because the thing most worth knowing
about this package is whether it can still read the file it exists to read.
"""
from __future__ import annotations

from pathlib import Path
import textwrap

import pytest

from yauvi_sources import SourceCache, SourceRegistry
from yauvi_sources.manifest import parse_manifest


def _workspace_root() -> Path | None:
    """The Protein Platform root, if these tests are running inside it."""
    here = Path(__file__).resolve()
    for directory in here.parents:
        if (directory / "catalogs" / "sources.yaml").is_file():
            return directory
    return None


@pytest.fixture(scope="session")
def workspace() -> Path:
    root = _workspace_root()
    if root is None:
        pytest.skip("not running inside the workspace; real-registry tests need catalogs/")
    return root


@pytest.fixture(scope="session")
def real_registry(workspace: Path) -> SourceRegistry:
    return SourceRegistry.load(workspace / "catalogs" / "sources.yaml")


SYNTHETIC_REGISTRY = textwrap.dedent(
    """
    schema_version: "1.0"
    catalog_id: test-registry
    updated_at: "2026-08-19"
    sources:
      - source_id: open_db
        display_name: An Open Database
        kind: sequence_db
        channel: conservation
        status: wired
        access: network_download
        url: https://example.invalid/open_db.fasta
        license_note: CC BY 4.0.
      - source_id: gated_db
        display_name: A Licence-Gated Database
        kind: sequence_db
        channel: essentiality
        status: configured_optional
        access: manual_download
        license_note: Registration required.
      - source_id: a_table
        display_name: A Predictor We Never Run
        kind: predictor
        channel: antigenicity
        status: table_only
        access: web_server_manual
      - source_id: a_binary
        display_name: A Command-Line Tool
        kind: aligner
        channel: infrastructure
        status: wired
        access: local_binary
      - source_id: a_heuristic
        display_name: An Internal Heuristic
        kind: internal_heuristic
        channel: localization
        status: wired
        access: internal
      - source_id: reachable_over_api
        display_name: Reachable But Deliberately Not Automated
        kind: annotation_db
        channel: antigenicity
        status: table_only
        access: network_api_and_download
    """
).strip()


@pytest.fixture
def registry_path(tmp_path: Path) -> Path:
    path = tmp_path / "sources.yaml"
    path.write_text(SYNTHETIC_REGISTRY, encoding="utf-8")
    return path


@pytest.fixture
def registry(registry_path: Path) -> SourceRegistry:
    return SourceRegistry.load(registry_path)


@pytest.fixture
def cache(tmp_path: Path) -> SourceCache:
    return SourceCache(tmp_path / "cache")


@pytest.fixture
def manifest():
    return parse_manifest(
        {
            "schema_version": "1.0",
            "module_id": "testmod",
            "requires": [
                {"source_id": "open_db", "role": "the main input", "required": True},
                {"source_id": "gated_db", "role": "optional reference", "required": False},
                {"source_id": "a_table", "role": "scoring table", "required": False},
                {"source_id": "a_binary", "role": "aligner", "required": True},
                {"source_id": "a_heuristic", "role": "fallback", "required": False},
            ],
        },
        origin="<fixture>",
    )
