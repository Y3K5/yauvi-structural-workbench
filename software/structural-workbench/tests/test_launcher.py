"""Launcher guards: refuse to attach to a server that is not the one intended.

Each check exists so the launcher never silently hands a user a page served by a
different build, a different analysis library, or a server running with
different reference-retrieval settings than they asked for. A wrong attachment
would show someone another workspace's evidence and look perfectly normal.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from yauvi_platform.structural_workbench import AnalysisError
from yauvi_structural_workbench.build_info import workspace_identity
from yauvi_structural_workbench.launcher import validate_running


def running_state(workspace: Path, **overrides) -> dict:
    state = {
        "application_id": "yauvi-structural-workbench",
        "version": "0.1.0.dev0",
        "build_id": "b" * 64,
        "interface_id": "i" * 64,
        "workspace_id": workspace_identity(workspace),
        "restart_required": False,
        "reference_fetch_enabled": False,
        "workspace_label": "Local workspace",
    }
    state.update(overrides)
    return state


@pytest.fixture
def expected() -> dict:
    return {
        "application_id": "yauvi-structural-workbench",
        "version": "0.1.0.dev0",
        "build_id": "b" * 64,
        "interface_id": "i" * 64,
    }


def test_a_matching_server_is_accepted(tmp_path, expected):
    validate_running(running_state(tmp_path), expected, tmp_path, False)


@pytest.mark.parametrize("field", ["application_id", "version", "build_id", "interface_id"])
def test_a_different_build_is_refused(tmp_path, expected, field):
    current = running_state(tmp_path, **{field: "something-else"})

    with pytest.raises(AnalysisError, match="different build is already running"):
        validate_running(current, expected, tmp_path, False)


def test_a_missing_identity_field_is_refused(tmp_path, expected):
    """An older server that reports nothing for a field must not be treated as a match."""
    current = running_state(tmp_path)
    del current["build_id"]

    with pytest.raises(AnalysisError, match="different build is already running"):
        validate_running(current, expected, tmp_path, False)


def test_a_server_for_another_library_is_refused(tmp_path, expected):
    other = tmp_path / "another-library"
    other.mkdir()
    current = running_state(other)

    with pytest.raises(AnalysisError, match="another analysis library"):
        validate_running(current, expected, tmp_path, False)


def test_a_server_whose_files_changed_is_refused(tmp_path, expected):
    current = running_state(tmp_path, restart_required=True)

    with pytest.raises(AnalysisError, match="files changed after this server started"):
        validate_running(current, expected, tmp_path, False)


def test_different_reference_retrieval_settings_are_refused(tmp_path, expected):
    """Attaching to a fetch-enabled server when fetching was not asked for."""
    current = running_state(tmp_path, reference_fetch_enabled=True)

    with pytest.raises(AnalysisError, match="reference-retrieval settings"):
        validate_running(current, expected, tmp_path, False)


def test_a_server_wanting_fetch_is_refused_when_running_without_it(tmp_path, expected):
    current = running_state(tmp_path, reference_fetch_enabled=False)

    with pytest.raises(AnalysisError, match="reference-retrieval settings"):
        validate_running(current, expected, tmp_path, True)


def test_a_differently_labelled_workspace_is_refused(tmp_path, expected):
    current = running_state(tmp_path, workspace_label="Someone else's workspace")

    with pytest.raises(AnalysisError, match="another workspace label"):
        validate_running(current, expected, tmp_path, False)


def test_surrounding_whitespace_in_a_label_is_not_a_mismatch(tmp_path, expected):
    current = running_state(tmp_path, workspace_label="Study one")

    validate_running(current, expected, tmp_path, False, label="  Study one  ")
