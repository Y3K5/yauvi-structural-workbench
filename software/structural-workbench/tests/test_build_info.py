"""Build and workspace identities: stable when nothing changed, different when it did.

These identities are what the launcher compares to decide whether a running
server is the one the user meant. If a digest were unstable across runs the
launcher would refuse valid servers; if it ignored a changed file it would
accept the wrong one.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from yauvi_structural_workbench import build_info
from yauvi_structural_workbench.build_info import (
    application_build,
    digest,
    workspace_identity,
)


def test_a_digest_does_not_depend_on_key_order(tmp_path):
    """Dict ordering is an implementation detail and must not change identity."""
    assert digest({"a": 1, "b": 2}) == digest({"b": 2, "a": 1})


def test_a_digest_is_stable_across_calls():
    assert digest({"file": "abc"}) == digest({"file": "abc"})


def test_a_changed_value_changes_the_digest():
    assert digest({"file": "abc"}) != digest({"file": "abd"})


def test_an_added_entry_changes_the_digest():
    assert digest({"a": 1}) != digest({"a": 1, "b": 2})


def test_workspace_identity_ignores_how_the_path_was_spelled(tmp_path):
    """The same folder reached two ways is the same workspace."""
    nested = tmp_path / "library" / "inner"
    nested.mkdir(parents=True)
    indirect = tmp_path / "library" / "inner" / ".." / "inner"

    assert workspace_identity(nested) == workspace_identity(indirect)
    assert workspace_identity(nested) == workspace_identity(str(nested))


def test_workspace_identity_separates_different_folders(tmp_path):
    first = tmp_path / "one"
    second = tmp_path / "two"
    first.mkdir()
    second.mkdir()

    assert workspace_identity(first) != workspace_identity(second)


def test_a_missing_application_package_is_reported_by_name(monkeypatch):
    """A build that cannot see one of its own packages must fail, not under-report."""
    monkeypatch.setattr(build_info, "PACKAGES",
                        build_info.PACKAGES + ("yauvi_absent_package",))

    with pytest.raises(RuntimeError, match="yauvi_absent_package"):
        application_build(assets={"index.html": b"<html></html>"})


def test_a_changed_browser_asset_changes_both_identities():
    """The interface identity tracks the pages; the build identity covers everything."""
    first = application_build(assets={"index.html": b"<html>one</html>"})
    second = application_build(assets={"index.html": b"<html>two</html>"})

    assert first["interface_id"] != second["interface_id"]
    assert first["build_id"] != second["build_id"]


def test_identical_assets_give_identical_identities():
    first = application_build(assets={"index.html": b"<html>same</html>"})
    second = application_build(assets={"index.html": b"<html>same</html>"})

    assert first["interface_id"] == second["interface_id"]
    assert first["build_id"] == second["build_id"]


def test_the_build_record_declares_what_it_does_not_cover():
    """The record must not be mistaken for scientific qualification."""
    record = application_build(assets={"index.html": b"<html></html>"})

    assert record["application_id"] == "yauvi-structural-workbench"
    assert record["channel"] == "Development preview"
    assert "qualification" in record["identity_scope"]
    assert record["distribution"] in {"Source checkout", "Installed package"}


def test_the_build_identity_is_not_the_interface_identity():
    """Two separate questions: did the pages change, did anything change."""
    record = application_build(assets={"index.html": b"<html></html>"})

    assert record["build_id"] != record["interface_id"]
