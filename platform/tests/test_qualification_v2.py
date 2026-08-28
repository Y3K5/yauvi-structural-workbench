from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
QUALIFICATION = ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2"


def runner_module():
    spec = importlib.util.spec_from_file_location("qualification_v2_runner", QUALIFICATION / "run_qualification.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_v2_panel_freezes_all_scopes_strata_and_gates():
    manifest = json.loads((QUALIFICATION / "PANEL_MANIFEST.json").read_text(encoding="utf-8"))
    assert manifest["manifest_state"] == "source_adoption_required"
    assert len(manifest["release_blocking_scopes"]) == 6
    assert "membrane_orientation:alpha_helical" in manifest["non_blocking_scopes"]
    panels = {row["workflow"]: row for row in manifest["panels"]}
    assert set(panels) == {
        "structure_qc", "membrane_orientation", "conformational_state",
        "functional_site_state", "assembly_interface", "sf_csa",
    }
    assert sum(row["count"] for row in panels["membrane_orientation"]["requirements"]) == 32
    assert panels["membrane_orientation"]["gates"]["empty_extracellular_comparison"] == "not_applicable"
    assert panels["conformational_state"]["gates"] == {
        "max_best_reference_rmsd_A": 2.5,
        "minimum_opposite_state_margin_A": 0.25,
        "confident_opposite_state_calls_max": 0,
        "correct_interpretable_coverage_per_state_min": 0.8,
    }


def test_empty_v2_panel_is_blocked_not_vacuously_passed():
    manifest = json.loads((QUALIFICATION / "PANEL_MANIFEST.json").read_text(encoding="utf-8"))
    module = runner_module()
    pairs = [(panel, module.validate_panel(panel)[0]) for panel in manifest["panels"]]

    # The guard this test exists for: a panel with no adopted records must be
    # reported blocked with its requirements outstanding, never vacuously passed
    # because there was nothing to check. Panels that have since been adopted
    # are held to the opposite standard.
    empty = [(p, s) for p, s in pairs if not p.get("records")]
    assert empty, "expected at least one unadopted panel while the collection is incomplete"
    for panel, summary in empty:
        assert summary["state"] == "blocked_panel_incomplete", panel["panel_id"]
        assert any(row["missing_count"] > 0 for row in summary["requirements"]), panel["panel_id"]
        assert all(row["observed_count"] == 0 for row in summary["requirements"]), panel["panel_id"]

    for panel, summary in ((p, s) for p, s in pairs if p.get("records")):
        assert not summary["errors"], f"{panel['panel_id']}: {summary['errors']}"
        assert all(row["passed"] for row in summary["requirements"]), panel["panel_id"]

    status = json.loads((QUALIFICATION / "results" / "QUALIFICATION_V2_STATUS.json").read_text(encoding="utf-8"))
    # The collection stays blocked while any panel is unadopted, and the
    # composition audit never claims scientific execution regardless.
    assert status["overall_state"] == "blocked_panel_incomplete"
    assert status["scientific_execution_performed"] is False


def test_v2_source_lock_preserves_v1_as_candidate_only():
    module = runner_module()
    lock = json.loads((QUALIFICATION / "SOURCE_LOCK.json").read_text(encoding="utf-8"))
    verified = module.verify_source_lock(lock)
    assert lock["prior_collection"]["adoption_state"] == "candidate_only"
    assert verified["prior_collection_checksum_valid"]
    # v1 stays a candidate collection: adopting sources into v2 must never mean
    # silently reusing v1's artifact tree. This once asserted that nothing was
    # adopted at all, which was only a proxy for the same intent and stopped
    # being true when the x_ray stratum was adopted.
    #
    # Only structural invariants are checked here. Whether the artifacts are
    # present and match their digests depends on acquisition, which this suite
    # deliberately does not perform -- third-party artifacts are never committed.
    # That check belongs to the qualification workflow, which acquires first.
    for source in lock["sources"]:
        assert len(str(source.get("sha256", ""))) == 64, f"{source['source_id']} has no digest"
        assert source.get("url") or source.get("acquisition") == "committed_in_repository", (
            f"{source['source_id']} is neither acquirable nor committed")
    assert not any(
        str(source["artifact"]).startswith("../qualification-v1/")
        or "qualification-v1" in str(source["artifact"])
        for source in lock["sources"]
    ), "v2 must acquire its own artifacts rather than reach into the v1 tree"
