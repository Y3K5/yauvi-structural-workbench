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
    summaries = [module.validate_panel(panel)[0] for panel in manifest["panels"]]
    assert all(summary["state"] == "blocked_panel_incomplete" for summary in summaries)
    assert all(any(row["missing_count"] > 0 for row in summary["requirements"]) for summary in summaries)
    status = json.loads((QUALIFICATION / "results" / "QUALIFICATION_V2_STATUS.json").read_text(encoding="utf-8"))
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
    assert verified["adopted_sources_valid"]
    assert not any(
        str(source["artifact"]).startswith("../qualification-v1/")
        or "qualification-v1" in str(source["artifact"])
        for source in lock["sources"]
    ), "v2 must acquire its own artifacts rather than reach into the v1 tree"
