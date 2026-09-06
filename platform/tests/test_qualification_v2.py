from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
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
    # Collection 2.4 moved membrane_orientation:beta_barrel to non-blocking:
    # under the 2.3 accuracy gate it passes 5 of 16 and cannot be qualified for
    # Mark 1. Both membrane strata are now non-blocking; the panel still executes.
    assert len(manifest["release_blocking_scopes"]) == 5
    assert "membrane_orientation" not in {s.split(":", 1)[0] for s in manifest["release_blocking_scopes"]}
    assert "membrane_orientation:alpha_helical" in manifest["non_blocking_scopes"]
    assert "membrane_orientation:beta_barrel" in manifest["non_blocking_scopes"]
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

    # Composed panels are checked for composition only. This suite never
    # acquires artifacts -- third-party files are not committed -- so every
    # record legitimately reports its artifact as absent here. Whether the bytes
    # are present and match is the qualification workflow's job, which acquires
    # first. Filtering that one error class is what keeps this test honest
    # offline instead of quietly requiring a warm working copy.
    artifact_error = "source artifact is missing or checksum-mismatched"
    shortfall = re.compile(r"requires \d+, observed \d+$")
    for panel, summary in ((p, s) for p, s in pairs if p.get("records")):
        # A panel may be adopted a stratum at a time. The membrane panel carries
        # its release-blocking beta_barrel stratum while alpha_helical is still
        # outstanding, so "has records" no longer implies "fully composed" and a
        # shortfall on an unadopted stratum is expected rather than a defect.
        composition_errors = [e for e in summary["errors"]
                              if artifact_error not in e and not shortfall.search(e)]
        assert not composition_errors, f"{panel['panel_id']}: {composition_errors}"

        # What must hold in every case: a requirement is either fully met or
        # untouched. A partially filled requirement would mean a stratum was
        # adopted half-way, which no split is allowed to be.
        for row in summary["requirements"]:
            assert row["observed_count"] in (0, row["count"]), (
                f"{panel['panel_id']}: {row['stratum']}/{row['split']} is partially adopted "
                f"({row['observed_count']} of {row['count']})")
        assert any(row["passed"] for row in summary["requirements"]), panel["panel_id"]

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


def summarizer_module():
    # It imports `run_execution` as a sibling, the way CI runs it: from inside
    # the qualification directory. Loading it by path has to put that directory
    # on sys.path or the import fails for a reason unrelated to the test.
    if str(QUALIFICATION) not in sys.path:
        sys.path.insert(0, str(QUALIFICATION))
    spec = importlib.util.spec_from_file_location(
        "qualification_v2_summarizer", QUALIFICATION / "summarize_execution.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_release_status_digest_chain_matches_the_files():
    # Adoption protocol rule 2: a digest nothing compares is a claim, not a
    # check. Three of the eight recorded here were stale on 2026-09-02 -- the
    # manifest had been revised twice and neither the composition audit nor the
    # execution summary had been regenerated against it.
    spec = importlib.util.spec_from_file_location(
        "verify_release_status_digests", ROOT / "tools" / "verify_release_status_digests.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.main([]) == 0, "RELEASE_STATUS.json records a digest that no longer matches its file"


def test_lock_health_digests_the_way_the_acquirer_does():
    """The rule that made a hand-run audit report nineteen false drifts.

    acquire_sources.py decompresses a .gz URL whose artifact path is not .gz
    before hashing. gzip embeds an mtime, so comparing the compressed bytes
    reports every validation report as drifted when none of them changed. The
    health check has to hash what the acquirer hashes or it is an alarm
    generator.
    """
    import gzip

    spec = importlib.util.spec_from_file_location(
        "verify_source_lock_health", ROOT / "tools" / "verify_source_lock_health.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    payload = b"><header>\nMKVLAA\n"
    # Two archives of identical content, written a second apart, differ in bytes.
    early = gzip.compress(payload, mtime=1)
    late = gzip.compress(payload, mtime=2)
    assert early != late, "gzip no longer embeds an mtime; this test's premise is gone"

    plain = hashlib.sha256(payload).hexdigest()
    # A .gz url whose artifact is not .gz: decompress, so both archives agree.
    assert module.digest_as_acquirer_would(early, "https://x/a.xml.gz", "sources/a.xml") == plain
    assert module.digest_as_acquirer_would(late, "https://x/a.xml.gz", "sources/a.xml") == plain
    # A .gz artifact is stored compressed, so the archive itself is the subject.
    assert module.digest_as_acquirer_would(early, "https://x/a.gz", "sources/a.gz") == \
        hashlib.sha256(early).hexdigest()


def test_status_prose_quotes_the_executed_numbers():
    # The digest chain above cannot reach this. A digest can be perfectly current
    # beside a sentence that is two collections stale, and on 2026-09-03 it was:
    # `scientific_execution_note` said membrane was 14/16 while
    # EXECUTION_SUMMARY.json had recorded 5/16 since collection 2.3 added the OPM
    # accuracy gate, and `external_benchmarks.membrane_orientation` in the same
    # file carried the correct figure the whole time.
    spec = importlib.util.spec_from_file_location(
        "verify_status_note_numbers", ROOT / "tools" / "verify_status_note_numbers.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.main([]) == 0, "RELEASE_STATUS.json prose disagrees with the executed evidence"

    # Rule 2: prove the check can fail. The exact text that stood for two
    # collections must be rejected, and so must a dropped panel figure.
    status = json.loads((ROOT / "yauvi-structural-workbench" / "RELEASE_STATUS.json").read_text(encoding="utf-8"))
    summary = json.loads((QUALIFICATION / "results" / "EXECUTION_SUMMARY.json").read_text(encoding="utf-8"))
    current = status["qualification_evidence"]["current_v2"]
    external = status["external_benchmarks"]["membrane_orientation"]

    superseded = current["scientific_execution_note_correction"]["superseded_text"]
    assert any("membrane" in p for p in module.check(superseded, external, summary))

    dropped = current["scientific_execution_note"].replace("assembly-context 16/16 with 6/6; ", "")
    assert any("assembly_interface" in p for p in module.check(dropped, external, summary))

    assert module.check(current["scientific_execution_note"], external, summary) == []


def test_non_blocking_scope_failure_does_not_gate_the_release():
    # Collection 2.4 made membrane orientation research-only, and until
    # 2026-09-02 the summarizer still failed the whole run on it: every CI job
    # on every push and the weekly schedule was red for a scope the manifest
    # says does not gate the release.
    manifest = json.loads((QUALIFICATION / "PANEL_MANIFEST.json").read_text(encoding="utf-8"))
    non_blocking = frozenset(manifest["non_blocking_scopes"])
    module = summarizer_module()

    assert module.release_blocking("membrane_orientation", ["beta_barrel"], non_blocking) is False
    assert module.release_blocking("membrane_orientation", ["alpha_helical"], non_blocking) is False
    assert module.release_blocking("structure_qc", ["x_ray", "cryo_em"], non_blocking) is True

    # conformational_state has one blocking stratum and one non-blocking one.
    # Executing the non-blocking stratum must not exempt the workflow, and a
    # panel that recorded no strata at all fails closed.
    assert module.release_blocking("conformational_state", ["other_proteins"], non_blocking) is False
    assert module.release_blocking(
        "conformational_state", ["abl_family", "other_proteins"], non_blocking
    ) is True
    assert module.release_blocking("conformational_state", [], non_blocking) is True


def test_execution_summary_reports_which_failures_gate_the_release():
    summary = json.loads(
        (QUALIFICATION / "results" / "EXECUTION_SUMMARY.json").read_text(encoding="utf-8")
    )
    # Both facts have to be readable without inferring either from the other:
    # a panel failed, and no panel that gates the release failed.
    assert summary["every_executed_panel_passed"] is False
    assert summary["every_executed_release_blocking_panel_passed"] is True
    assert summary["release_blocking_panels_failed"] == []
    assert summary["non_blocking_panels_failed"] == ["membrane_orientation"]
    # The gate this summary must never be able to close, whatever it executed.
    assert summary["all_release_blocking_scopes_qualified"] is False
    assert summary["second_machine_reproduction"] == "not_recorded"

    # Required-case totals are derived from the manifest, not typed. The
    # 2026-09-01 ABL revision moved this from 114 to 110 and the committed
    # summary kept quoting 114.
    manifest = json.loads((QUALIFICATION / "PANEL_MANIFEST.json").read_text(encoding="utf-8"))
    required = sum(
        requirement["count"] for panel in manifest["panels"] for requirement in panel["requirements"]
    )
    assert summary["cases_required"] == required


def _cross_machine_module():
    spec = importlib.util.spec_from_file_location(
        "verify_cross_machine_reproduction",
        ROOT / "tools" / "verify_cross_machine_reproduction.py",
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _runner_summary(platform_name, machine, python, panels):
    """One runner's EXECUTION_SUMMARY.json, reduced to what the combiner reads."""
    return {
        "source": f"{platform_name}-{machine}-py{python}",
        "summary": {
            "recorded_on": [{"platform": platform_name, "machine": machine, "python": python}],
            "panels": [
                {
                    "workflow": workflow,
                    "release_blocking": blocking,
                    "stratum_state": "passed" if failed == 0 else "failed",
                    "cases": {"passed": passed, "failed": failed, "total": passed + failed},
                    "controls": {"passed": 1, "total": 1},
                }
                for workflow, blocking, passed, failed in panels
            ],
        },
    }


# The panel shape the four adopted blocking panels share: everything passes.
_CLEAN = [("structure_qc", True, 16, 0), ("conformational_state", True, 14, 0)]


def test_cross_machine_reproduction_needs_two_architectures_not_two_pythons():
    """Rule 2, on the gate that would otherwise close itself.

    Six green runners are not six machines. The qualification matrix varies
    Python three ways against two operating systems, and counting runners would
    report six independent reproductions of what are really two environments.
    A different interpreter on the same platform and architecture does not
    exercise the libm, BLAS or floating-point-contraction differences this gate
    exists to catch.
    """
    module = _cross_machine_module()

    # Three Pythons, one OS, one architecture. Three runners, one environment.
    one_environment = [
        _runner_summary("Linux", "x86_64", version, _CLEAN) for version in ("3.10", "3.11", "3.12")
    ]
    thin = module.build_report(one_environment, min_environments=2)
    assert thin["environment_count"] == 1, thin["environments_observed"]
    assert thin["runners_usable"] == 3
    assert thin["every_release_blocking_panel_reproduced"] is False
    assert sorted(thin["release_blocking_panels_under_minimum"]) == [
        "conformational_state", "structure_qc"
    ]

    # Add a second architecture and the same evidence now reproduces.
    two_environments = one_environment + [
        _runner_summary("Darwin", "arm64", version, _CLEAN) for version in ("3.10", "3.11", "3.12")
    ]
    wide = module.build_report(two_environments, min_environments=2)
    assert wide["environment_count"] == 2
    assert wide["environments_observed"] == ["Darwin/arm64", "Linux/x86_64"]
    assert wide["every_release_blocking_panel_reproduced"] is True
    assert wide["release_blocking_panels_under_minimum"] == []


def test_cross_machine_reproduction_fails_when_environments_disagree():
    """Both green is not the same as both agreeing.

    A panel that passes 16/16 on Linux and 14/16 on macOS has not reproduced,
    and the failure mode worth guarding against is a gate that reads only the
    verdict and calls two different measurements a match.
    """
    module = _cross_machine_module()

    disagreeing = [
        _runner_summary("Linux", "x86_64", "3.12", [("structure_qc", True, 16, 0)]),
        _runner_summary("Darwin", "arm64", "3.12", [("structure_qc", True, 14, 2)]),
    ]
    report = module.build_report(disagreeing, min_environments=2)
    assert report["environment_count"] == 2, "two environments were present"
    assert report["release_blocking_panels_disagreed"] == ["structure_qc"]
    assert report["every_release_blocking_panel_reproduced"] is False
    panel = report["panels"][0]
    assert panel["agreed_across_environments"] is False
    assert panel["reproduced"] is False

    # The subtler case: both report "passed", but not the same passed. Counting
    # verdicts would call this reproduced; comparing outcomes does not.
    same_verdict_different_counts = [
        _runner_summary("Linux", "x86_64", "3.12", [("structure_qc", True, 16, 0)]),
        _runner_summary("Darwin", "arm64", "3.12", [("structure_qc", True, 15, 0)]),
    ]
    subtle = module.build_report(same_verdict_different_counts, min_environments=2)
    assert subtle["panels"][0]["passed_on_every_environment"] is True
    assert subtle["panels"][0]["reproduced"] is False
    assert subtle["release_blocking_panels_disagreed"] == ["structure_qc"]


def test_cross_machine_reproduction_does_not_gate_on_a_non_blocking_scope():
    """Membrane orientation is expected to differ across machines.

    Its drift deltas are the reason the matrix exists, and collection 2.4 made
    the scope research-only. A cross-machine disagreement there is the finding,
    not a broken release -- the same split summarize_execution.py enforces.
    """
    module = _cross_machine_module()

    report = module.build_report(
        [
            _runner_summary("Linux", "x86_64", "3.12",
                            _CLEAN + [("membrane_orientation", False, 5, 11)]),
            _runner_summary("Darwin", "arm64", "3.12",
                            _CLEAN + [("membrane_orientation", False, 4, 12)]),
        ],
        min_environments=2,
    )
    membrane = next(p for p in report["panels"] if p["workflow"] == "membrane_orientation")
    assert membrane["agreed_across_environments"] is False, "the disagreement is still recorded"
    assert membrane["release_blocking"] is False
    assert report["release_blocking_panels_disagreed"] == []
    assert report["every_release_blocking_panel_reproduced"] is True


def test_cross_machine_reproduction_catches_one_environment_disagreeing_with_itself():
    """Two Pythons on one machine reporting different counts is a finding.

    Collapsing an environment to a single representative outcome would hide it,
    so the disagreement is recorded before the per-environment reduction.
    """
    module = _cross_machine_module()

    report = module.build_report(
        [
            _runner_summary("Linux", "x86_64", "3.10", [("structure_qc", True, 16, 0)]),
            _runner_summary("Linux", "x86_64", "3.12", [("structure_qc", True, 15, 1)]),
            _runner_summary("Darwin", "arm64", "3.12", [("structure_qc", True, 16, 0)]),
        ],
        min_environments=2,
    )
    panel = report["panels"][0]
    assert panel["internally_inconsistent_environments"] == ["Linux/x86_64"]
    assert panel["reproduced"] is False
    assert report["release_blocking_panels_disagreed"] == ["structure_qc"]


def test_cross_machine_reproduction_rejects_a_summary_spanning_two_runtimes():
    """A summary carrying two runtimes was not produced by one run.

    summarize_execution.py sets counts_are_single_machine false in exactly that
    case. Such a document cannot speak for an environment, and silently letting
    it count as one would manufacture a second machine out of a merged local
    directory.
    """
    module = _cross_machine_module()

    merged = {
        "source": "merged",
        "summary": {
            "recorded_on": [
                {"platform": "Linux", "machine": "x86_64", "python": "3.12"},
                {"platform": "Darwin", "machine": "arm64", "python": "3.12"},
            ],
            "panels": [{"workflow": "structure_qc", "release_blocking": True,
                        "stratum_state": "passed",
                        "cases": {"passed": 16, "failed": 0, "total": 16},
                        "controls": {"passed": 1, "total": 1}}],
        },
    }
    report = module.build_report([merged], min_environments=2)
    assert report["runners_usable"] == 0
    assert report["unusable_runners"][0]["why"].startswith("summary spans 2 runtimes")
    assert report["every_release_blocking_panel_reproduced"] is False


def test_cross_machine_missing_evidence_is_not_a_failed_reproduction(tmp_path):
    """Exit 2, not 1, when the evidence simply did not arrive.

    Uploads in qualification.yml are continue-on-error because the artifact
    service is not evidence about the software -- on run #46 a CreateArtifact
    call timed out through all five retries on a runner whose panels had all
    passed. If a dropped upload came back as a red reproduction gate, the
    workflow would state something false about the science for an infrastructure
    reason. verify_source_lock_health.py draws the same line, exiting 0 when a
    provider is merely unreachable.
    """
    import io
    import tarfile

    module = _cross_machine_module()

    def artifact(directory, platform_name, machine, python, cases_passed):
        directory.mkdir(parents=True, exist_ok=True)
        blob = json.dumps(_runner_summary(
            platform_name, machine, python,
            [("structure_qc", True, cases_passed, 16 - cases_passed)],
        )["summary"]).encode("utf-8")
        with tarfile.open(directory / "execution-evidence.tar.gz", "w:gz") as archive:
            info = tarfile.TarInfo("results/EXECUTION_SUMMARY.json")
            info.size = len(blob)
            archive.addfile(info, io.BytesIO(blob))

    # One environment's upload arrived. Nothing disagrees; there is just not
    # enough to conclude with.
    thin = tmp_path / "thin"
    artifact(thin / "execution-evidence-py3.12-ubuntu-latest", "Linux", "x86_64", "3.12", 16)
    assert module.main(["--evidence-dir", str(thin),
                        "--json-out", str(thin / "report.json")]) == 2

    # Both arrived and they disagree. That is a finding, and it is louder.
    split = tmp_path / "split"
    artifact(split / "execution-evidence-py3.12-ubuntu-latest", "Linux", "x86_64", "3.12", 16)
    artifact(split / "execution-evidence-py3.12-macos-latest", "Darwin", "arm64", "3.12", 14)
    assert module.main(["--evidence-dir", str(split),
                        "--json-out", str(split / "report.json")]) == 1

    # Both arrived and they agree.
    whole = tmp_path / "whole"
    artifact(whole / "execution-evidence-py3.12-ubuntu-latest", "Linux", "x86_64", "3.12", 16)
    artifact(whole / "execution-evidence-py3.12-macos-latest", "Darwin", "arm64", "3.12", 16)
    assert module.main(["--evidence-dir", str(whole),
                        "--json-out", str(whole / "report.json")]) == 0

    # The archives are what CI uploads, so the reader has to survive them.
    report = json.loads((whole / "report.json").read_text(encoding="utf-8"))
    assert report["environments_observed"] == ["Darwin/arm64", "Linux/x86_64"]
    assert report["every_release_blocking_panel_reproduced"] is True
