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



# Compact synthetic protocol: independent cases, a negative control, and one
# optional experimental panel. No real qualification expectations are changed.
def _contract_fixture():
    identity = {k: 'a' * 64 for k in ('code_sha256', 'manifest_sha256', 'source_lock_sha256', 'protocols_sha256')}
    identity['execution_panels'] = {'structure_qc': 'b' * 64}
    manifest = {'release_blocking_scopes': ['structure_qc:x_ray'], 'panels': [
        {'workflow': 'structure_qc', 'requirements': [{'stratum': 'x_ray', 'count': 2}],
         'records': [{'record_id': x, 'stratum': 'x_ray'} for x in ('case-a', 'case-b')],
         'controls': [{'record_id': 'control-a'}]}]}
    def runner(system, machine):
        return {'source': system, 'summary': {'recorded_on': [{'platform': system, 'machine': machine, 'python': '3.12'}],
            'panels': [{'workflow': 'structure_qc', 'stratum_state': 'passed', 'evidence_identity': identity.copy(),
                'execution_panel_sha256': 'b' * 64, 'cases': {'passed': 2, 'failed': 0, 'total': 2},
                'controls': {'passed': 1, 'total': 1}, 'coverage': {'unmet': []},
                'record_verdicts': {kind: [{'record_id': x, 'passed': True,
                  'checks': [{'check': 'independent-reference', 'required': True, 'passed': True}]} for x in ids]
                  for kind, ids in [('cases', ['case-a', 'case-b']), ('controls', ['control-a'])]}}]}}
    return manifest, identity, [runner('Linux', 'x86_64'), runner('Darwin', 'arm64')]


def _report(change=None):
    manifest, identity, rows = _contract_fixture()
    if change:
        change(manifest, rows)
    return _cross_machine_module().build_report(rows, manifest=manifest, identity=identity)


def test_cross_machine_complete_identity_bound_evidence_passes():
    report = _report()
    assert report['exit_code'] == 0
    assert report['every_release_blocking_panel_reproduced'] is True


def _fail_panel(row):
    p = row['summary']['panels'][0]
    p['stratum_state'] = 'failed'
    p['cases'].update(passed=1, failed=1)
    p['record_verdicts']['cases'][0]['passed'] = False
    p['record_verdicts']['cases'][0]['checks'][0]['passed'] = False


def test_matching_failures_never_pass_aggregate():
    report = _report(lambda m, rows: [_fail_panel(r) for r in rows])
    assert report['release_blocking_panels_disagreed'] == []
    assert report['release_blocking_panels_failed'] == ['structure_qc']
    assert report['every_release_blocking_panel_reproduced'] is False
    assert report['exit_code'] == 1


def test_different_case_outcomes_fail():
    report = _report(lambda m, rows: _fail_panel(rows[0]))
    assert report['release_blocking_panels_disagreed'] == ['structure_qc']
    assert report['exit_code'] == 1


def test_python_versions_do_not_supply_second_environment():
    def change(m, rows):
        rows[1]['summary']['recorded_on'][0].update(platform='Linux', machine='x86_64', python='3.11')
    assert _report(change)['exit_code'] == 2


def test_missing_expected_panel_fails_closed():
    def change(m, rows):
        m['release_blocking_scopes'].append('sf_csa:enzymes')
        m['panels'].append({'workflow': 'sf_csa', 'records': [], 'requirements': [{'stratum':'enzymes','count':16}]})
    r = _report(change)
    assert 'sf_csa' in r['release_blocking_panels_incomplete']
    assert r['exit_code'] == 2


def test_corrupt_duplicate_missing_and_mixed_evidence_is_rejected():
    mutations = [
        lambda p: p['cases'].update(passed=99),
        lambda p: p['record_verdicts']['cases'].append(p['record_verdicts']['cases'][0]),
        lambda p: p['record_verdicts']['cases'].pop(),
        lambda p: p['record_verdicts']['cases'][0].update(record_id='wrong'),
        lambda p: p['evidence_identity'].update(code_sha256='c'*64),
        lambda p: p['evidence_identity'].update(source_lock_sha256='c'*64),
        lambda p: p.update(execution_panel_sha256='c'*64),
        lambda p: p.pop('evidence_identity'),
        lambda p: p['record_verdicts']['cases'][0].update(passed=False),
    ]
    for mutation in mutations:
        r = _report(lambda m, rows: mutation(rows[0]['summary']['panels'][0]))
        assert r['exit_code'] == 2
        assert not r['every_release_blocking_panel_reproduced']


def test_duplicate_or_missing_runner_is_incomplete():
    assert _report(lambda m, rows: rows.append(rows[0]))['exit_code'] == 2
    assert _report(lambda m, rows: rows.pop())['exit_code'] == 2


def test_control_failure_is_not_hidden_by_passing_cases():
    def change(m, rows):
        for row in rows:
            p = row['summary']['panels'][0]
            p['controls']['passed'] = 0
            p['record_verdicts']['controls'][0]['passed'] = False
            p['record_verdicts']['controls'][0]['checks'][0]['passed'] = False
    assert _report(change)['exit_code'] == 1


def test_archive_reader_never_extracts_unsafe_members(tmp_path):
    import io, tarfile
    module = _cross_machine_module()
    with tarfile.open(tmp_path/'bad.tar.gz', 'w:gz') as tar:
        info = tarfile.TarInfo('../EXECUTION_SUMMARY.json'); info.size = 2
        tar.addfile(info, io.BytesIO(b'{}'))
    rows = module.load_summaries(tmp_path)
    assert rows[0]['error']
    assert not (tmp_path.parent/'EXECUTION_SUMMARY.json').exists()
