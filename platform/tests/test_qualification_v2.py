from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import re
import subprocess
import sys

import pytest
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
    # because there was nothing to check.
    #
    # It used to require the manifest to contain such a panel, and asserted so
    # while one did. sf_csa was the last, and collection 2.11 adopted it -- at
    # which point the assertion failed for the best possible reason and the guard
    # would otherwise have been deleted along with the failure. A check that
    # retires itself the moment everything is adopted stops protecting anything
    # exactly when the stakes are highest, so it is exercised on a synthetic
    # emptied panel instead. Real requirements, real shape, no records.
    for panel in manifest["panels"]:
        emptied = copy.deepcopy(panel)
        emptied["records"] = []
        emptied["controls"] = []
        summary = module.validate_panel(emptied)[0]
        assert summary["state"] == "blocked_panel_incomplete", panel["panel_id"]
        assert any(row["missing_count"] > 0 for row in summary["requirements"]), panel["panel_id"]
        assert all(row["observed_count"] == 0 for row in summary["requirements"]), panel["panel_id"]

    # And nothing in the real manifest is empty any more, which is the state
    # collection 2.11 reached. If a future collection adds an unadopted panel,
    # the loop above still covers it.
    assert all(p.get("records") for p, _ in pairs), \
        "an unadopted panel reappeared: " + str([p["panel_id"] for p, _ in pairs if not p.get("records")])

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


def _report(change=None, drafts=None):
    manifest, identity, rows = _contract_fixture()
    if change:
        change(manifest, rows)
    return _cross_machine_module().build_report(
        rows, manifest=manifest, identity=identity, drafts=drafts)


def _add_draft_panel(manifest, identity, rows, workflow='sf_csa'):
    """A release-blocking scope the manifest names but has not populated.

    This is the exact shape that discarded seven runners: the workflow executes
    the panel from its adoption draft, because reproduction has to exist before
    adoption can be granted, and the manifest therefore holds no records for it.
    """
    manifest['release_blocking_scopes'].append(workflow + ':curated')
    manifest['panels'].append({'workflow': workflow, 'requirements': [], 'records': [], 'controls': []})
    identity['execution_panels'][workflow] = 'c' * 64
    draft = {'workflow': workflow, 'draft_state': 'not_adopted',
             'records': [{'record_id': 'draft-a', 'stratum': 'curated'},
                         {'record_id': 'draft-b', 'stratum': 'curated'}],
             'controls': []}
    for row in rows:
        row['summary']['panels'].append({
            'workflow': workflow, 'stratum_state': 'passed',
            'evidence_identity': identity.copy(), 'execution_panel_sha256': 'c' * 64,
            'cases': {'passed': 2, 'failed': 0, 'total': 2},
            'controls': {'passed': 0, 'total': 0}, 'coverage': {'unmet': []},
            'record_verdicts': {'cases': [{'record_id': x, 'passed': True,
                'checks': [{'check': 'independent-reference', 'required': True, 'passed': True}]}
                for x in ('draft-a', 'draft-b')], 'controls': []}})
    return {workflow: draft}


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


def _acquirer_module():
    spec = importlib.util.spec_from_file_location(
        "acquire_sources",
        ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2" / "acquire_sources.py",
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fake_response(body: bytes, content_type: str):
    import email.message, io

    headers = email.message.Message()
    headers["Content-Length"] = str(len(body))
    headers["Content-Type"] = content_type

    class Response(io.BytesIO):
        def __init__(self):
            super().__init__(body)
            self.headers = headers

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.close()

    return Response()


def test_a_definitive_provider_answer_is_not_retried(tmp_path, monkeypatch):
    """What made a broken acquisition step unreadable in CI.

    Every failure -- a 404, a truncated read, a landing page where a gzip stream
    was expected -- collapsed into "6 attempts failed", after fifty seconds of
    backoff per artifact. A status that will not change on the seventh ask has to
    be reported, once, with its number.
    """
    import urllib.error

    module = _acquirer_module()

    def refuse(*_args, **_kwargs):
        raise urllib.error.HTTPError("https://example.invalid/x.fasta.gz", 404, "Not Found", None, None)

    monkeypatch.setattr(module.urllib.request, "urlopen", refuse)
    calls = []
    monkeypatch.setattr(module.time, "sleep", lambda seconds: calls.append(seconds))

    with pytest.raises(module.Unacquirable) as caught:
        module.fetch("https://example.invalid/x.fasta.gz", tmp_path / "x.fasta", "0" * 64)

    assert "404" in str(caught.value)
    assert calls == [], "a terminal status must not consume the retry budget"


def test_a_landing_page_served_for_a_gz_url_names_what_arrived(tmp_path, monkeypatch):
    """A URL that addresses a preview page rather than a file says so.

    Decompressing HTML raises "Not a gzipped file", which reads as a corrupt
    download and invites a retry. It is neither: the URL names a page.
    """
    module = _acquirer_module()
    page = b"<!doctype html><html><body>record preview</body></html>"
    monkeypatch.setattr(module.urllib.request, "urlopen",
                        lambda *a, **k: _fake_response(page, "text/html"))
    monkeypatch.setattr(module.time, "sleep", lambda seconds: None)

    with pytest.raises(module.Unacquirable) as caught:
        module.fetch("https://example.invalid/x.fasta.gz", tmp_path / "x.fasta", "0" * 64)

    message = str(caught.value)
    assert "markup" in message and "text/html" in message
    assert not (tmp_path / "x.fasta").exists(), "a refused fetch must leave no artifact behind"


def test_a_transient_status_still_exhausts_the_retry_budget(tmp_path, monkeypatch):
    """The truncation case the retry loop exists for is unchanged."""
    import urllib.error

    module = _acquirer_module()
    monkeypatch.setattr(module.urllib.request, "urlopen", lambda *a, **k: (_ for _ in ()).throw(
        urllib.error.HTTPError("https://example.invalid/x.fasta.gz", 502, "Bad Gateway", None, None)))
    slept = []
    monkeypatch.setattr(module.time, "sleep", lambda seconds: slept.append(seconds))

    with pytest.raises(RuntimeError):
        module.fetch("https://example.invalid/x.fasta.gz", tmp_path / "x.fasta", "0" * 64)

    assert len(slept) == 5, "a busy provider must still be retried"


def test_no_tracked_file_publishes_a_local_home_path():
    """The screen that reads the resolved tree rather than the change.

    51 occurrences of one home directory reached public history through a
    generated `output_dir` field, in a class the publishing protocol had already
    recorded once. Nobody wrote those lines, so no diff review could catch them.
    """
    spec = importlib.util.spec_from_file_location(
        "verify_no_local_paths_in_evidence", ROOT / "tools" / "verify_no_local_paths_in_evidence.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    disclosed, exempt = module.findings(
        ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2"
        / "results" / "execution-structqc" / "EXECUTION_STATUS.json"
    )
    assert disclosed == [], f"recorded evidence still publishes {sorted(set(disclosed))}"

    # The whole-tree pass needs a git checkout. The private edit tree is not one,
    # and the check that matters there is the per-file one above.
    if subprocess.run(["git", "rev-parse", "--git-dir"], cwd=ROOT,
                      capture_output=True).returncode != 0:
        pytest.skip("not a git checkout; the published tree is checked in CI")
    assert module.main([]) == 0, "a tracked file publishes a local home path"


def test_a_throttled_provider_is_given_the_window_it_asked_for(tmp_path, monkeypatch):
    """A 150 MB serial acquisition gets throttled, and 20 seconds is not the ask.

    The backoff tops out below the window a repository usually names, so the whole
    retry budget could be spent re-asking inside a throttle that had not lifted.
    """
    import email.message
    import urllib.error

    module = _acquirer_module()
    headers = email.message.Message()
    headers["Retry-After"] = "45"
    monkeypatch.setattr(module.urllib.request, "urlopen", lambda *a, **k: (_ for _ in ()).throw(
        urllib.error.HTTPError("https://example.invalid/x.fasta.gz", 429, "Too Many Requests",
                               headers, None)))
    slept = []
    monkeypatch.setattr(module.time, "sleep", lambda seconds: slept.append(seconds))

    with pytest.raises(RuntimeError):
        module.fetch("https://example.invalid/x.fasta.gz", tmp_path / "x.fasta", "0" * 64)

    assert slept and slept[0] == 45, f"honoured the provider's window, got {slept}"
    assert max(slept) <= module.MAX_RETRY_AFTER, "one artifact must not hold the job open"


def test_a_mirror_carries_the_traffic_and_the_archive_stays_the_record():
    """Mirrors are tried first, the citable archive last.

    Seven runners re-acquiring a 146 MB set on every push aimed roughly a
    gigabyte a day at a preservation archive, which began returning gateway
    timeouts. The location a manifest cites and the location CI hammers should
    not be the same one.
    """
    module = _acquirer_module()
    entry = {
        "artifact": "sources/proteomes/UP000000579.fasta",
        "url": "https://zenodo.org/records/22652863/files/UP000000579.fasta.gz",
        "mirrors": ["https://huggingface.co/datasets/x/y/resolve/abc123/UP000000579.fasta.gz"],
    }
    assert module.candidates(entry) == [entry["mirrors"][0], entry["url"]]

    # A bare string is accepted, and a mirror repeating the archive is not tried twice.
    assert module.candidates({"url": "https://a/x", "mirrors": "https://b/x"}) == \
        ["https://b/x", "https://a/x"]
    assert module.candidates({"url": "https://a/x", "mirrors": ["https://a/x"]}) == ["https://a/x"]
    assert module.candidates({"url": "https://a/x"}) == ["https://a/x"]
    assert module.candidates({}) == []


def test_a_failing_archive_falls_through_to_its_mirror(tmp_path, monkeypatch):
    """And the run still says the archive failed, rather than reporting silence."""
    import gzip
    import urllib.error

    module = _acquirer_module()
    payload = b">sp|TEST\nMKVLAA\n"
    digest = hashlib.sha256(payload).hexdigest()

    def urlopen(request, *_args, **_kwargs):
        url = request.full_url if hasattr(request, "full_url") else request
        if "zenodo" in url:
            raise urllib.error.HTTPError(url, 504, "Gateway Time-out", None, None)
        return _fake_response(gzip.compress(payload), "application/gzip")

    monkeypatch.setattr(module.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(module.time, "sleep", lambda seconds: None)

    entry = {
        "artifact": "sources/proteomes/x.fasta",
        "sha256": digest,
        "mirrors": ["https://huggingface.co/datasets/x/y/resolve/abc/x.fasta.gz"],
        "url": "https://zenodo.org/records/1/files/x.fasta.gz",
    }
    served_by, attempted = module.acquire(entry, tmp_path / "x.fasta")

    assert "huggingface" in served_by
    assert (tmp_path / "x.fasta").read_bytes() == payload
    assert attempted == [], "the mirror is tried first, so nothing failed before it"

    # With the order reversed, the archive's failure is carried into the report.
    entry_archive_first = {**entry, "mirrors": [], "url": entry["url"]}
    entry_archive_first["mirrors"] = []
    with pytest.raises(RuntimeError) as caught:
        module.acquire(entry_archive_first, tmp_path / "y.fasta")
    assert "504" in str(caught.value)


def test_a_mirror_serving_different_bytes_is_refused(tmp_path, monkeypatch):
    """A mirror can never weaken the lock; it is held to the same digest."""
    import gzip

    module = _acquirer_module()
    monkeypatch.setattr(module.urllib.request, "urlopen",
                        lambda *a, **k: _fake_response(gzip.compress(b"different\n"),
                                                       "application/gzip"))
    monkeypatch.setattr(module.time, "sleep", lambda seconds: None)

    entry = {
        "artifact": "sources/proteomes/x.fasta",
        "sha256": hashlib.sha256(b">sp|TEST\nMKVLAA\n").hexdigest(),
        "mirrors": ["https://huggingface.co/datasets/x/y/resolve/abc/x.fasta.gz"],
        "url": "https://zenodo.org/records/1/files/x.fasta.gz",
    }
    with pytest.raises(RuntimeError) as caught:
        module.acquire(entry, tmp_path / "x.fasta")
    assert "digest mismatch" in str(caught.value)
    assert not (tmp_path / "x.fasta").exists()


def test_health_reports_every_promise_the_lock_makes():
    """Health checks all candidates; a live mirror must not mask a dead archive."""
    spec = importlib.util.spec_from_file_location(
        "verify_source_lock_health", ROOT / "tools" / "verify_source_lock_health.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    entry = {"url": "https://archive/x.gz", "mirrors": ["https://mirror/x.gz"]}
    assert module.locked_urls(entry) == ["https://archive/x.gz", "https://mirror/x.gz"]
    assert module.locked_urls({"mirrors": ["https://only-a-mirror/x.gz"]}) == \
        ["https://only-a-mirror/x.gz"]


def test_the_qualification_cache_key_is_the_lock_itself():
    """Caching acquired sources is only safe if the key moves with the lock.

    A key that outlived a lock change would restore artifacts belonging to a
    different lock, and acquisition skips what is already present -- so the run
    would quietly execute against superseded sources.
    """
    workflow = (ROOT / ".github" / "workflows" / "qualification.yml").read_text(encoding="utf-8")
    assert "actions/cache/restore@v4" in workflow
    assert "hashFiles('yauvi-structural-workbench/benchmarks/qualification-v2/SOURCE_LOCK.json')" \
        in workflow, "the cache key must be derived from the lock"
    # As a key, not as the comment explaining why it is absent.
    assert not any(line.strip().startswith("restore-keys:") for line in workflow.splitlines()), \
        "a partial key match would restore artifacts belonging to a different lock"
    # The save must survive a red run, or the matrix can never bank what it fetched.
    assert "actions/cache/save@v4" in workflow
    assert "if: always() && steps.sources-cache.outputs.cache-hit != 'true'" in workflow


def test_an_unadopted_panel_does_not_void_the_runners_that_carry_it():
    """The deadlock that discarded seven usable runners.

    sf-csa executed from its adoption draft, so its 16 record IDs were absent
    from the manifest, so `_records` raised -- and because that validation sat
    inside the per-runner load, one unadopted panel threw out the runner's whole
    summary. All seven runners were unusable and the five adopted panels lost
    their reproduction with it. Adoption became unreachable through the only
    second machine the project has.
    """
    module = _cross_machine_module()
    manifest, identity, rows = _contract_fixture()
    drafts = _add_draft_panel(manifest, identity, rows)
    report = module.build_report(rows, manifest=manifest, identity=identity, drafts=drafts)

    assert report['runners_usable'] == len(rows), report['unusable_runners']
    assert report['unusable_runners'] == []
    adopted = {p['workflow']: p for p in report['panels']}
    assert adopted['structure_qc']['reproduced'] is True, "the adopted panel still reproduces"


def test_a_draft_panel_is_scored_but_never_counts_as_qualified():
    """Reproduced and still not adopted are different states, and stay different."""
    module = _cross_machine_module()
    manifest, identity, rows = _contract_fixture()
    drafts = _add_draft_panel(manifest, identity, rows)
    report = module.build_report(rows, manifest=manifest, identity=identity, drafts=drafts)

    drafted = {p['workflow']: p for p in report['draft_panels']}
    assert 'sf_csa' in drafted, "the draft panel is reported"
    assert drafted['sf_csa']['reproduced_as_draft'] is True
    assert drafted['sf_csa']['adopted'] is False
    assert drafted['sf_csa']['coverage_assessed'] is False, "a draft carries no requirements"
    assert 'reproduced' not in drafted['sf_csa'], "must not read as a qualified scope"

    # And the release bar does not shrink to whatever happens to be adopted.
    assert report['release_blocking_panels_unadopted'] == ['sf_csa']
    assert report['every_release_blocking_panel_reproduced'] is False
    assert report['exit_code'] == 2
    assert 'sf_csa' not in {p['workflow'] for p in report['panels']}


def test_one_unusable_panel_does_not_silence_the_others():
    """Per-panel isolation, without losing the fact that something was rejected."""
    module = _cross_machine_module()
    manifest, identity, rows = _contract_fixture()
    identity['execution_panels']['rogue'] = 'd' * 64
    manifest['panels'].append({'workflow': 'rogue', 'requirements': [], 'records': [], 'controls': []})
    for row in rows:
        row['summary']['panels'].append({
            'workflow': 'rogue', 'stratum_state': 'passed',
            'evidence_identity': identity.copy(), 'execution_panel_sha256': 'd' * 64,
            'cases': {'passed': 1, 'failed': 0, 'total': 1}, 'controls': {'passed': 0, 'total': 0},
            'coverage': {'unmet': []},
            'record_verdicts': {'cases': [{'record_id': 'nobody-declared-me', 'passed': True,
                'checks': [{'check': 'x', 'required': True, 'passed': True}]}], 'controls': []}})
    report = module.build_report(rows, manifest=manifest, identity=identity, drafts={})

    assert report['runners_usable'] == len(rows), "the runner survives its rogue panel"
    assert [e['workflow'] for e in report['unusable_panels']] == ['rogue'] * len(rows)
    assert {p['workflow'] for p in report['panels']} >= {'structure_qc'}


def test_a_blocking_panel_short_of_its_declared_coverage_cannot_reproduce():
    """The coverage check has teeth, and is not vacuously satisfied.

    `release_blocking_scopes` entries are `workflow:scope-label`, but this checker
    rebuilt `workflow:stratum` and tested it for membership. Nothing ever matched,
    so `requirements` was always empty, `coverage_complete` always False, and no
    release-blocking panel could reach `reproduced: True` by construction.

    Fixing that must not replace an impossible check with an empty one. A panel
    that declares a requirement it does not meet still fails.
    """
    def drop_a_record(manifest, rows):
        panel = manifest['panels'][0]
        panel['requirements'] = [{'stratum': 'x_ray', 'count': 3}]   # declares 3
        # records still hold 2, so the declared coverage is not met
    report = _report(drop_a_record)
    assert report['panels'][0]['complete'] is False
    assert report['panels'][0]['reproduced'] is False
    assert report['release_blocking_panels_incomplete'] == ['structure_qc']
    assert report['every_release_blocking_panel_reproduced'] is False


def test_declared_coverage_is_read_from_the_panel_not_from_the_scope_label():
    """A blocking panel whose requirements are met does reproduce.

    The companion to the test above: with the same fixture and its declared
    coverage satisfied, the panel reaches `reproduced`. Before the fix this was
    unreachable no matter what the evidence said.
    """
    report = _report()
    panel = report['panels'][0]
    assert panel['workflow'] == 'structure_qc'
    assert panel['complete'] is True
    assert panel['reproduced'] is True
    assert report['every_release_blocking_panel_reproduced'] is True
    assert report['exit_code'] == 0


def test_a_non_blocking_panel_is_not_held_to_release_coverage():
    """Coverage gates the release scopes, not everything that runs.

    membrane_orientation declares 32 records across four requirements and carries
    16; it is non-blocking and research-only from collection 2.4, so its coverage
    must not be assessed as if it gated a release it does not gate.
    """
    def add_non_blocking(manifest, rows):
        manifest['panels'].append({'workflow': 'membrane_orientation',
                                   'requirements': [{'stratum': 'beta_barrel', 'count': 99}],
                                   'records': [], 'controls': []})
    report = _report(add_non_blocking)
    entry = [p for p in report['panels'] if p['workflow'] == 'membrane_orientation']
    assert entry == [] or entry[0]['release_blocking'] is False
    assert report['every_release_blocking_panel_reproduced'] is True, \
        "an unmet requirement on a non-blocking panel must not gate the release"


def test_the_sf_csa_gates_are_falsifiable():
    """Protocol rule 2, run rather than filed.

    `sf_csa_gate_falsification.py` is the evidence that sf-csa's gates can fail,
    and nothing invoked it — no CI step, no test. It described the panel as it
    stood on 2026-08-31 and went stale the next day, when the RBH defect was
    repaired and `rbh_asserted_rejected` was retired. Ten days later it reported
    four failures against a panel that had simply moved on, and that only came to
    light because adoption was being considered.

    A rule-2 harness nobody runs is a claim, not a check. This runs it.
    """
    module = ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2" / "sf_csa_gate_falsification.py"
    assert module.is_file(), module
    result = subprocess.run([sys.executable, str(module)], capture_output=True, text=True,
                            cwd=module.parent)
    assert result.returncode == 0, (
        "sf-csa has an unfalsifiable or misbehaving gate:\n"
        + result.stderr.strip() + "\n" + result.stdout.strip()[-1500:])
    assert "behaved as declared" in result.stdout
