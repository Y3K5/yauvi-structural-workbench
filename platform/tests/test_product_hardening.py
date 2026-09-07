"""Adversarial product contracts; scientific references are tested separately."""
from __future__ import annotations
import hashlib
import importlib.util
import json
import threading
import time
from pathlib import Path
import pytest
from yauvi_platform.structural_workbench import StructuralAnalysisStore, AnalysisError
from yauvi_structural_workbench.jobs import JobManager
from yauvi_structural_workbench.server import Handler

ROOT=Path(__file__).resolve().parents[2]

def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'tools'/f'{name}.py')
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj

def case(tmp_path):
    store=StructuralAnalysisStore(tmp_path)
    store.create('qc',analysis_type='structure_qc',question='What evidence is available?')
    for role,name in [('structure','model.pdb'),('provenance','provenance.json'),('validation_report','validation.json')]:
        store.add_file('qc',role=role,path=ROOT/'structqc/examples'/name)
    return store

def test_case_revisions_are_preserved(tmp_path):
    store=StructuralAnalysisStore(tmp_path)
    initial=store.create('qc',analysis_type='structure_qc',question='A research question')
    store.update_parameters('qc',{})
    archive=store.cases_root/'qc/revisions'/f"{initial['revision_sha256']}.json"
    assert json.loads(archive.read_text())==initial
    assert len(list(archive.parent.glob('*.json')))==2

def test_failed_attempt_is_not_a_cache_hit_or_overwritten(tmp_path,monkeypatch):
    store=case(tmp_path)
    monkeypatch.setattr(store,'_run_registered',lambda *a,**k: ([],2))
    first=store.run('qc');path=store.cases_root/'qc/runs'/first['run_id']/'ANALYSIS_RUN.json'
    original=path.read_bytes()
    monkeypatch.setattr(store,'_run_registered',lambda *a,**k: ([],0))
    second=store.run('qc')
    assert first['status']=='failed' and second['status']=='completed'
    assert first['run_id']!=second['run_id'] and path.read_bytes()==original
    assert store.run('qc')['run_id']==second['run_id']

def test_missing_completed_artifact_refuses_cache(tmp_path,monkeypatch):
    store=case(tmp_path);monkeypatch.setattr(store,'_run_registered',lambda *a,**k:([],0))
    result=store.run('qc')
    (store.cases_root/'qc/runs'/result['run_id']/'REPORT_DATA.json').write_text('corrupted')
    with pytest.raises(AnalysisError,match='checksum'):store.run('qc')

def test_unexpected_engine_error_records_failed_attempt(tmp_path,monkeypatch):
    store=case(tmp_path)
    def fail(*a,**k):raise OSError('simulated unavailable executable')
    monkeypatch.setattr(store,'_run_registered',fail)
    result=store.run('qc')
    assert result['status']=='failed' and result['exit_code']==2

def test_portable_projection_preserves_original_and_accounts_for_changes(tmp_path):
    exporter=module('export_execution_evidence');source=tmp_path/'source';source.mkdir()
    doc={'cases':[{'record_id':'case-a','output_dir':'/Users/researcher/work/case-a','passed':True}]}
    path=source/'EXECUTION_STATUS.json';path.write_text(json.dumps(doc));before=path.read_bytes()
    report=exporter.export(source,tmp_path/'portable')
    assert path.read_bytes()==before
    assert report['files'][0]['replacements'][0]['replacement']=='case-a'
    assert b'/Users/' not in (tmp_path/'portable/EXECUTION_STATUS.json').read_bytes()
    assert report['publication_authorized'] is False

def test_projection_blocks_unhandled_paths_without_partial_output(tmp_path):
    exporter=module('export_execution_evidence');source=tmp_path/'source';source.mkdir()
    (source/'EXECUTION_STATUS.json').write_text(json.dumps({'unexpected':'/Users/researcher/private'}))
    with pytest.raises(ValueError,match='unresolved'):exporter.export(source,tmp_path/'portable')
    assert not (tmp_path/'portable').exists()

def test_guard_rejects_foreign_host_origin_and_missing_session():
    from types import SimpleNamespace
    handler=object.__new__(Handler);handler.server=SimpleNamespace(server_port=8947,token='test-token')
    for headers in ({'Host':'evil.example:8947'}, {'Host':'127.0.0.1:8947','Origin':'https://evil.example'},
                    {'Host':'127.0.0.1:8947','Sec-Fetch-Site':'cross-site'}, {'Host':'127.0.0.1:8947'}):
        handler.headers=headers
        with pytest.raises(PermissionError):handler._guard(True)
    handler.headers={'Host':'127.0.0.1:8947','Origin':'http://127.0.0.1:8947','X-Yauvi-Token':'test-token'}
    handler._guard(True)

def test_jobs_cancel_persist_and_recover_interrupted(tmp_path,monkeypatch):
    store=StructuralAnalysisStore(tmp_path);store.create('qc',analysis_type='structure_qc',question='Test')
    started=threading.Event()
    def run(aid,*,cancel_event,**callbacks):
        started.set();assert cancel_event.wait(3)
        return {'status':'cancelled','run_id':'run-'+'a'*16}
    monkeypatch.setattr(store,'run',run)
    jobs=JobManager(store)
    try:
        row=jobs.submit('qc');assert started.wait(3)
        jobs.cancel(row['job_id']);jobs._queue.join()
        assert jobs.list()[0]['state']=='cancelled'
        with pytest.raises(AnalysisError,match='already owns'):JobManager(store)
    finally:jobs.close()
    path=store.root/'jobs'/f"{row['job_id']}.json";record=json.loads(path.read_text());record['state']='running';path.write_text(json.dumps(record))
    restored=JobManager(store)
    try:assert restored.list()[0]['state']=='interrupted'
    finally:restored.close()

def test_experimental_job_requires_explicit_acknowledgment(tmp_path):
    store=StructuralAnalysisStore(tmp_path);store.create('membrane',analysis_type='membrane_orientation',question='Explore')
    jobs=JobManager(store)
    try:
        with pytest.raises(AnalysisError,match='experimental'):jobs.submit('membrane')
    finally:jobs.close()


def test_completed_cache_requires_checksum_coverage(tmp_path,monkeypatch):
    store=case(tmp_path);monkeypatch.setattr(store,'_run_registered',lambda *a,**k:([],0))
    result=store.run('qc')
    path=store.cases_root/'qc/runs'/result['run_id']/'CHECKSUMS.json'
    path.write_text(json.dumps({'files':{}}))
    with pytest.raises(AnalysisError,match='coverage'):store.run('qc')


def test_archive_corruption_cannot_satisfy_cache(tmp_path,monkeypatch):
    store=case(tmp_path);monkeypatch.setattr(store,'_run_registered',lambda *a,**k:([],0))
    result=store.run('qc')
    (store.cases_root/'qc/runs'/result['run_id']/'RAW_EVIDENCE.zip').write_bytes(b'broken archive')
    with pytest.raises(AnalysisError,match='checksum'):store.run('qc')


def test_execution_limit_terminates_owned_process_group(tmp_path):
    import sys
    store=StructuralAnalysisStore(tmp_path)
    with pytest.raises(TimeoutError,match='execution limit'):
        store._execute([sys.executable,'-c','import time; time.sleep(30)'],cwd=tmp_path,
                       log_path=tmp_path/'log.txt',timeout_seconds=.05)
    assert (tmp_path/'log.txt').exists()


def test_body_rejects_oversize_partial_and_malformed_requests():
    from email.message import Message
    from io import BytesIO
    from types import SimpleNamespace
    for size,body in [('1048577',b'{}'),('5',b'{}'),('2',b'[]'),('1',b'{')]:
        handler=object.__new__(Handler);headers=Message();headers['Content-Length']=size;headers['Content-Type']='application/json'
        handler.headers=headers;handler.rfile=BytesIO(body);handler.connection=SimpleNamespace(settimeout=lambda n:None)
        with pytest.raises(ValueError):handler._body()


def test_report_assembly_failure_cannot_leave_completed_run(tmp_path,monkeypatch):
    store=case(tmp_path);monkeypatch.setattr(store,'_run_registered',lambda *a,**k:([],0))
    def fail(manifest,record,run_dir):
        (run_dir/'REPORT_DATA.json').write_text('{"status":"completed"}')
        raise OSError('simulated report write failure')
    monkeypatch.setattr(store,'_render_report',fail)
    record=store.run('qc');root=store.cases_root/'qc/runs'/record['run_id']
    assert record['status']=='failed' and record['exit_code']==2
    assert (root/'PARTIAL_REPORT_DATA.json').is_file()
    assert store.snapshot('qc')['report'] is None
    assert json.loads((root/'ANALYSIS_RUN.json').read_text())['status']=='failed'


def test_live_progress_is_durable_and_snapshots_do_not_share_events(tmp_path, monkeypatch):
    from types import SimpleNamespace
    store = case(tmp_path)
    def run(aid, *, cancel_event, on_process, on_progress):
        on_progress('checking', 'Checking inputs')
        on_process(SimpleNamespace(args=['python', '-m', 'structqc.cli', '/private/input']))
        on_progress('report', 'Assembling evidence')
        return {'status': 'completed', 'run_id': 'run-' + 'a' * 16}
    monkeypatch.setattr(store, 'run', run)
    jobs = JobManager(store)
    try:
        job = jobs.submit('qc'); jobs._queue.join()
        row = jobs.list()[0]
        assert row['state'] == 'completed'
        assert row['created_at'] <= row['started_at'] <= row['finished_at']
        assert [e['stage'] for e in row['events']] == ['queued', 'starting', 'checking', 'engine', 'report', 'completed']
        assert '/private/input' not in json.dumps(row)
        saved = json.loads((store.root/'jobs'/f"{job['job_id']}.json").read_text())
        assert saved['events'] == row['events']
        row['events'].clear()
        assert jobs.list()[0]['events'] == saved['events']
    finally:
        jobs.close()


def test_history_keeps_run_inputs_separate_from_editable_case(tmp_path, monkeypatch):
    store = case(tmp_path)
    monkeypatch.setattr(store, '_run_registered', lambda *a, **k: ([], 0))
    first = store.run('qc')
    assert store.snapshot('qc')['report_matches_case']
    store.update_parameters('qc', {'chain': 'A'})
    second = store.run('qc')
    old = store.snapshot('qc', run_id=first['run_id'])
    assert old['analysis']['parameters'] == {'chain': 'A'}
    assert old['run_inputs']['parameters'] == {}
    assert old['run']['run_id'] == first['run_id']
    assert old['report_matches_case'] is False
    assert old['analysis']['latest_run_id'] == second['run_id']
    assert len(old['available_artifacts']) == 5
    with pytest.raises(AnalysisError, match='does not belong'):
        store.snapshot('qc', run_id='run-'+'f'*16)


def test_progress_reports_verified_cache_and_report_assembly(tmp_path, monkeypatch):
    store = case(tmp_path)
    monkeypatch.setattr(store, '_run_registered', lambda *a, **k: ([], 0))
    events = []
    record = store.run('qc', on_progress=lambda stage, message: events.append(stage))
    assert events == ['checking', 'identity', 'analysis', 'report']
    events.clear()
    assert store.run('qc', on_progress=lambda stage, message: events.append(stage)) == record
    assert events == ['checking', 'identity', 'reused']


def test_every_json_template_download_preserves_its_bytes():
    from io import BytesIO
    from yauvi_platform.structural_workbench import analysis_definitions, template_artifact
    templates = {role['template_id'] for definition in analysis_definitions()
                 for role in definition['inputs'] if role.get('template_id')}
    for template in templates:
        filename, kind, content = template_artifact(template)
        handler = object.__new__(Handler)
        handler.wfile = BytesIO()
        headers = {}
        handler.send_response = lambda status: headers.update(status=status)
        handler.send_header = lambda key, value: headers.update({key: value})
        handler.end_headers = lambda: None
        handler._send(200, content, kind, attachment=filename)
        assert handler.wfile.getvalue() == content
        assert json.loads(content)
        assert headers['Content-Type'] == 'application/json'
        assert headers['Content-Length'] == str(len(content))
        assert filename in headers['Content-Disposition']


def test_projection_detects_current_account_without_hardcoding_it(monkeypatch):
    exporter = module('export_execution_evidence')
    monkeypatch.setattr(exporter.Path, 'home', lambda: Path('/home/example-account'))
    with pytest.raises(ValueError, match='unresolved'):
        exporter.project({'unexpected_account': 'example-account'})
