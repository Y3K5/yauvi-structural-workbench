"""Invented public-profile route and native-format custody controls."""
import hashlib
from io import BytesIO
import json
from types import SimpleNamespace
import zipfile

import pytest

from yauvi_platform.structural_workbench.function_investigations import FunctionInvestigations
from yauvi_platform.structural_workbench.function_records import record
from yauvi_platform.structural_workbench.store import AnalysisError
from yauvi_structural_workbench.server import Handler
from yauvi_structural_workbench.function_record_routes import get
from yauvi_structural_workbench.build_info import ui_assets


@pytest.fixture
def local(tmp_path, monkeypatch):
    monkeypatch.delenv('YAUVI_PROTEIN_FUNCTION_EVIDENCE', raising=False)
    key = 'example|EXAMPLE'
    row = {'key': key, 'id': 'EXAMPLE', 'species': 'example', 'strain': 'Invented reference',
           'length': 20, 'sequence_sha256': '0' * 64, 'name': 'Invented unresolved protein', 'primary_group': 'unresolved'}
    data = {'schema': 'functional-protein-catalog-v1', 'dataset_id': 'a' * 64, 'source_ledger_sha256': 'b' * 64,
            'records': [row], 'groups': [{'id': 'unresolved', 'label': 'Unresolved'}],
            'species': [{'id': 'example', 'label': 'Invented organism'}],
            'evidence': {key: {'location': {'status': 'unresolved', 'zone': 'unplaced'},
                               'structure': {'status': 'not_joined', 'available': False}}}}
    path = tmp_path / 'catalogue.json'
    path.write_text(json.dumps(data))
    config = tmp_path / 'config.json'
    config.write_text(json.dumps({'schema': 'function-research-config-v1', 'catalogue_path': str(path),
                                 'catalogue_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}))
    monkeypatch.setenv('YAUVI_FUNCTION_RESEARCH_CONFIG', str(config))
    service = FunctionInvestigations(tmp_path / 'workspace')
    value = record(key)
    request = {'key': key, 'record_sha256': value['record_sha256'], 'question': 'Invented evidence question',
               'action_id': 'domains', 'hypotheses': [{'statement': 'Invented proposed explanation',
                'alternative': 'Invented alternative', 'distinguishing_evidence': 'Review exact source evidence',
                'assertion_ids': []}], 'notes': 'Invented route check', 'origin_study': None}
    return service, request


def handler(local, path):
    h = object.__new__(Handler)
    h.path = path
    h.command = 'POST'
    h.server = SimpleNamespace(store=local[0].studies.analysis, ui_assets=ui_assets(), jobs=None)
    h.sent = []
    h._send = lambda *args, **kw: h.sent.append((args, kw))
    return h


def test_post_consumes_parsed_body_once_and_put_cannot_save(local):
    h = handler(local, '/api/protein-functions/investigations')
    reads = []
    h._body = lambda: reads.append(True) or local[1]
    h._mutate(['api', 'protein-functions', 'investigations'])
    assert len(reads) == 1 and h.sent[0][0][0] == 201
    assert len(local[0].list(local[1]['key'])) == 1
    h.command = 'PUT'
    with pytest.raises(ValueError):
        h._mutate(['api', 'protein-functions', 'investigations'])
    assert len(local[0].list(local[1]['key'])) == 1


@pytest.mark.parametrize('query', ['q=a&q=b', 'gap=outside', 'unexpected=x'])
def test_invalid_get_filters_refused(local, query):
    h = handler(local, '/api/protein-functions/catalogue?' + query)
    with pytest.raises((ValueError, AnalysisError)):
        get(h, ['api', 'protein-functions', 'catalogue'])


@pytest.mark.parametrize('headers', [
    {'Host': 'outside.invalid:8968'},
    {'Host': '127.0.0.1:8968', 'Origin': 'http://outside.invalid'},
    {'Host': '127.0.0.1:8968', 'Sec-Fetch-Site': 'cross-site'},
    {'Host': '127.0.0.1:8968'},
])
def test_investigation_uses_native_host_origin_and_token_guards(local, headers):
    h = handler(local, '/api/protein-functions/investigations')
    h.headers = headers
    h.server.server_port = 8968
    h.server.token = 'invented-local-session'
    with pytest.raises(PermissionError):
        h._guard(mutating=True)
    assert local[0].list(local[1]['key']) == []


def test_native_restart_revision_and_bundle_routes_keep_exact_proposals(local):
    service, request = local
    old = service.save(request)
    later = service.save({**request, 'study_id': old['study_id'], 'expected_revision_sha256': old['revision_sha256'],
                          'question': 'Revised invented question'})
    reopened = FunctionInvestigations(service.studies.analysis.workspace)
    assert reopened.load(old['study_id'], old['revision_sha256'])['study'] == old
    assert reopened.load(old['study_id'])['study'] == later
    h = handler(local, '/api/studies/' + old['study_id'] + '/bundle?revision=' + old['revision_sha256'])
    assert get(h, ['api', 'studies', old['study_id'], 'bundle'])
    z = zipfile.ZipFile(BytesIO(h.sent[0][0][1]))
    assert json.loads(z.read('STUDY_RECORD.json')) == old
    for name, digest in json.loads(z.read('CHECKSUMS.json'))['files'].items():
        assert hashlib.sha256(z.read(name)).hexdigest() == digest
    assert not old['evidence'] and not old['scientifically_approved']


def test_bound_origin_requires_exact_native_subject(local):
    service, request = local
    origin = service.studies.create('invented-origin', title='Invented origin', question='Invented question', scope='Software fixture')
    with pytest.raises(AnalysisError):
        service.save({**request, 'origin_study': {'study_id': origin['study_id'], 'revision_sha256': origin['revision_sha256']}})
    origin['subjects'].append({'subject_id': 'protein-reference', 'kind': 'sequence', 'label': 'Invented protein',
                               'reference_key': request['key'], 'sequence_sha256': '0' * 64, 'artifact_refs': []})
    with service.studies._transaction():
        origin = service.studies._save(origin)
    bound = service.save({**request, 'origin_study': {'study_id': origin['study_id'], 'revision_sha256': origin['revision_sha256']}})
    assert service.load(bound['study_id'])['packet']['plan']['origin_study']['revision_sha256'] == origin['revision_sha256']


def test_public_custody_refuses_general_case_link_and_noninvestigation_export(local):
    service, _ = local
    origin = service.studies.create('invented-origin', title='Invented origin', question='Invented question', scope='Software fixture')
    with pytest.raises(AnalysisError):
        service.studies.bundle(origin['study_id'])
    origin['case_links'] = [{'analysis_id': 'not-admitted'}]
    with service.studies._transaction():
        service.studies._save(origin)
    with pytest.raises(AnalysisError):
        service.studies.snapshot(origin['study_id'])


def test_public_assets_have_local_links_and_use_shared_renderer(local):
    assets = ui_assets()
    assert b'/function-records/' in assets['index.html']
    assert b'/protein-functions.js' in assets['function-records.html']
    assert b'WorkbenchProteinFunctions' in assets['function-records-host.js']
    h = handler(local, '/function-records/')
    h._get(['function-records'])
    assert h.sent[0][0] == (200, assets['function-records.html'], 'text/html')
