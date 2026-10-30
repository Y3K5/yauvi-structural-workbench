"""Portable invented controls for distributions without private adapters."""
import hashlib
import json
import socket
import subprocess
from types import SimpleNamespace

import pytest

from yauvi_platform.structural_workbench import StructuralAnalysisStore
from yauvi_structural_workbench.capabilities import CapabilityView
from yauvi_structural_workbench.server import Handler
from yauvi_structural_workbench.build_info import ui_assets


@pytest.fixture
def core(tmp_path):
    store = StructuralAnalysisStore(tmp_path / 'invented-library')
    store._definitions = {'structure_qc': store._definitions['structure_qc']}
    store.source_root = tmp_path / 'invented-packages'
    store.installed_mode = False
    root = store.source_root / 'structqc/src/structqc'
    root.mkdir(parents=True)
    for name in ['core.py', 'cli.py']:
        (root / name).write_text("raise RuntimeError('Invented engine must never execute')\n")
    return store


def hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob('*') if p.is_file()}


def test_optional_methods_are_explicitly_unavailable(core):
    rows = CapabilityView(core).catalogue()['capabilities']
    assert len(rows) == 3
    assert rows[0]['source_state'] == 'configured'
    assert rows[0]['scientific_state'] == 'not_qualified_by_this_view'
    for row in rows[1:]:
        assert row['registration'] == 'unavailable'
        assert row['source_state'] == 'not_registered'
        assert 'identity_sha256' not in row


def test_inspection_neither_executes_nor_changes_invented_library(core, monkeypatch):
    before = hashes(core.workspace.parent)
    def forbidden(*args, **kwargs):
        raise AssertionError('Metadata inspection executed work')
    for target, name in [(subprocess, 'run'), (subprocess, 'Popen'),
                         (socket.socket, 'connect'), (socket, 'getaddrinfo'),
                         (core, 'run'), (core, 'preflight')]:
        monkeypatch.setattr(target, name, forbidden)
    view = CapabilityView(core)
    for row in view.catalogue()['capabilities']:
        view.detail(row['id'])
    assert hashes(core.workspace.parent) == before
    assert core.list() == []


def test_changed_source_requires_refresh(core):
    view = CapabilityView(core)
    old = view.detail('structure_qc')['identity_sha256']
    (core.source_root / 'structqc/src/structqc/core.py').write_text('# Invented changed source\n')
    current = view.detail('structure_qc', expected_sha256=old)
    assert current['source_state'] == 'changed_since_view'
    assert current['identity_sha256'] != old


def test_input_entries_do_not_claim_a_validation_pass(core):
    core.create('invented-case', analysis_type='structure_qc', question='Invented input accounting')
    path = core.workspace / 'invented.pdb'
    path.write_text('REMARK INVENTED SOFTWARE INPUT; NO BIOLOGICAL STRUCTURE\n')
    core.add_file('invented-case', role='structure', path=path)
    row = CapabilityView(core).detail('structure_qc', analysis_id='invented-case')
    assert row['input_state']['state'] == 'present_unchecked'
    assert row['input_state']['recorded_check']['state'] == 'not_recorded'


def request(core, path, origin=None):
    handler = object.__new__(Handler)
    handler.server = SimpleNamespace(store=core, allow_reference_fetch=False,
                                     server_port=8947, token='invented-token')
    handler.path = path
    handler.headers = {'Host': '127.0.0.1:8947'}
    if origin: handler.headers['Origin'] = origin
    handler._send = lambda code, body, *args, **kwargs: setattr(handler, 'response', (code, body))
    handler._dispatch(False)
    return handler.response


def test_api_inspection_stays_local_and_omits_absolute_paths(core):
    code, value = request(core, '/api/capabilities')
    assert code == 200
    assert value['privacy']['public_reference_fetch_enabled'] is False
    assert str(core.workspace.parent) not in json.dumps(value)
    assert request(core, '/api/capabilities', 'https://invented.invalid')[0] == 403


@pytest.mark.parametrize('path', [
    '/api/capabilities/unknown',
    '/api/capabilities/structure_qc?expected_sha256=invalid',
    '/api/capabilities/structure_qc?analysis_id=../../invented',
    '/api/capabilities?analysis_id=invented',
    '/api/capabilities/structure_qc?analysis_id=a&analysis_id=b',
])
def test_api_rejects_invalid_or_ambiguous_requests(core, path):
    assert request(core, path)[0] == 400


def test_native_host_is_served_and_included_in_build_identity():
    assets = ui_assets()
    assert b'/capabilities.js' in assets['index.html']
    assert b'id="capability-view"' in assets['index.html']
    assert b"/api/capabilities" in assets['capabilities.js']
