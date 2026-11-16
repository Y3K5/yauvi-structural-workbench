"""Domain-aware membrane orientation through the Workbench: registry entry, run, report and viewer files.

Invented coordinates only; nothing is downloaded. Software checks only: these tests do not
show that any real protein's membrane pose is correct.
"""
from __future__ import annotations

import hashlib
import json
import math
import threading
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pytest

from yauvi_platform.structural_workbench.store import StructuralAnalysisStore, analysis_definitions, tool_readiness

TYPE = 'domain_membrane_orientation'
THREE = dict(L='LEU', S='SER', Y='TYR', G='GLY', D='ASP', N='ASN', T='THR', E='GLU', K='LYS', Q='GLN', A='ALA')


def invented_barrel(n=12):
    """Idealised 12-strand barrel along +z with long loops on +z and a polar domain below (same layout as memorient's tests)."""
    R = 4.8 / (2 * math.sin(math.pi / n)); ca, seq = [], []
    for i in range(60):
        a = math.radians(100 * i); ca.append([2.3 * math.cos(a) + 3, 2.3 * math.sin(a), -45 + 1.5 * i * .4]); seq.append('EKSQ'[i % 4])
    for k in range(n):
        th = 2 * math.pi * k / n; up = k % 2 == 0
        zs = [-12.8 + 3.2 * j for j in range(9)]
        for j, z in enumerate(zs if up else zs[::-1]):
            out = j % 2 == 0; r = R + (0.9 if out else -0.9)
            ca.append([r * math.cos(th), r * math.sin(th), z]); seq.append('Y' if out and abs(z) > 9 else 'L' if out else 'S')
        th2 = 2 * math.pi * (k + .5) / n
        if up:
            for m in range(10):
                h = 8 * math.sin(math.pi * (m + 1) / 11); ca.append([(R + 2) * math.cos(th2), (R + 2) * math.sin(th2), 14 + h]); seq.append('GDNT'[m % 4])
        else:
            for m in range(2):
                ca.append([(R + 1) * math.cos(th2), (R + 1) * math.sin(th2), -15 - m]); seq.append('GN'[m])
    return np.array(ca), ''.join(seq)


def soluble_bundle():
    ca, seq = [], ''
    for h in range(4):
        for i in range(20):
            a = math.radians(100 * i); ca.append([2.3 * math.cos(a) + 10 * (h % 2), 2.3 * math.sin(a) + 10 * (h // 2), 1.5 * i]); seq += 'EKAQ'[i % 4]
    return np.array(ca), seq


def pdb_bytes(ca, seq, chain='A') -> bytes:
    return ("\n".join(f"ATOM  {i + 1:5d}  CA  {THREE[a]} {chain}{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
                      for i, ((x, y, z), a) in enumerate(zip(ca, seq))) + "\nEND\n").encode()


def make_case(store, analysis_id, raw, parameters=None):
    store.create(analysis_id, analysis_type=TYPE, question='Invented software check', subject_id='invented')
    upload = store.begin_ingest(analysis_id, role='structure', file_name='invented.pdb', size=len(raw),
                                expected_sha256=hashlib.sha256(raw).hexdigest())
    store.ingest_chunk(upload['upload_id'], 0, raw); store.finalize_ingest(upload['upload_id'])
    if parameters is not None:
        store.update_parameters(analysis_id, parameters)


def test_registered_last_with_full_contract_and_ready():
    definitions = analysis_definitions()
    assert definitions[-1]['analysis_type'] == TYPE and sum(row['analysis_type'] == TYPE for row in definitions) == 1
    d = definitions[-1]
    assert d['readiness'] == 'experimental' and d['claim_ceiling'] and d['non_claim']
    assert [r['role'] for r in d['inputs']] == ['structure', 'validation_report'] and d['inputs'][0]['accepted_extensions'] == ['.pdb', '.ent']
    assert {p['name'] for p in d['parameters']} == {'chain', 'skip_n'}


def test_readiness_row(tmp_path):
    rows = [r for r in tool_readiness(tmp_path) if r['analysis_type'] == TYPE]
    assert len(rows) == 1 and rows[0]['state'] == 'experimental' and len(rows[0]['package_source_sha256']) == 64


def test_run_report_viewer_files_and_reuse(tmp_path):
    store = StructuralAnalysisStore(tmp_path / 'ws')
    make_case(store, 'barrel', pdb_bytes(*invented_barrel()), {'chain': 'A', 'skip_n': 0})
    pre = store.preflight('barrel')
    assert pre['valid'] and any(c['name'] == 'evidence:domain_orientation_input' and c['ok'] for c in pre['checks'])
    run = store.run('barrel')
    assert run['status'] == 'completed' and run['steps'] == [{'module_id': TYPE, 'exit_code': 0}]
    assert run['scientific_qualification'] == 'not_independently_qualified' and run['publication_authorized'] is False
    snap = store.snapshot('barrel')
    assert snap['oriented_structure']['file_name'] == 'ORIENTED_STRUCTURE.pdb'
    assert snap['membrane_layer']['file_name'] == 'MEMBRANE_LAYER.json'
    report = snap['report']
    assert [t['path'] for t in report['tables']] == ['outputs/memorient/RESIDUE_ZONES.tsv']
    assert 'independent reference orientation not compared in this run' in report['missing_evidence']
    html = store.artifact_path('barrel', run['run_id'], 'REPORT.html').read_text()
    assert 'RESIDUE_ZONES.tsv' in html
    layer = json.loads(store.artifact_path('barrel', run['run_id'], 'outputs/memorient/MEMBRANE_LAYER.json').read_text())
    assert layer['membrane_spanning'] is True and layer['membrane_slab']['half_thickness'] > 0
    assert store.run('barrel')['run_id'] == run['run_id']          # identical inputs: verified reuse, no new run


def test_no_membrane_core_is_scientifically_incomplete_without_viewer_files(tmp_path):
    store = StructuralAnalysisStore(tmp_path / 'ws')
    make_case(store, 'soluble', pdb_bytes(*soluble_bundle()))
    run = store.run('soluble')
    assert run['status'] == 'scientifically_incomplete'
    snap = store.snapshot('soluble')
    assert snap['oriented_structure'] is None and snap['membrane_layer'] is None


@pytest.mark.parametrize('parameters, message', [({'chain': 'ZZ'}, 'one character'), ({'chain': 'Z'}, 'No usable C-alpha'),
                                                 ({'skip_n': 5000}, 'whole number')])
def test_bad_settings_block_preflight(tmp_path, parameters, message):
    store = StructuralAnalysisStore(tmp_path / 'ws')
    make_case(store, 'bad', pdb_bytes(*invented_barrel()), parameters)
    pre = store.preflight('bad')
    check = next(c for c in pre['checks'] if c['name'] == 'evidence:domain_orientation_input')
    assert not pre['valid'] and not check['ok'] and message in check['detail']


def test_negative_skip_is_refused_when_parameters_are_saved(tmp_path):
    from yauvi_platform.structural_workbench.store import AnalysisError
    store = StructuralAnalysisStore(tmp_path / 'ws')
    with pytest.raises(AnalysisError, match='non-negative'):
        make_case(store, 'neg', pdb_bytes(*invented_barrel()), {'skip_n': -1})


def test_viewer_routes_serve_the_oriented_model_and_layer(tmp_path):
    from yauvi_structural_workbench.server import WorkbenchServer
    server = WorkbenchServer(('127.0.0.1', 0), tmp_path / 'ws')
    thread = threading.Thread(target=server.serve_forever, daemon=True); thread.start()
    base = f'http://127.0.0.1:{server.server_port}'
    try:
        make_case(server.store, 'barrel-http', pdb_bytes(*invented_barrel()))
        run = server.store.run('barrel-http')
        with urlopen(Request(f"{base}/api/analyses/barrel-http/oriented-structure/{run['run_id']}"), timeout=30) as r:
            assert b'MEMORIENT DOMAIN-AWARE ESTIMATE' in r.read()
        with urlopen(Request(f"{base}/api/analyses/barrel-http/membrane-layer/{run['run_id']}"), timeout=30) as r:
            assert json.loads(r.read())['context'] == 'domain_aware_estimate'
    finally:
        server.shutdown(); server.server_close(); thread.join(timeout=3)


def test_job_requires_experimental_acknowledgement(store, manager):
    store.add_case('invented', analysis_type=TYPE)
    with pytest.raises(Exception, match='experimental'):
        manager.submit('invented')


def test_interface_marks_the_type_experimental_and_draws_helical_slabs():
    app = (Path(__file__).resolve().parents[1] / 'src' / 'yauvi_structural_workbench' / 'ui' / 'app.js').read_text()
    assert app.count("'domain_membrane_orientation'") == 3
    assert "layer.membrane_spanning===true" in app
