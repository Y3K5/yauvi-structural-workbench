"""Invented public-distribution controls for record listing and QC custody."""
import hashlib
from importlib.resources import files
from pathlib import Path

import pytest

from yauvi_platform.structural_workbench import AnalysisError, StructuralAnalysisStore
from yauvi_platform.structural_workbench.biological_case import BiologicalCaseStore
from yauvi_structural_workbench.example import create_example


def custody(root):
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            if p.is_file() else None for p in root.rglob('*')}


@pytest.mark.parametrize('with_revisions', [False, True])
def test_unwritten_folders_do_not_block_committed_analysis_listing(tmp_path, with_revisions):
    store = StructuralAnalysisStore(tmp_path / 'invented-library')
    store.create('saved', analysis_type='structure_qc', question='Invented saved case')
    empty = store.cases_root / 'interrupted'
    empty.mkdir()
    if with_revisions:
        (empty / 'revisions').mkdir()
    before = custody(store.workspace)
    assert [row['analysis_id'] for row in store.list()] == ['saved']
    with pytest.raises(AnalysisError):
        store.load('interrupted')
    assert custody(store.workspace) == before


def test_biological_listing_ignores_empty_folders_but_keeps_invalid_record_visible(tmp_path):
    store = BiologicalCaseStore(tmp_path / 'invented-library')
    (store.root / 'interrupted' / 'revisions').mkdir(parents=True)
    damaged = store.root / 'damaged'
    damaged.mkdir()
    (damaged / 'case.json').write_text('{')
    before = custody(store.root)
    rows = store.list()
    assert len(rows) == 1 and rows[0]['id'] == 'damaged'
    assert rows[0]['state'] == 'invalid' and rows[0]['error']
    with pytest.raises(ValueError):
        store.load('damaged')
    assert custody(store.root) == before


def test_existing_qc_source_and_bundled_reference_have_verified_local_bytes():
    root = Path(__file__).resolve().parents[3]
    source = root / 'software/structqc/examples/reference.fasta'
    marker = root / 'software/structqc/src/structqc/py.typed'
    for path in [source, marker]:
        assert not (getattr(path.stat(), 'st_flags', 0) & 0x40000000), 'Local fixture is offloaded'
    expected = 'daf3cb9543fe62dfb2104efcf1ea54cd710eedfdbf1b5c3b1001bc6afd06cc26'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == expected
    assert marker.read_bytes() == b'\n'
    bundled = files('yauvi_structural_workbench').joinpath('examples/reference.fasta')
    assert hashlib.sha256(bundled.read_bytes()).hexdigest() == expected


@pytest.fixture
def completed(tmp_path):
    workspace = tmp_path / 'invented-library'
    create_example(workspace, 'invented-qc')
    store = StructuralAnalysisStore(workspace)
    result = store.run('invented-qc')
    assert result['status'] == 'completed'
    return store, result['run_id']


def test_old_qc_output_retains_original_revision_after_parameter_edit(completed):
    store, run_id = completed
    original = store.snapshot('invented-qc', run_id=run_id)
    directory = store.artifact_path('invented-qc', run_id, 'REPORT_DATA.json').parent
    before = custody(directory)
    store.update_parameters('invented-qc', {'chain': 'A'})
    reopened = store.snapshot('invented-qc', run_id=run_id)
    assert not reopened['report_matches_case']
    assert reopened['run']['input_revision_sha256'] == original['run']['input_revision_sha256']
    assert reopened['run_inputs'] == original['run_inputs']
    assert custody(directory) == before


def test_unrelated_case_cannot_claim_qc_run_or_input(completed):
    store, run_id = completed
    source = store.load('invented-qc')
    store.create('unrelated', analysis_type='structure_qc', question='Invented unrelated case')
    before = custody(store.workspace)
    with pytest.raises(AnalysisError, match='does not belong'):
        store.snapshot('unrelated', run_id=run_id)
    with pytest.raises(AnalysisError, match='not attached'):
        store.input_path('unrelated', source['inputs'][0]['sha256'])
    assert custody(store.workspace) == before
