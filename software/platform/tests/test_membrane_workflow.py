"""Actual synthetic Workbench custody/execution/import, without MD or network."""
import json
from pathlib import Path
import zipfile
import pytest
from yauvi_platform.structural_workbench.store import StructuralAnalysisStore
from yauvi_platform.structural_workbench.membrane_workflow import recorded_manifest,unpack,runtime

PDB='ATOM      1  CA  ALA A   1       0.000   0.000  -1.000  1.00 20.00           C\nATOM      2  CA  ALA A   2       0.000   0.000   1.000  1.00 20.00           C\nEND\n'


@pytest.mark.parametrize('kind,size',[('plane',[100,100]),('sphere',[100]),('cylinder',[100,180]),('ellipsoid',[100,120,150])])
def test_actual_placement_and_portable_import(tmp_path,kind,size):
    store=StructuralAnalysisStore(tmp_path/'workspace');store.create('surface',analysis_type='membrane_surface_placement',question='Invented geometry check')
    source=tmp_path/'source.pdb';source.write_text(PDB)
    spec=tmp_path/'spec.json';spec.write_text(json.dumps({'synthetic':True,'surface':{'kind':kind,'size':size,'thickness':30},'placement':{'axis':[0,0,1],'u':0,'v':0,'spin':0,'side':1,'depth':0}}))
    store.add_file('surface',role='structure',path=source);store.add_file('surface',role='placement_spec',path=spec)
    assert store.preflight('surface')['valid']
    result=store.run('surface')
    assert result['status']=='completed',result
    path,doc=recorded_manifest(store,'surface',result['run_id'])
    assert doc['stage']=='synthetic' and doc['membrane_geometry']['kind']==kind
    assert doc['artifacts']['source']['sha256']==next(r['sha256'] for r in store.load('surface')['inputs'] if r['role']=='structure')
    archive=tmp_path/'bundle.zip'
    with zipfile.ZipFile(archive,'w') as out:
        out.write(path,'MEMBRANE_RUN.json')
        for ref in doc['artifacts'].values():out.write(path.parent/ref['path'],ref['path'])
    store.create('imported',analysis_type='membrane_result_import',question='Reopen invented result')
    store.add_file('imported',role='membrane_bundle',path=archive)
    assert store.preflight('imported')['valid']
    imported=store.run('imported');assert imported['status']=='completed',imported
    _,after=recorded_manifest(store,'imported',imported['run_id']);assert after==doc
    (path.parent/'scene.json').write_text('{}')
    with pytest.raises(ValueError,match='checksum'):recorded_manifest(store,'surface',result['run_id'])


def test_missing_runtime_is_held_without_installation(monkeypatch):
    monkeypatch.setenv('YAUVI_MEMORIENT_GMX','nonexistent-memorient-binary')
    with pytest.raises(ValueError,match='held'):runtime()


def test_zip_traversal_and_unlisted_payload_are_held(tmp_path):
    archive=tmp_path/'bad.zip'
    with zipfile.ZipFile(archive,'w') as out:
        out.writestr('MEMBRANE_RUN.json',json.dumps({'artifacts':{'evil':{'path':'../escape'}}}))
        out.writestr('../escape','do not extract')
    with pytest.raises(ValueError,match='inside'):unpack(archive,tmp_path/'extract')
    assert not (tmp_path/'escape').exists()
