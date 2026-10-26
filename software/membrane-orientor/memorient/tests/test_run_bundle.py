"""Invented portable bundle and cache controls; no biological acquisitions."""
import hashlib
import json
from pathlib import Path
import shutil
import pytest
from memorient.run_bundle import write_manifest,load_manifest,artifact_path


def test_relocated_synthetic_bundle_retains_hashes(tmp_path):
    root=tmp_path/'original';root.mkdir();(root/'scene.json').write_text('{"synthetic":true}')
    path=write_manifest(root,stage='synthetic',artifacts={'scene':'scene.json'},sources=[],assumptions=['invented'])
    original=load_manifest(path)
    shutil.copytree(root,tmp_path/'moved')
    assert load_manifest(tmp_path/'moved/MEMBRANE_RUN.json')==original
    (tmp_path/'moved/scene.json').write_text('{}')
    with pytest.raises(ValueError,match='checksum'):load_manifest(tmp_path/'moved/MEMBRANE_RUN.json')


@pytest.mark.parametrize('name',['../secret','/absolute','C:/secret','a\\b','bad\nheader','bad\rheader','bad\"name','\x00bad'])
def test_paths_cannot_escape_bundle(tmp_path,name):
    with pytest.raises(ValueError):artifact_path(tmp_path,name)


def test_minimum_cannot_be_claimed_by_stage_alone(tmp_path):
    with pytest.raises(ValueError,match='acceptance'):
        write_manifest(tmp_path,stage='minimized',artifacts={})
    (tmp_path/'MEMBRANE_RUN.json').write_text(json.dumps({'schema_version':'memorient.run.v1','stage':'prepared',
        'units':'nm','artifacts':{},'equilibration':'completed'}))
    with pytest.raises(ValueError,match='Equilibration'):load_manifest(tmp_path/'MEMBRANE_RUN.json')


def test_metadata_survives_empty_and_partial_cache(tmp_path,monkeypatch):
    from memorient import starter_catalogue as module
    home=tmp_path/'catalogue';home.mkdir();cache=tmp_path/'cache';cache.mkdir()
    raw=b'REMARK invented cache fixture\n';sha=hashlib.sha256(raw).hexdigest()
    (home/'metadata.json').write_text(json.dumps({'organisms':[],'sources':[],
        'proteins':[{'id':'invented','view_pdb_sha256':sha},{'id':'absent','view_pdb_sha256':'0'*64}]}))
    monkeypatch.setattr(module,'HOME',home)
    module.configure(cache)
    try:
        assert len(module.catalogue_index()['proteins'])==2
        with pytest.raises(ValueError,match='unavailable'):module.protein_record('invented')
        (cache/(sha+'.pdb')).write_bytes(raw);module.load_catalogue.cache_clear()
        assert module.protein_record('invented')['pdb']==raw.decode()
        (cache/(sha+'.pdb')).write_text('changed');module.load_catalogue.cache_clear()
        assert module.catalogue_index()['proteins'][0]['availability']=='checksum_mismatch'
    finally:module.configure()


@pytest.mark.parametrize('charge',[float('nan'),float('inf'),None,True,'0'])
def test_minimized_charge_must_be_finite_and_numeric(charge):
    from memorient.run_bundle import validate_minimum
    record={'minimization':{'converged':True,'maximum_force_kj_mol_nm':10.,'potential_energy_kj_mol':-1.},
            'charge_e':charge,'geometry_checks':{'passed':True}}
    with pytest.raises(ValueError,match='acceptance'):validate_minimum(record)
