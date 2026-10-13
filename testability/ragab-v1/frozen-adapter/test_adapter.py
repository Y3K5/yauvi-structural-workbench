from pathlib import Path
import copy, os
import numpy as np
import pytest
from adapter import *

@pytest.fixture
def reference_arrays():
    p=ROOT/'prepared'/'REF_OPEN'/'state-input'
    keys=read(p/'MAPPING.json')['keys']
    arrays=verify_mapping([p/'query.cif',p/'open.cif',p/'closed.cif'],keys)
    return arrays,keys

def test_reference_self_and_independent_quaternion(reference_arrays):
    arrays,_=reference_arrays
    for a in arrays:
        assert state._superpose(a,a)[1]<NUMERICAL_TOLERANCE_A
    assert abs(state._superpose(arrays[1],arrays[2])[1]-horn_rmsd(arrays[1],arrays[2]))<NUMERICAL_TOLERANCE_A

def test_rigid_transform_invariance(reference_arrays):
    arrays,_=reference_arrays; a,b=arrays[1:]
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]],float)
    assert abs(state._superpose(a@rotation+[123,19,-46],b)[1]-state._superpose(a,b)[1])<NUMERICAL_TOLERANCE_A

def test_relative_lid_motion_retained(reference_arrays):
    arrays,keys=reference_arrays; a=arrays[0].copy(); moved=a.copy()
    moved[np.array([k[0]=='B' for k in keys])]+=[8,3,-1]
    assert state._superpose(moved,a)[1]>0.5

def test_pair_reordering_and_renaming_by_identity(tmp_path):
    cfg=read(ROOT/'CONFIG.json')['protomers']['REF_OPEN']; source=ROOT/'sources'/'6SML.cif'
    original,meta=pair(source,cfg)
    # Reorder and rename the complete source while updating the declared author IDs.
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from Bio.PDB import MMCIFIO
    d=MMCIF2Dict(str(source)); d['_atom_site.auth_asym_id']=[{'A':'Y','B':'X'}.get(v,v) for v in d['_atom_site.auth_asym_id']]
    idx=list(reversed(range(len(d['_atom_site.id']))))
    for key,value in list(d.items()):
        if key.startswith('_atom_site.'): d[key]=[value[i] for i in idx]
    path=tmp_path/'renamed.cif'; writer=MMCIFIO(); writer.set_dict(d); writer.save(str(path))
    changed={**cfg,'ragA_auth':'X','ragB_auth':'Y','candidate_ragA_auth':['X']}
    other,_=pair(path,changed)
    a,b=residue_dict(original),residue_dict(other)
    assert a.keys()==b.keys()
    for k in a:
        for atom in a[k]: np.testing.assert_allclose(atom.coord,b[k][atom.name].coord,atol=1e-6)

def test_wrong_pair_rejected():
    cfg=read(ROOT/'CONFIG.json')['protomers']['Q1']
    with pytest.raises(ValueError,match='pairing contradicts'):
        pair(ROOT/'sources'/'6SMQ.cif',{**cfg,'ragA_auth':'B'})

def test_wrong_polymer_rejected():
    cfg=read(ROOT/'CONFIG.json')['protomers']['Q1']
    with pytest.raises(ValueError,match='polymer identities'):
        pair(ROOT/'sources'/'6SMQ.cif',{**cfg,'ragA_auth':'D'})

def test_missing_mapped_coordinate_rejected(tmp_path):
    p=ROOT/'prepared'/'REF_OPEN'/'state-input'; s=load(p/'query.cif')
    s['A'].detach_child(next(iter(s['A'])).id)
    path=tmp_path/'missing.cif'; save_structure(path,s)
    with pytest.raises(ValueError,match='mapped coordinates'):
        verify_mapping([path,p/'open.cif'],read(p/'MAPPING.json')['keys'])

def test_checksum_mismatch_fails_engine():
    p=ROOT/'prepared'/'REF_OPEN'/'state-input'; evidence=read(p/'query-qc'/'STRUCTURE_EVIDENCE.json')
    evidence['coordinate']['sha256']='0'*64
    with pytest.raises(state.InputError,match='checksum'):
        state.analyze(evidence,read(p/'references.json'),reference_base=p,structure_path=p/'query.cif')

def test_missing_freesasa_fails_before_fallback(monkeypatch):
    monkeypatch.setenv('PATH','')
    with pytest.raises(RuntimeError,match='fallback forbidden'): ensure_freesasa()

def test_signature_positions_excluded(reference_arrays):
    _,keys=reference_arrays
    assert not any(k[0]=='A' and (103<=k[1]<=114 or 838<=k[1]<=841) for k in keys)

def test_no_nonprotein_or_hydrogen_in_pairs():
    for q in ['REF_OPEN','REF_CLOSED','Q1','Q2']:
        m=load(ROOT/'prepared'/q/'pair.cif')
        assert [c.id for c in m]==['A','B']
        assert all(r.id[0]==' ' for r in m.get_residues())
        assert all(a.element not in ('H','D') for a in m.get_atoms())

def test_reference_engine_aliases():
    for q,expected in [('REF_OPEN','active_like'),('REF_CLOSED','inactive_like')]:
        p=ROOT/'prepared'/q/'state-input'
        evidence=read(p/'query-qc'/'STRUCTURE_EVIDENCE.json')
        result=state.analyze(evidence,read(p/'references.json'),reference_base=p,structure_path=p/'query.cif')
        assert result['frame_metrics'][0]['call']==expected
