"""Adversarial identity/assembly/residue tests, with synthetic deposited records."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from yauvi_platform.structural_workbench.biological_case import BiologicalCase, BiologicalCaseStore


def digest(value):
    return hashlib.sha256(value).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value))


@pytest.fixture
def case_dir(tmp_path):
    directory = tmp_path / 'biological-cases' / 'synthetic'
    directory.mkdir(parents=True)
    columns = ['id', 'group_PDB', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id',
               'label_asym_id', 'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x',
               'Cartn_y', 'Cartn_z', 'occupancy', 'auth_seq_id', 'auth_asym_id', 'pdbx_PDB_model_num']
    cif = '''data_TEST
_entry.id TEST
_entity_poly.entity_id 1
_entity_poly.type polypeptide(L)
_entity_poly.pdbx_seq_one_letter_code_can AGST
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
_pdbx_poly_seq_scheme.auth_seq_num
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.pdb_ins_code
_pdbx_poly_seq_scheme.pdb_strand_id
A 1 1 ALA 10 10 . Z
A 1 2 GLY 10 10 A Z
A 1 3 SER ? 11 . Z
A 1 4 THR 12 12 . Z
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.method_details
_pdbx_struct_assembly.oligomeric_details
_pdbx_struct_assembly.oligomeric_count
1 author_defined_assembly ? dimeric 2
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(1,2)' A
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.vector[3]
1 'identity operation' 1 0 0 0 1 0 0 0 1 0 0 0
2 'point symmetry operation' -1 0 0 0 -1 0 0 0 1 20 0 0
loop_
'''
    cif += ''.join('_atom_site.' + c + '\n' for c in columns)
    cif += '''1 ATOM C CA . ALA A 1 1 . 1 2 3 1 10 Z 1
2 ATOM C CA A GLY A 1 2 A 2 3 4 0.6 10 Z 1
3 ATOM C CA B GLY A 1 2 A 5 6 7 0.4 10 Z 1
4 ATOM C CA . THR A 1 4 . 3 4 5 1 12 Z 1
5 ATOM C CA . ALA A 1 1 . 8 9 10 1 10 Z 2
6 ATOM C CA . GLY A 1 2 A 9 10 11 1 10 Z 2
7 ATOM C CA . THR A 1 4 . 10 11 12 1 12 Z 2
'''
    (directory / 'structure.cif').write_text(cif)
    protein = {'primaryAccession': 'SYNTHETIC', 'sequence': {'value': 'AGST'},
               'organism': {'scientificName': 'Synthetic software fixture'}}
    dump(directory / 'protein.json', protein)
    sources = [{'id': sid, 'path': filename, 'label': sid, 'retrieved_at': '2026-09-20',
                'sha256': digest((directory/filename).read_bytes())} for sid, filename in [
                    ('coordinates', 'structure.cif'), ('protein', 'protein.json')]]
    manifest = {'schema_version': '1.0', 'id': 'synthetic', 'title': 'Synthetic case',
        'protein': {'accession': 'SYNTHETIC', 'organism': protein['organism']['scientificName'],
                    'sequence': 'AGST', 'sequence_sha256': digest(b'AGST'), 'source_id': 'protein'},
        'sources': sources, 'structures': [{'id': 'structure', 'entry_id': 'TEST', 'source_id': 'coordinates',
            'sequence_mappings': [{'entity_id': '1', 'label_start': 1, 'reference_start': 1, 'length': 4}]}],
        'assertions': [], 'evidence': [], 'membranes': []}
    dump(directory/'case.json', manifest)
    return directory


def edit(directory, fn):
    p = directory/'case.json'; d = json.loads(p.read_text()); fn(d); dump(p,d)


def binding(directory, assembly='asu'):
    d=json.loads((directory/'case.json').read_text())
    return {'structure_id':'structure','coordinate_sha256':d['sources'][0]['sha256'],
            'sequence_sha256':d['protein']['sequence_sha256'],'model_id':'1','assembly_id':assembly}


def source(directory, sid, value):
    path=directory/(sid+'.json'); dump(path,value)
    row={'id':sid,'path':path.name,'label':sid,'retrieved_at':'2026-09-20','sha256':digest(path.read_bytes())}
    edit(directory,lambda d:d['sources'].append(row))
    return row


def test_generated_copies_preserve_exact_residue_and_operator_identity(case_dir):
    case=BiologicalCase(case_dir)
    view=case.view('structure','1','1')
    assert len(view['chains'])==2 and len(view['atoms'])==6
    a,b=view['chains']
    assert a['auth_asym_id']=='Z' and a['source_label_asym_id']=='A'
    assert a['id']!=b['id'] and a['residues'][1]['insertion_code']=='A'
    assert a['residues'][0]['auth_seq_id']==a['residues'][1]['auth_seq_id']=='10'
    assert a['residues'][0]['id']!=a['residues'][1]['id']
    assert a['residues'][2]['sequence_position']==3 and not a['residues'][2]['observed']
    assert a['residues'][2]['atom_indices']==[]
    assert a['residues'][1]['alternate_locations']==['A','B']
    assert a['residues'][1]['selected_altloc']=='A'
    for chain in view['chains']:
        for residue in chain['residues']:
            for i in residue['atom_indices']:
                assert view['atoms'][i]['properties']['residue_id']==residue['id']
                assert view['atoms'][i]['chain']==chain['id']
    assert [view['atoms'][3][k] for k in ('x','y','z')]==[19,-2,3]


def test_models_are_never_concatenated_or_assumed(case_dir):
    case=BiologicalCase(case_dir)
    assert case.summary()['structures'][0]['models']==['1','2']
    assert case.view('structure','2','asu')['atoms'][0]['x']==8
    for model,assembly in [('0','asu'),('1','unknown'),('1','')]:
        with pytest.raises(ValueError):case.view('structure',model,assembly)


@pytest.mark.parametrize('fault',['coordinates','sequence','source_sequence','mapping','overlap','entry','source_path','source_url'])
def test_corrupted_or_ambiguous_identity_is_rejected(case_dir,fault):
    if fault=='coordinates':
        with (case_dir/'structure.cif').open('a') as f:f.write('\n# changed\n')
    elif fault=='sequence':edit(case_dir,lambda d:d['protein'].update(sequence='AAAA'))
    elif fault=='source_sequence':edit(case_dir,lambda d:d['protein'].update(sequence='AAAA',sequence_sha256=digest(b'AAAA')))
    elif fault=='mapping':edit(case_dir,lambda d:d['structures'][0]['sequence_mappings'][0].update(reference_start=2))
    elif fault=='overlap':edit(case_dir,lambda d:d['structures'][0]['sequence_mappings'].append(d['structures'][0]['sequence_mappings'][0]))
    elif fault=='entry':edit(case_dir,lambda d:d['structures'][0].update(entry_id='OTHER'))
    elif fault=='source_path':edit(case_dir,lambda d:d['sources'][0].update(path='../structure.cif'))
    elif fault=='source_url':edit(case_dir,lambda d:d['sources'][0].update(url='javascript:alert(1)'))
    with pytest.raises(ValueError):BiologicalCase(case_dir)


def test_unmapped_positions_remain_unmapped(case_dir):
    edit(case_dir,lambda d:d['structures'][0].update(sequence_mappings=[]))
    view=BiologicalCase(case_dir).view('structure','1','asu')
    assert all(r['sequence_position'] is None and r['mapping_state']=='unmapped' for r in view['chains'][0]['residues'])


def test_conflicting_localization_assertions_are_retained(case_dir):
    assertions=[{'id':name,'property':'native_accessibility','value':value,'state':state,'basis':'curated',
        'source_id':'protein','scope':{'sequence_sha256':digest(b'AGST')},'interpretation_limit':'Synthetic fixture, no biological claim'}
        for name,value,state in [('positive',True,'supported'),('negative',False,'contradicted'),('unknown',None,'unknown')]]
    edit(case_dir,lambda d:d.update(assertions=assertions))
    case=BiologicalCase(case_dir)
    assert case.summary()['assertions']==assertions
    assert case.view('structure','1','asu')['membrane'] is None


def add_membrane(case_dir):
    b=binding(case_dir)
    source(case_dir,'placement',{'binding':b,'frame':'deposited_coordinates','normal':[0,0,1],
        'center':[0,0,0],'half_thickness':15,'sidedness':'unknown'})
    edit(case_dir,lambda d:d['membranes'].append({'binding':b,'source_id':'placement'}))


def test_placement_is_limited_to_its_coordinate_model_and_assembly(case_dir):
    add_membrane(case_dir)
    case=BiologicalCase(case_dir)
    assert case.view('structure','1','asu')['membrane']['sidedness']=='unknown'
    assert case.view('structure','1','1')['membrane'] is None
    assert case.view('structure','2','asu')['membrane'] is None


def test_soluble_case_refuses_bilayer_even_if_placement_is_present(case_dir):
    add_membrane(case_dir)
    edit(case_dir,lambda d:d['assertions'].append({'id':'soluble','property':'membrane_applicability',
        'value':False,'state':'not_applicable','basis':'curated','source_id':'protein',
        'scope':{'sequence_sha256':digest(b'AGST')},'interpretation_limit':'Synthetic control'}))
    with pytest.raises(ValueError,match='soluble'):BiologicalCase(case_dir)


def add_engine(case_dir):
    b=binding(case_dir)
    row=source(case_dir,'memorient',{'summary':{'host_antibody_accessible':True,'localization':'unknown'}})
    source(case_dir,'run',{'binding':b,'output_sha256':row['sha256'],'method':'memorient.orient_structure',
        'method_source_sha256':'a'*64})
    edit(case_dir,lambda d:d['evidence'].append({'id':'membrane','engine':'memorient','source_id':'memorient',
        'run_source_id':'run','binding':b}))


def test_legacy_exposure_is_raw_measurement_not_native_accessibility(case_dir):
    add_engine(case_dir)
    evidence=BiologicalCase(case_dir).view('structure','1','asu')['evidence'][0]
    assert evidence['summary']['host_antibody_accessible'] is True
    assert evidence['native_accessibility']=='unknown'


@pytest.mark.parametrize('field,value',[('coordinate_sha256','a'*64),('sequence_sha256','b'*64),('assembly_id','absent'),('model_id','2')])
def test_wrong_engine_scope_cannot_join(case_dir,field,value):
    add_engine(case_dir)
    edit(case_dir,lambda d:d['evidence'][0]['binding'].update({field:value}))
    with pytest.raises(ValueError):BiologicalCase(case_dir)


def test_store_refuses_traversal_and_symlink_cases(case_dir):
    store=BiologicalCaseStore(case_dir.parent.parent)
    assert store.list()[0]['state']=='ready'
    with pytest.raises(ValueError):store.load('../synthetic')
    (case_dir.parent/'alias').symlink_to(case_dir,target_is_directory=True)
    with pytest.raises(ValueError):store.load('alias')
    assert len(store.list())==1
