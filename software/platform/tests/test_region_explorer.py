"""Generic synthetic region, identity, review and relocated replay checks."""
import copy
import json
import shutil
import pytest
from test_biological_case import case_dir, edit, dump, digest
from yauvi_platform.structural_workbench.biological_case import BiologicalCase
from yauvi_platform.structural_workbench.regions import view_binding, select_region, measure_distance
from yauvi_platform.structural_workbench import evidence_tools as tools
from yauvi_platform.structural_workbench.region_workflow import build_case, main, auto_mappings
from yauvi_platform.structural_workbench import StructuralAnalysisStore


def annotate(directory):
    p=directory/'protein.json';reference=json.loads(p.read_text())
    reference['features']=[{'type':'Active site','location':{'start':{'value':2},'end':{'value':3}},'description':'Synthetic mapping test'},
        {'type':'Binding site','location':{'start':{'value':1,'modifier':'UNKNOWN'},'end':{'value':1}}}]
    dump(p,reference)
    edit(directory,lambda d:d['sources'][1].update(sha256=digest(p.read_bytes())))


def test_regions_keep_exact_missing_and_ambiguous_coverage(case_dir):
    annotate(case_dir);view=BiologicalCase(case_dir).view('structure','1','asu');region=view['regions'][0]
    assert region['observed_positions']==[2] and region['missing_positions']==[3]
    assert region['coverage']=={'observed':1,'expected':2} and region['atom_indices']==[1]
    assert view['regions'][1]['atom_indices']==[] and view['regions'][1]['mapping_state']=='unmapped'
    assert view['recorded_states']==[]


@pytest.mark.parametrize('field',['coordinate_sha256','reference_record_sha256','sequence_sha256','case_manifest_sha256','source_identities_sha256','model_id','assembly_id','frame','conformer_selection','mapping_sha256','method_sha256','atom_identity_sha256'])
def test_wrong_bound_identity_rejects_region_and_distance(case_dir,field):
    annotate(case_dir);view=BiologicalCase(case_dir).view('structure','1','asu');binding=view_binding(view);binding[field]='wrong'
    with pytest.raises(ValueError):select_region(view,view['regions'][0]['id'],binding)
    with pytest.raises(ValueError):measure_distance(view,0,1,binding)


def test_alternate_selection_cannot_reuse_atoms_or_legacy_evidence(case_dir):
    from test_biological_case import add_engine
    annotate(case_dir);add_engine(case_dir);case=BiologicalCase(case_dir)
    a=case.view('structure','1','asu','A');b=case.view('structure','1','asu','B');blank=case.view('structure','1','asu','blank')
    assert a['atoms'][1]['x']==2 and b['atoms'][1]['x']==5
    assert blank['regions'][0]['atom_indices']==[]
    assert a['evidence']==b['evidence']==blank['evidence']==[]
    with pytest.raises(ValueError):select_region(b,a['regions'][0]['id'],view_binding(a))
    with pytest.raises(ValueError):case.view('structure','1','asu','absent')


def add_component(directory):
    path=directory/'structure.cif';s=path.read_text()
    s += '8 HETATM O O1 . LIG B 2 . . 1 2 7 1 20 L 1 O1 LIG 20\n'
    s += '9 HETATM H H1 . LIG B 2 . . 1 2 3 1 20 L 1 H1 LIG 20\n'
    s += '\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\nLIG O1 O\nLIG O2 O\nLIG H1 H\n'
    path.write_text(s);edit(directory,lambda d:d['sources'][0].update(sha256=digest(path.read_bytes())))


def test_contacts_use_heavy_atoms_and_retain_missing_component_atoms(case_dir):
    add_component(case_dir);view=BiologicalCase(case_dir).view('structure','1','asu')
    r=next(r for r in view['regions'] if r['basis']=='measured_contact')
    assert r['missing_heavy_atom_names']==['O2'] and r['coverage']=={'observed':1,'expected':2}
    assert r['measurement']['minimum_distance_A']==pytest.approx(11**.5)
    assert r['measurement']['cutoff_A']==5 and r['measurement']['contact_residue_count']==3
    assert all(view['atoms'][p['ligand_atom_index']]['elem']!='H' for p in r['measurement']['closest_pairs'])
    pair=r['measurement']['closest_pairs'][0]
    assert measure_distance(view,pair['protein_atom_index'],pair['ligand_atom_index'],view_binding(view))['distance_A']==pair['distance_A']


def test_generic_workflow_mapping_and_report(case_dir,tmp_path):
    annotate(case_dir);out=tmp_path/'result'
    assert main(['--structure',str(case_dir/'structure.cif'),'--reference-record',str(case_dir/'protein.json'),'--out',str(out)])==0
    assert json.loads((out/'REGION_VIEW.json').read_text())['regions'][0]['observed_positions']==[2]
    assert (out/'REGIONS.tsv').exists()
    assert auto_mappings({'_entity_poly.entity_id':['1'],'_entity_poly.pdbx_seq_one_letter_code_can':['AAAA']},'AAAAAAAA')==[]
    with pytest.raises(ValueError):build_case(case_dir/'structure.cif',case_dir/'protein.json',tmp_path/'bad',model=2)


def test_registered_analysis_exports_regions_and_retains_incomplete_qc(case_dir,tmp_path):
    annotate(case_dir);store=StructuralAnalysisStore(tmp_path/'workbench')
    store.create('regions',analysis_type='region_explorer',question='Synthetic software mapping')
    store.add_file('regions',role='structure',path=case_dir/'structure.cif')
    store.add_file('regions',role='reference_record',path=case_dir/'protein.json')
    assert store.preflight('regions')['valid']
    result=store.run('regions');assert result['status']=='scientifically_incomplete'
    snapshot=store.snapshot('regions')
    assert any(d['path'].endswith('/REGION_VIEW.json') for d in snapshot['report']['documents'])
    extraction=store.artifact_path('regions',result['run_id'],'generated/region_reference/REFERENCE_EXTRACTION.json')
    assert json.loads(extraction.read_text())['source_record_sha256']==digest((case_dir/'protein.json').read_bytes())
    assert not any('No reference FASTA' in w for w in snapshot['report'].get('warnings',[]))
    import zipfile
    with zipfile.ZipFile(store.artifact_path('regions',result['run_id'],'RAW_EVIDENCE.zip')) as archive:
        assert 'outputs/region_explorer/REGION_VIEW.json' in archive.namelist()


def test_review_preserves_original_is_append_only_and_rejects_wrong_targets(case_dir,tmp_path):
    annotate(case_dir);case=BiologicalCase(case_dir);selection=tools.selectors(case)[0];view=tools.selected_view(case,selection)
    original=(case_dir/'case.json').read_bytes();target=view['regions'][0]
    record={'binding':view_binding(view),'adjudication':'unresolved','reviewer':'Synthetic reviewer','rationale':'Synthetic review requiring further evidence','target_id':target['id'],'target_sha256':tools.digest(target),'source_ids':['protein']}
    a=tools.review(case,selection,record,tmp_path/'reviews');b=tools.review(case,selection,record,tmp_path/'reviews')
    assert a==b and len(list((tmp_path/'reviews').iterdir()))==1 and (case_dir/'case.json').read_bytes()==original
    wrong=copy.deepcopy(record);wrong['binding']['assembly_id']='other'
    with pytest.raises(ValueError):tools.review(case,selection,wrong,tmp_path/'reviews')
    with pytest.raises(ValueError):tools.review(case,selection,record,case_dir/'reviews')
    wrong=copy.deepcopy(record);wrong['source_ids']=['unlocked']
    with pytest.raises(ValueError):tools.review(case,selection,wrong,tmp_path/'reviews')


def test_relocated_replay_recalculates_and_detects_substantive_differences(case_dir,tmp_path):
    annotate(case_dir);add_component(case_dir);bundle=tmp_path/'bundle';tools.bundle(BiologicalCase(case_dir),bundle)
    moved=tmp_path/'independent-location';shutil.copytree(bundle,moved);assert tools.replay(moved)['state']=='passed'
    expected=json.loads((moved/'EXPECTED_VIEWS.json').read_text());expected[0]['atoms'][0]['x']+=1;dump(moved/'EXPECTED_VIEWS.json',expected)
    manifest=json.loads((moved/'REPLAY_MANIFEST.json').read_text());manifest['files']['EXPECTED_VIEWS.json']=digest((moved/'EXPECTED_VIEWS.json').read_bytes());dump(moved/'REPLAY_MANIFEST.json',manifest)
    assert tools.replay(moved)['state']=='different'


@pytest.mark.parametrize('fault',['missing_file','checksum','runtime','method','symlink'])
def test_replay_fails_clearly_on_missing_or_wrong_identity(case_dir,tmp_path,fault):
    bundle=tmp_path/'bundle';tools.bundle(BiologicalCase(case_dir),bundle)
    if fault=='missing_file':(bundle/'case/protein.json').unlink()
    if fault=='checksum':(bundle/'case/protein.json').write_text('{}')
    if fault=='symlink':
        (bundle/'case/protein.json').unlink();(bundle/'case/protein.json').symlink_to(case_dir/'protein.json')
    if fault in ('runtime','method'):
        manifest=json.loads((bundle/'REPLAY_MANIFEST.json').read_text());manifest['runtime' if fault=='runtime' else 'method_locks']={};dump(bundle/'REPLAY_MANIFEST.json',manifest)
    with pytest.raises(ValueError):tools.replay(bundle)


def test_proper_assembly_motion_preserves_pair_distances(case_dir):
    view=BiologicalCase(case_dir).view('structure','1','1')
    binding=view_binding(view)
    a=measure_distance(view,0,1,binding)['distance_A']
    b=measure_distance(view,3,4,binding)['distance_A']
    assert a==pytest.approx(b,abs=1e-12)
    # Chain-copy atom identity is part of the binding, not just sequence number.
    assert view['atoms'][0]['chain']!=view['atoms'][3]['chain']


def test_review_individual_observation_pointer_is_exact(case_dir,tmp_path):
    annotate(case_dir);case=BiologicalCase(case_dir);selection=tools.selectors(case)[0];view=tools.selected_view(case,selection)
    target=view['regions'][0]
    record={'binding':view_binding(view),'adjudication':'unresolved','reviewer':'Synthetic reviewer','rationale':'Coverage needs evidence',
        'target_id':target['id'],'target_pointer':'/coverage','target_sha256':tools.digest(target['coverage']),'source_ids':['protein']}
    result=tools.review(case,selection,record,tmp_path/'reviews')
    assert json.loads(__import__('pathlib').Path(result['file']).read_text())['original_target']==target['coverage']
    record['target_pointer']='/absent'
    with pytest.raises(ValueError):tools.review(case,selection,record,tmp_path/'reviews')


def test_declared_chemical_difference_stays_excluded_from_exact_mapping():
    doc={'_entity_poly.entity_id':['1'],'_entity_poly.pdbx_seq_one_letter_code_can':['AGDT'],
         '_pdbx_poly_seq_scheme.entity_id':['1'],'_pdbx_poly_seq_scheme.seq_id':['3'],
         '_pdbx_poly_seq_scheme.pdb_strand_id':['A'],'_pdbx_poly_seq_scheme.mon_id':['MOD'],
         '_struct_ref_seq_dif.pdbx_seq_db_accession_code':['SYNTHETIC'],
         '_struct_ref_seq_dif.seq_num':['3'],'_struct_ref_seq_dif.pdbx_pdb_strand_id':['A'],
         '_struct_ref_seq_dif.mon_id':['MOD'],'_struct_ref_seq_dif.db_mon_id':['SER']}
    assert auto_mappings(doc,'AGST','OTHER')==[]
    assert auto_mappings(doc,'AGST','SYNTHETIC')==[
        {'entity_id':'1','label_start':1,'reference_start':1,'length':2},
        {'entity_id':'1','label_start':4,'reference_start':4,'length':1}]
    assert auto_mappings(doc,'ASST','SYNTHETIC')==[]  # Unreported mismatch is not ignored.


def test_atom_payload_mutation_cannot_receive_old_region_evidence(case_dir):
    annotate(case_dir);view=BiologicalCase(case_dir).view('structure','1','asu');binding=view_binding(view)
    view['atoms'][0]['x']+=1
    with pytest.raises(ValueError):select_region(view,view['regions'][0]['id'],binding)
    with pytest.raises(ValueError):measure_distance(view,0,1,binding)


def test_replay_reports_missing_dependency(case_dir,tmp_path,monkeypatch):
    bundle=tmp_path/'bundle';tools.bundle(BiologicalCase(case_dir),bundle)
    def missing():raise __import__('importlib.metadata',fromlist=['PackageNotFoundError']).PackageNotFoundError('synthetic-dependency')
    monkeypatch.setattr(tools,'runtime',missing)
    with pytest.raises(ValueError,match='dependency is absent'):tools.replay(bundle)


def test_bundled_region_example_runs_without_downloads(tmp_path):
    from yauvi_structural_workbench.example import create_region_example
    create_region_example(tmp_path,'example-regions')
    store=StructuralAnalysisStore(tmp_path)
    assert store.preflight('example-regions')['valid']
    result=store.run('example-regions')
    assert result['status']=='scientifically_incomplete'
    case=BiologicalCase(store.artifact_path('example-regions',result['run_id'],'outputs/region_explorer/CASE/case.json').parent)
    view=case.view('coordinates','1','asu')
    assert view['regions'][0]['missing_positions']==[3]
    assert view['regions'][-1]['missing_heavy_atom_names']==['O2']
    assert view['recorded_states']==[]


def test_static_state_requires_exact_frame_and_conformer_scope(case_dir):
    from test_biological_case import binding
    state={'id':'state','property':'conformation','value':{'label':'Source-declared synthetic static snapshot','construct':'synthetic construct'},
        'state':'supported','basis':'curated','source_id':'protein','scope':{**binding(case_dir),'frame':'deposited_coordinates','conformer_selection':'auto'},
        'interpretation_limit':'Synthetic source scope, no reaction history or probability'}
    edit(case_dir,lambda d:d['assertions'].append(state))
    case=BiologicalCase(case_dir)
    assert case.view('structure','1','asu')['recorded_states']==[state]
    assert case.view('structure','2','asu')['recorded_states']==[]
    assert case.view('structure','1','asu','blank')['recorded_states']==[]
    edit(case_dir,lambda d:d['assertions'][0]['scope'].update(frame='other_coordinates'))
    assert BiologicalCase(case_dir).view('structure','1','asu')['recorded_states']==[]
