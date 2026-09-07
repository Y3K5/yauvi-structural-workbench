"""Public-import scientific identity and no-overwrite boundaries, without network."""
import json
from types import SimpleNamespace
import pytest
from yauvi_platform.structural_workbench import StructuralAnalysisStore, AnalysisError
from yauvi_platform.structural_workbench.protein_import import ProteinImportStore, parse_uniprot_link


def fixture(tmp_path, *, no_model=False, sequence='AG', pae=None, duplicate=False):
    store=StructuralAnalysisStore(tmp_path)
    store.create('qc',analysis_type='structure_qc',question='Inspect a known protein')
    pdb=("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C  \n"
         "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00 85.00           C  \nEND\n")
    af={'uniprotAccession':'P02931','modelEntityId':'AF-P02931-F1','sequence':sequence,
        'sequenceStart':1,'sequenceEnd':len(sequence),'latestVersion':6,
        'pdbUrl':'https://alphafold.ebi.ac.uk/files/AF-P02931-F1-model_v6.pdb',
        'paeDocUrl':'https://alphafold.ebi.ac.uk/files/AF-P02931-F1-predicted_aligned_error_v6.json'}
    calls=[]
    def get(url,**kwargs):
        calls.append(url)
        if '/uniprotkb/' in url:
            doc={'primaryAccession':'P02931','sequence':{'value':'AG'},'proteinDescription':{'recommendedName':{'fullName':{'value':'Outer membrane protein F'}}},'organism':{'scientificName':'Escherichia coli'}}
        elif '/api/prediction/' in url:
            if no_model:return None,'not_found'
            doc=[af,af] if duplicate else [af]
        elif url.endswith('.pdb'):doc=pdb
        else:doc=[{'predicted_aligned_error':pae if pae is not None else [[0,1],[1,0]]}]
        content=doc.encode() if isinstance(doc,str) else json.dumps(doc).encode()
        return SimpleNamespace(content=content,url=url,headers={'X-UniProt-Release':'test-release'}),''
    importer=ProteinImportStore(store,get=get)
    return importer,store,calls,af


@pytest.mark.parametrize('value',['P02931','p02931','https://www.uniprot.org/uniprotkb/P02931/entry','https://rest.uniprot.org/uniprotkb/P02931.json'])
def test_links(value):assert parse_uniprot_link(value)=='P02931'


@pytest.mark.parametrize('value',['https://uniprot.org.evil.test/uniprotkb/P02931','http://127.0.0.1/P02931','https://x@www.uniprot.org/uniprotkb/P02931','https://www.uniprot.org:444/uniprotkb/P02931','P02931-2','https://www.uniprot.org/uniprotkb?query=P02931','../../private','https://www.uniprot.org/uniprotkb/P02931%2fentry'])
def test_reject_unsafe_and_ambiguous_links(value):
    with pytest.raises(ValueError):parse_uniprot_link(value)


def test_acquisition_separate_from_adoption_and_exact_provenance(tmp_path):
    importer,store,calls,af=fixture(tmp_path)
    original=store.load('qc');record=importer.lookup('P02931');assert not store.load('qc')['inputs']
    prepared=importer.prepare(record['import_id']);assert not store.load('qc')['inputs']
    result=importer.adopt(record['import_id'],'qc',original['revision_sha256'])
    inputs={r['role']:r for r in result['analysis']['inputs']}
    assert set(inputs)=={'structure','reference_fasta','provenance','pae'}
    provenance=json.loads(store.input_path('qc',inputs['provenance']['sha256']).read_text())
    assert provenance['class']=='predicted' and provenance['coordinate_sha256']==inputs['structure']['sha256']
    assert provenance['release']=='6' and result['analysis']['parameters']['chain']=='A'
    assert result['analysis']['protein_import']['metadata'][0]['artifact']['release']=='test-release'
    assert 'wwPDB or local validation report' in result['missing']
    assert len(calls)==4
    importer.prepare(record['import_id']);assert len(calls)==4


def test_missing_model_retains_sequence_and_explicit_gap(tmp_path):
    importer,store,calls,_=fixture(tmp_path,no_model=True)
    record=importer.prepare(importer.lookup('P02931')['import_id'])
    assert [f['role'] for f in record['files']]==['reference_fasta']
    assert record['warnings'] and 'provenance' not in record


@pytest.mark.parametrize('kwargs',[{'sequence':'AA'},{'duplicate':True}])
def test_mismatched_or_ambiguous_model_not_selected(tmp_path,kwargs):
    importer,store,calls,_=fixture(tmp_path,**kwargs)
    record=importer.prepare(importer.lookup('P02931')['import_id'])
    assert [f['role'] for f in record['files']]==['reference_fasta']
    assert not any(url.endswith('.pdb') for url in calls)


def test_bad_pae_not_adopted(tmp_path):
    importer,store,calls,_=fixture(tmp_path,pae=[[0]])
    record=importer.prepare(importer.lookup('P02931')['import_id'])
    assert 'pae' not in [f['role'] for f in record['files']]
    assert record['provenance'] and 'PAE' in record['warnings'][0]


def test_external_model_url_not_requested(tmp_path):
    importer,store,calls,af=fixture(tmp_path);af['pdbUrl']='https://evil.test/model.pdb'
    record=importer.prepare(importer.lookup('P02931')['import_id'])
    assert not any('evil.test' in url for url in calls)
    assert record['warnings']


def test_stale_revision_and_conflicting_inputs_unchanged(tmp_path):
    importer,store,calls,_=fixture(tmp_path)
    record=importer.prepare(importer.lookup('P02931')['import_id']);initial=store.load('qc')
    store.update_parameters('qc',{'chain':'B'})
    with pytest.raises(AnalysisError,match='changed'):importer.adopt(record['import_id'],'qc',initial['revision_sha256'])
    assert not store.load('qc')['inputs']
    p=tmp_path/'different.fasta';p.write_text('>other\nAA\n');store.add_file('qc',role='reference_fasta',path=p)
    before=store.load('qc')
    with pytest.raises(AnalysisError,match='differs'):importer.adopt(record['import_id'],'qc',before['revision_sha256'])
    assert store.load('qc')==before


def test_corrupt_cache_refused_before_any_adoption(tmp_path):
    importer,store,calls,_=fixture(tmp_path);record=importer.prepare(importer.lookup('P02931')['import_id'])
    acquisition=record['files'][-1]['acquisition'];path=importer.sources.artifact_path(acquisition['acquisition_id']);path.write_text('corrupt')
    with pytest.raises(Exception,match='integrity|checksum|verify|manifest'):
        importer.adopt(record['import_id'],'qc',store.load('qc')['revision_sha256'])
    assert not store.load('qc')['inputs']


def test_sf_roles_and_metadata_filled_without_inventing_science(tmp_path):
    importer,store,_,_=fixture(tmp_path);store.create('sf',analysis_type='sf_csa',question='Inspect')
    record=importer.prepare(importer.lookup('P02931')['import_id']);result=importer.adopt(record['import_id'],'sf',store.load('sf')['revision_sha256'])
    assert set(r['role'] for r in result['analysis']['inputs'])=={'query_structure','query_fasta','provenance'}
    params=result['analysis']['parameters'];assert params['accession']=='P02931' and params['organism']=='Escherichia coli'
    assert not params.get('mechanism_group') and not params.get('database_pack')


def test_requires_per_request_public_opt_in():
    from yauvi_structural_workbench.server import Handler
    handler=object.__new__(Handler);handler.server=SimpleNamespace(store=None)
    handler._body=lambda:{'link':'P02931'}
    with pytest.raises(PermissionError,match='Retrieve protein'):handler._mutate(['api','proteins'])


def test_membrane_snapshot_exposes_oriented_coordinate_artifact(tmp_path):
    store=StructuralAnalysisStore(tmp_path)
    case=store.create('mem',analysis_type='membrane_orientation',question='Orient a model')
    run_id='run-'+'a'*16
    run_dir=store.cases_root/'mem'/'runs'/run_id/'outputs'/'memorient'
    run_dir.mkdir(parents=True)
    (run_dir/'ORIENTED_STRUCTURE.pdb').write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C  \\nEND\\n')
    (run_dir/'MEMBRANE_LAYER.json').write_text('{"membrane_slab":{"half_thickness":12.0}}')
    manifest=store.load('mem');manifest['latest_run_id']=run_id
    store._commit(store._case_dir('mem'),manifest)
    snapshot=store.snapshot('mem')
    assert snapshot['oriented_structure']['run_id']==run_id
    assert snapshot['oriented_structure']['file_name']=='ORIENTED_STRUCTURE.pdb'
    assert snapshot['membrane_layer']['run_id']==run_id
    assert snapshot['membrane_layer']['file_name']=='MEMBRANE_LAYER.json'
    assert store.artifact_path('mem',run_id,'outputs/memorient/ORIENTED_STRUCTURE.pdb').is_file()
    assert store.artifact_path('mem',run_id,'outputs/memorient/MEMBRANE_LAYER.json').is_file()
