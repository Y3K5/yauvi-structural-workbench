"""Failability checks for context-bound public references, entirely offline."""
import json
from types import SimpleNamespace

import pytest

from yauvi_platform.structural_workbench import StructuralAnalysisStore
from yauvi_platform.structural_workbench.protein_import import ProteinImportStore, parse_uniprot_link
from yauvi_platform.structural_workbench.reference_context import proteome_identifier, validate_pdb

PID = 'UP000000625'


def setup(tmp_path):
    store = StructuralAnalysisStore(tmp_path)
    store.create('qc', analysis_type='structure_qc', question='Context checks')
    entry = {'primaryAccession':'P02931', 'sequence':{'value':'AG'},
        'organism':{'scientificName':'Test organism', 'taxonId':83333},
        'entryAudit':{'entryVersion':7,'sequenceVersion':2},
        'uniProtKBCrossReferences':[
            {'database':'Proteomes','id':PID,'properties':[{'key':'Component','value':'Chromosome'}]},
            {'database':'PDB','id':'1ABC','properties':[{'key':'Chains','value':'A=1-2'},
                {'key':'Method','value':'X-ray'},{'key':'Resolution','value':'2.1 A'}]}]}
    meta = {'id':PID,'taxonomy':{'taxonId':83333,'scientificName':'Test organism'},
        'strain':'test strain','proteomeType':'Other proteome','proteinCount':2,
        'genomeAssembly':{'assemblyId':'GCA_000005845.2'}}
    cif = b'data_1ABC\n_entry.id 1ABC\n_exptl.method "X-RAY DIFFRACTION"\n_atom_site.Cartn_x 1.0\n'
    payload = {'fasta':b'>sp|P02931|TEST OX=83333\nAG\n>tr|P12345|OTHER OX=83333\nAA\n',
               'release':'2026_03','cif':cif}
    calls=[]
    def get(url, **kwargs):
        calls.append(url)
        if '/api/prediction/' in url:
            return None, 'not_found'
        if '/uniprotkb/stream?' in url:
            body = payload['fasta']
        elif '/proteomes/' in url:
            body = json.dumps(meta).encode()
        elif url.endswith('.cif'):
            body = payload['cif']
        elif url.endswith('.pdb'):
            return None, 'not_found'
        else:
            body = json.dumps(entry).encode()
        release = payload['release'] if '/stream?' in url else '2026_03'
        return SimpleNamespace(content=body,url=url,headers={'X-UniProt-Release':release}), ''
    importer = ProteinImportStore(store, get=get)
    return importer, store, entry, meta, payload, calls


@pytest.mark.parametrize('link',[PID,PID.lower(),'https://www.uniprot.org/proteomes/'+PID,
    'https://rest.uniprot.org/proteomes/'+PID+'.json'])
def test_proteome_link_recognition(link):
    assert proteome_identifier(link) == PID


@pytest.mark.parametrize('link',['https://www.uniprot.org/proteomes/'+PID+'?query=all',
    'https://www.uniprot.org.evil.test/proteomes/'+PID,'https://user@www.uniprot.org/proteomes/'+PID,
    'http://www.uniprot.org/proteomes/'+PID,'UP00000564','UP000000625/../../private',
    'https://www.uniprot.org/proteomes/'+PID+'#fragment'])
def test_invalid_or_ambiguous_proteome_links_never_fetch(tmp_path,link):
    importer,*_,calls=setup(tmp_path)
    with pytest.raises(ValueError):importer.lookup(link)
    assert calls == []


def test_single_proteome_lookup_is_metadata_only_until_explicit_selection(tmp_path):
    importer,store,_,_,_,calls=setup(tmp_path)
    record=importer.prepare(importer.lookup(PID)['import_id'])
    assert len(calls)==1 and '/proteomes/' in calls[0]
    assert record['context']['proteome_type']=='Other proteome'
    assert record['context']['strain']=='test strain'
    assert not record['files'] and not store.load('qc')['inputs']
    record=importer.retrieve(record['import_id'],'uniprot.proteome:'+PID)
    row=record['retrievals'][0]
    assert row['checks']['record_count']==2
    assert 'includeIsoform=false' in calls[-1] and 'proteome%3A'+PID in calls[-1]
    assert row['acquisition']['artifact']['retrieved_at']
    assert row['metadata']['artifact']['sha256']==row['parent_metadata_sha256']
    with pytest.raises(ValueError,match='whole proteome'):
        importer.adopt(record['import_id'],'qc',store.load('qc')['revision_sha256'])
    assert not store.load('qc')['inputs']


def test_protein_retains_versions_relationships_and_exact_membership(tmp_path):
    importer,store,_,_,_,calls=setup(tmp_path)
    record=importer.prepare(importer.lookup('P02931')['import_id'])
    assert record['context']['entry_audit']['sequenceVersion']==2
    assert record['resources'][1]['context']['chain_coverage']=='A=1-2'
    assert not any('/stream?' in u or 'files.rcsb' in u for u in calls)
    record=importer.retrieve(record['import_id'],'uniprot.proteome:'+PID)
    assert record['retrievals'][0]['checks']['protein_membership']=='exact accession and sequence matched'
    count=len(calls)
    importer.retrieve(record['import_id'],'uniprot.proteome:'+PID)
    assert len(calls)==count
    exported=importer.export_record(record['import_id'])
    assert exported['metadata'][0]['artifact']['release']=='2026_03'
    assert not store.load('qc')['inputs']


@pytest.mark.parametrize('fault',['taxon','metadata_id','count','truncated','absent','changed_sequence',
    'wrong_fasta_taxon','duplicate','isoform','html','release','unknown_release'])
def test_proteome_mismatches_fail_without_adoption(tmp_path,fault):
    importer,store,entry,meta,payload,calls=setup(tmp_path)
    record=importer.prepare(importer.lookup('P02931')['import_id'])
    if fault=='taxon':meta['taxonomy']['taxonId']=9606
    elif fault=='metadata_id':meta['id']='UP000005640'
    elif fault=='count':meta['proteinCount']=3
    elif fault=='truncated':payload['fasta']=payload['fasta'].split(b'>tr')[0]
    elif fault=='absent':payload['fasta']=payload['fasta'].replace(b'P02931',b'P02932')
    elif fault=='changed_sequence':payload['fasta']=payload['fasta'].replace(b'\nAG\n',b'\nAA\n')
    elif fault=='wrong_fasta_taxon':payload['fasta']=payload['fasta'].replace(b'OX=83333',b'OX=9606')
    elif fault=='duplicate':payload['fasta']=payload['fasta'].replace(b'P12345',b'P02931')
    elif fault=='isoform':payload['fasta']=payload['fasta'].replace(b'P12345',b'P12345-2')
    elif fault=='html':payload['fasta']=b'<html>Service unavailable</html>'
    elif fault=='release':payload['release']='2026_04'
    else:payload['release']='unknown'
    with pytest.raises(ValueError,match='retrieval stopped'):
        importer.retrieve(record['import_id'],'uniprot.proteome:'+PID)
    saved=importer.load(record['import_id'])
    assert saved['retrievals'][-1]['state']=='failed'
    assert 'acquisition' not in saved['retrievals'][-1]
    assert not store.load('qc')['inputs']


def test_experimental_reference_requires_selection_and_preserves_mapping_gap(tmp_path):
    importer,store,_,_,_,calls=setup(tmp_path)
    record=importer.prepare(importer.lookup('P02931')['import_id'])
    before=len(calls)
    with pytest.raises(ValueError,match='listed'):
        importer.retrieve(record['import_id'],'pdb.coordinates:2XYZ')
    assert len(calls)==before
    record=importer.retrieve(record['import_id'],'pdb.coordinates:1ABC')
    row=record['retrievals'][-1]
    assert row['checks']['entry_id']=='1ABC'
    assert row['checks']['chain_mapping']=='requires review'
    assert row['context']['chain_coverage']=='A=1-2'
    assert not any(f['role']=='structure' for f in record['files'])
    assert not store.load('qc')['inputs']


@pytest.mark.parametrize('body',[b'<html>Error</html>',b'data_other\n_entry.id 2XYZ\n_atom_site.Cartn_x 1\n'])
def test_wrong_pdb_response_is_not_accepted(tmp_path,body):
    importer,store,_,_,payload,_=setup(tmp_path)
    record=importer.prepare(importer.lookup('P02931')['import_id'])
    payload['cif']=body
    with pytest.raises(ValueError):importer.retrieve(record['import_id'],'pdb.coordinates:1ABC')
    assert importer.load(record['import_id'])['retrievals'][-1]['state']=='failed'


def test_unavailable_legacy_pdb_is_recorded(tmp_path):
    importer,*_=setup(tmp_path)
    record=importer.prepare(importer.lookup('P02931')['import_id'])
    with pytest.raises(ValueError,match='not_found'):
        importer.retrieve(record['import_id'],'pdb.legacy:1ABC')
    assert importer.export_record(record['import_id'])['retrievals'][-1]['reason']=='not_found'


def test_legacy_pdb_identity_checked():
    body=('HEADER'.ljust(62)+'1ABC\nATOM      1  CA\n').encode()
    assert validate_pdb(body,'1ABC',legacy=True)['format']=='legacy PDB'
    with pytest.raises(ValueError):validate_pdb(body,'2XYZ',legacy=True)


def test_tampered_metadata_and_exported_files_refused(tmp_path):
    importer,*_=setup(tmp_path)
    record=importer.prepare(importer.lookup(PID)['import_id'])
    metadata=record['metadata'][0]
    importer.sources.artifact_path(metadata['acquisition_id']).write_bytes(b'{}')
    with pytest.raises(Exception,match='checksum|integrity|manifest'):
        importer.retrieve(record['import_id'],'uniprot.proteome:'+PID)
    with pytest.raises(Exception,match='checksum|integrity|manifest'):
        importer.export_record(record['import_id'])


def test_related_fetch_requires_explicit_opt_in():
    from yauvi_structural_workbench.server import Handler
    handler=object.__new__(Handler);handler.server=SimpleNamespace(store=None)
    handler._body=lambda:{'resource_id':'uniprot.proteome:'+PID}
    with pytest.raises(PermissionError):handler._mutate(['api','proteins','protein_'+'a'*24,'retrieve'])
