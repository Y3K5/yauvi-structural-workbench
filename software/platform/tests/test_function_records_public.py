"""Invented data only: public record, custody and missingness software controls."""
import hashlib,json,socket,subprocess,zipfile
from copy import deepcopy
from io import BytesIO
import pytest
from yauvi_platform.structural_workbench.function_records import record,browse,bundle,derive
from yauvi_platform.structural_workbench.function_investigations import FunctionInvestigations
from yauvi_platform.structural_workbench.store import AnalysisError,_canonical

@pytest.fixture
def functions(tmp_path,monkeypatch):
 monkeypatch.delenv('YAUVI_PROTEIN_FUNCTION_EVIDENCE',raising=False)
 data={'schema':'functional-protein-catalog-v1','dataset_id':'a'*64,'source_ledger_sha256':'b'*64,'groups':[{'id':'attachment','label':'Attachment'},{'id':'envelope','label':'Envelope'}],'species':[{'id':'a','label':'Invented microbe A'},{'id':'b','label':'Invented microbe B'}],'records':[],'evidence':{},'profiles':{},'source_relationships':[],'source_references':[]}
 for i in range(4):
  key=f'a|P{i}';data['records'].append({'key':key,'id':f'P{i}','name':'Invented protein','species':'a','strain':'invented-A','primary_group':'attachment' if i<2 else 'envelope','length':20,'sequence_sha256':str(i)*64,'review_category':'held','holds':['source unresolved'],'epitope_state':'unevaluated'})
  data['evidence'][key]={'structure':{'status':'not_joined','available':False},'location':{'status':'unknown','zone':'unplaced'}}
 path=tmp_path/'catalogue.json';path.write_text(json.dumps(data));cfg=tmp_path/'config.json';cfg.write_text(json.dumps({'schema':'function-research-config-v1','catalogue_path':str(path),'catalogue_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}));monkeypatch.setenv('YAUVI_FUNCTION_RESEARCH_CONFIG',str(cfg))
 return FunctionInvestigations(tmp_path/'workspace'),data,path,cfg

def proposal(key='a|P0',**extra):
 r=record(key)
 return {'key':key,'record_sha256':r['record_sha256'],'question':'Invented functional question','action_id':'domains','hypotheses':[{'statement':'Invented proposed role','alternative':'An alternative role','distinguishing_evidence':'A source-qualified distinguishing observation','assertion_ids':[r['assertions'][0]['id']]}],'notes':'Invented software check','origin_study':None,**extra}


def supplement(functions,monkeypatch,kind='refinement',changes=None):
 _,data,path,_=functions;r=data['records'][0]
 row={**deepcopy(r),'detail_group':'named_family','detail_level':'family_hint','assignments':[{'text':'Invented family','level':'family_hint'}],'inventory':[{'dimension':'Function','state':'not_joined','detail':'Invented missing activity'}],'limit':'Family is not activity'}
 if kind=='duf':row.update(assessment='cached_function_inference',proposed_role='Invented reference annotation',channels={'function_comment':True},remaining_unknowns=['No native function established'])
 if kind=='bam':row.update(group_label='Invented assembly',mechanism='Proposed scaffold role',hold='Reference strain only',assignment_state='inferred_reference_supported',measured_w83_function=False)
 if changes:row.update(changes)
 schemas={'refinement':'protein-function-refinement-v1','duf':'local-duf-function-review-v1','bam':'bam-domain-function-review-v1'}
 doc={'schema':schemas[kind],'base_dataset_id':data['dataset_id'],'records':{r['key']:row}}
 source=path.parent/'extra.json';source.write_text(json.dumps(doc));cfg={'schema':'protein-function-evidence-config-v1','dataset_id':data['dataset_id'],'catalogue_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'sources':{kind:{'path':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()}}};cp=path.parent/'extra-config.json';cp.write_text(json.dumps(cfg));monkeypatch.setenv('YAUVI_PROTEIN_FUNCTION_EVIDENCE',str(cp));return source,cp,doc


def test_every_member_has_record_with_explicit_gaps_and_stable_identity(functions):
 service,data,path,_=functions;before=path.read_bytes()
 for row in data['records']:
  r=record(row['key']);identity=r.pop('record_sha256');assert identity==hashlib.sha256(_canonical(r)).hexdigest();assert r['original_record']==row and r['function']['state']=='annotation_only';assert r['partners']==[];assert any(u['dimension']=='partners' and u['state']=='not_joined' for u in r['unknowns']);assert r['binding']['supplement_state']['state']=='not_configured';assert r['identity']['sequence_sha256']==row['sequence_sha256']
 assert path.read_bytes()==before and service.studies.list()==[]


def test_exact_optional_review_joins_without_new_function_assignment(functions,monkeypatch):
 supplement(functions,monkeypatch);r=record('a|P0');assert r['binding']['supplement_state']['state']=='current';assert r['retained_reviews']['refinement']['detail_level']=='family_hint';assert r['original_record']['review_category']=='held';assert r['function']['state']=='annotation_only';a=next(a for a in r['assertions'] if a['predicate']=='annotation_refinement');assert a['source_kind']=='annotation_organization' and a['source']['json_pointer']=='/records/a|P0';assert browse(gap='specificity')['matching']==1


@pytest.mark.parametrize('change',[{'sequence_sha256':'f'*64},{'key':'outside'},{'length':999},{'id':'different-accession'},{'species':'b'}])
def test_misbound_review_is_unavailable_not_negative(functions,monkeypatch,change):
 service,_,_,_=functions;supplement(functions,monkeypatch,changes=change);r=record('a|P0');assert r['binding']['supplement_state']['state']=='unavailable' and not r['retained_reviews'];assert any(x['dimension']=='review_coverage' for x in r['unknowns'])
 with pytest.raises(AnalysisError):service.save(proposal())
 assert service.studies.list()==[]


def test_changed_optional_bytes_hold_save_but_original_snapshot_reads(functions,monkeypatch):
 service,_,_,_=functions;source,cfg,_=supplement(functions,monkeypatch);s=service.save(proposal());old=service.load(s['study_id']);source.write_text('{}');after=service.load(s['study_id']);assert after['packet']==old['packet'] and after['source_state']=='changed';assert record('a|P0')['binding']['supplement_state']['state']=='unavailable'
 with pytest.raises(AnalysisError):service.save(proposal(study_id=s['study_id'],expected_revision_sha256=s['revision_sha256']))
 assert service.studies.load(s['study_id'])==s


def test_source_scope_and_counterpart_material_are_not_functional_transfer(functions,monkeypatch):
 _,data,path,cfg=functions;data=deepcopy(data);data['profiles']['a|P0']={'role':'Reported binding role','join':'Nonidentical source-strain counterpart','caveat':'No native target assay','study_context':{'connections':[{'study_strains':'Source strain X','unknowns':['Target counterpart remains unresolved']} ]},'sources':[{'id':'x'}],'connections':[{'partner_name':'Partner X','relation':'binding','level':'direct assay','claim':'Source material bound a fragment','limit':'Not the target strain','source_ids':['x']} ]};path.write_text(json.dumps(data));c=json.loads(cfg.read_text());c['catalogue_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();cfg.write_text(json.dumps(c));r=record('a|P0');assert r['partners'][0]['limit']=='Not the target strain';a=next(a for a in r['assertions'] if a['predicate']=='reported_relationship');assert a['scope']==data['profiles']['a|P0']['study_context'];assert r['function']['state']=='annotation_only' and r['reference_transfer']=='Nonidentical source-strain counterpart'


@pytest.mark.parametrize('kwargs',[{'gap':'outside'},{'species':'outside'},{'group':'outside'},{'page':'-1'},{'page':'1'},{'page':'9999999'},{'q':'x'*201}])
def test_invalid_filters_refused(functions,kwargs):
 with pytest.raises(AnalysisError):browse(**kwargs)


def test_browse_keeps_original_groups_and_missing_coverage(functions):
 assert browse()['total']==4 and browse(gap='function')['matching']==0 and browse(gap='partners')['matching']==4 and browse(gap='structure')['matching']==4;assert browse(group='attachment')['matching']==2 and browse(q='P0')['matching']==1 and browse(q='not-present')['matching']==0
 with pytest.raises(AnalysisError):record('outside')


def test_plan_revision_and_exact_export_keep_sources_and_prior_revisions(functions):
 service,_,path,_=functions;before=path.read_bytes();s=service.save(proposal());old=deepcopy(s);request=proposal(study_id=s['study_id'],expected_revision_sha256=s['revision_sha256'],question='Revised invented question');s=service.save(request);assert s['revision']==2 and service.load(s['study_id'],old['revision_sha256'])['study']==old;z=zipfile.ZipFile(BytesIO(service.studies.bundle(s['study_id'])));assert b'Protein function record and unknown-function investigation' in z.read('STUDY_REPORT.html');p=json.loads(z.read('objects/'+s['protein_function']['packet_sha256']));assert p['plan']['hypotheses'][0]['status']=='user_proposed_not_measured';assert not s['evidence'] and not s['scientifically_approved'];assert path.read_bytes()==before
 for n,h in json.loads(z.read('CHECKSUMS.json'))['files'].items():assert hashlib.sha256(z.read(n)).hexdigest()==h
 with pytest.raises(AnalysisError):service.save(request)
 assert service.studies.load(s['study_id'])==s


@pytest.mark.parametrize('change',[{'record_sha256':'f'*64},{'action_id':'run_predictor'},{'question':''},{'hypotheses':[]},{'hypotheses':[{'statement':'A','alternative':'B','distinguishing_evidence':'C','assertion_ids':['invented-unsupported']}]},{'hypotheses':[{'statement':'A','alternative':'B','distinguishing_evidence':'C','assertion_ids':[],'status':'measured'}]},{'notes':'x'*4001},{'unexpected':'field'},{'expected_revision_sha256':'f'*64}])
def test_invalid_plans_refuse_before_any_native_write(functions,change):
 service,_,_,_=functions
 with pytest.raises(AnalysisError):service.save(proposal(**change))
 assert service.studies.list()==[]


def test_catalogue_change_keeps_stored_investigation_readable(functions):
 service,_,path,_=functions;s=service.save(proposal());before=service.load(s['study_id']);path.write_text('{}');after=service.load(s['study_id']);assert after['packet']==before['packet'] and after['source_state']=='unavailable'


def test_valid_new_evidence_does_not_rewrite_old_investigation(functions,monkeypatch):
 service,_,_,_=functions;source,cfg,doc=supplement(functions,monkeypatch);s=service.save(proposal());doc['records']['a|P0']['limit']='Revised source scope';source.write_text(json.dumps(doc));c=json.loads(cfg.read_text());c['sources']['refinement']['sha256']=hashlib.sha256(source.read_bytes()).hexdigest();cfg.write_text(json.dumps(c));assert service.load(s['study_id'])['source_state']=='changed'
 with pytest.raises(AnalysisError):service.save(proposal(study_id=s['study_id'],expected_revision_sha256=s['revision_sha256']))
 assert service.studies.load(s['study_id'])==s


def test_report_escapes_user_text_and_performs_no_scientific_dispatch(functions,monkeypatch):
 service,_,_,_=functions
 def forbid(*a,**k):raise AssertionError('No scientific process or network allowed')
 monkeypatch.setattr(subprocess,'run',forbid);monkeypatch.setattr(socket,'create_connection',forbid);s=service.save(proposal(question='<script>invented</script>'));report=zipfile.ZipFile(BytesIO(service.studies.bundle(s['study_id']))).read('STUDY_REPORT.html');assert b'&lt;script&gt;invented&lt;/script&gt;' in report and b'<script>' not in report


def test_complete_source_graph_is_joined_and_duplicate_profile_edge_not_counted_twice(functions,monkeypatch):
 _,data,path,cfg=functions;data=deepcopy(data);edge={'id':'source-edge','a':'source-a','b':'host-matrix','relation':'binding','level':'direct assay','claim':'A source material bound matrix','limit':'Reference material only','source_ids':['paper']};data['source_relationships']=[edge];data['profiles']['a|P0']={'atlas_id':'source-a','join':'A nonidentical counterpart','connections':[{**edge,'partner_name':'Matrix protein'}]};path.write_text(json.dumps(data));c=json.loads(cfg.read_text());c['catalogue_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();cfg.write_text(json.dumps(c));r=record('a|P0');assert len(r['partners'])==1 and r['partners'][0]['name']=='Matrix protein';a=next(a for a in r['assertions'] if a['predicate']=='reported_relationship');assert a['source']['json_pointer']=='/source_relationships/0';assert browse(gap='partners')['matching']==3;assert not record('a|P1')['partners']


@pytest.mark.parametrize('kind',['duf','bam'])
def test_original_review_methods_are_retained_without_running_them(functions,monkeypatch,kind):
 source,cfg,doc=supplement(functions,monkeypatch,kind);doc.update(policy={'threshold':'Original descriptive screen only'},tools={'engine':'Historical engine version'},execution={'status':'recorded predecessor execution'},methods=['Original correspondence method'],evidence=[{'scope':'Reference strain, not selected target','limit':'No native target assay'}]);source.write_text(json.dumps(doc));c=json.loads(cfg.read_text());c['sources'][kind]['sha256']=hashlib.sha256(source.read_bytes()).hexdigest();cfg.write_text(json.dumps(c));r=record('a|P0');assert r['function']['state']=='reference_supported_hypothesis' and r['original_record']['review_category']=='held'
 if kind=='duf':assert r['review_method_context'][kind]['policy']==doc['policy'] and r['review_method_context'][kind]['execution']==doc['execution']
 else:assert r['review_method_context'][kind]['evidence']==doc['evidence']

@pytest.mark.parametrize('registry,change',[('groups','duplicate'),('species','duplicate'),('groups','unsafe'),('species','unsafe')])
def test_duplicate_or_unsafe_registry_ids_refuse_before_record_generation(functions,registry,change):
 _,data,path,cfg=functions;data=deepcopy(data)
 if change=='duplicate':data[registry].append(deepcopy(data[registry][0]))
 else:data[registry][0]['id']='unsafe" attribute'
 path.write_text(json.dumps(data));config=json.loads(cfg.read_text());config['catalogue_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();cfg.write_text(json.dumps(config))
 with pytest.raises(AnalysisError):record('a|P0')
