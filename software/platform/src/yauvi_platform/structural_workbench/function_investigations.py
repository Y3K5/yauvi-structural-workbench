"""Proposed investigations over injected native custody; no scientific execution."""
from datetime import datetime, timezone
import hashlib,json,secrets
from .store import AnalysisError,_canonical,_serialized
from .function_record_store import _text
from .function_records import record,ACTIONS,INVESTIGATION

class FunctionInvestigations:
    def __init__(self,workspace,*,studies=None):
        if studies is None:
            from .function_record_store import FunctionRecordStore
            studies=FunctionRecordStore(workspace)
        self.studies=studies
    def _transaction(self):return self.studies._transaction()
    def list(self,key):
        return [{'study_id':s['study_id'],'title':s['title'],'question':s['question'],'revision':s['revision'],'revision_sha256':s['revision_sha256']} for s in self.studies.list() if s.get('protein_function',{}).get('schema')==INVESTIGATION and s['protein_function']['key']==key]
    def load(self,study_id,revision=None):
        s=self.studies.snapshot(study_id,revision_sha256=revision)['study'];f=s.get('protein_function',{})
        if f.get('schema')!=INVESTIGATION:raise AnalysisError('Not a protein-function investigation')
        _,raw=self.studies.artifact(study_id,f['packet_sha256']);packet=json.loads(raw);state='unavailable'
        try:state='current' if record(f['key'])['record_sha256']==packet['function_record']['record_sha256'] else 'changed'
        except (AnalysisError,OSError,ValueError,KeyError):pass
        return {'study':s,'packet':packet,'source_state':state}
    @_serialized
    def save(self,request):
        if not isinstance(request,dict) or set(request)-{'key','record_sha256','study_id','expected_revision_sha256','question','action_id','hypotheses','notes','origin_study'}:raise AnalysisError('Use a bounded function-investigation request')
        current=record(request.get('key'))
        if request.get('record_sha256')!=current['record_sha256']:raise AnalysisError('Protein evidence changed; reopen before saving')
        if current['binding']['supplement_state']['state']=='unavailable':raise AnalysisError('Review unavailable optional evidence before saving a new investigation')
        question=_text(request.get('question'),'Investigation question');action=request.get('action_id')
        if not isinstance(action,str) or action not in {x['id'] for x in ACTIONS}:raise AnalysisError('Choose a recorded investigation action')
        notes=request.get('notes','')
        if not isinstance(notes,str) or len(notes)>4000:raise AnalysisError('Investigation notes exceed the local limit')
        hypotheses=request.get('hypotheses');qualified=[];ids={a['id'] for a in current['assertions']}
        if not isinstance(hypotheses,list) or not 1<=len(hypotheses)<=5:raise AnalysisError('Record one to five competing-hypothesis plans')
        for h in hypotheses:
            if not isinstance(h,dict) or set(h)!={'statement','alternative','distinguishing_evidence','assertion_ids'}:raise AnalysisError('Declare a hypothesis, alternative, distinguishing evidence and exact support links')
            support=h['assertion_ids']
            if not isinstance(support,list) or len(support)>100 or any(not isinstance(x,str) or x not in ids for x in support) or len(set(support))!=len(support):raise AnalysisError('Hypothesis support must name retained assertion identities')
            qualified.append({**{k:_text(h[k],k) for k in ['statement','alternative','distinguishing_evidence']},'assertion_ids':support,'status':'user_proposed_not_measured'})
        origin=request.get('origin_study')
        if origin is not None:
            if not isinstance(origin,dict) or set(origin)!={'study_id','revision_sha256'}:raise AnalysisError('Origin requires an exact native revision')
            s=self.studies.snapshot(origin['study_id'],revision_sha256=origin['revision_sha256'])['study'];identity=current['identity']
            if not any(x.get('reference_key')==identity['key'] and x.get('sequence_sha256')==identity['sequence_sha256'] for x in s['subjects']):raise AnalysisError('Origin study does not contain this exact reference')
        sid=request.get('study_id')
        if sid:
            s=self.studies.load(sid);f=s.get('protein_function',{})
            if f.get('schema')!=INVESTIGATION or f['key']!=current['identity']['key'] or request.get('expected_revision_sha256')!=s['revision_sha256']:raise AnalysisError('Investigation revision or protein differs; reopen before saving')
            old=self.load(sid)['packet']
            if old['function_record']['record_sha256']!=current['record_sha256']:raise AnalysisError('Saved function evidence changed; create a reviewed successor investigation')
            if old['plan'].get('origin_study')!=origin:raise AnalysisError('Keep the exact original study context')
        else:
            if request.get('expected_revision_sha256'):raise AnalysisError('Expected investigation is unavailable')
            sid='protein-function-'+secrets.token_hex(8);s={'schema_version':'1.0','contract_id':'study_record','study_id':sid,'title':('Function investigation · '+current['name'])[:160],'question':question,'scope':'Retained protein evidence and proposed investigation only; no measured-function update or scientific execution.','subjects':[{'subject_id':'protein-reference','kind':'sequence','label':current['name'][:160],'reference_key':current['identity']['key'],'sequence_sha256':current['identity']['sequence_sha256'],'artifact_refs':[]}],'artifacts':[],'evidence':[],'case_links':[],'scientifically_approved':False}
        packet={'schema':INVESTIGATION,'function_record':current,'plan':{'question':question,'action_id':action,'hypotheses':qualified,'notes':notes,'origin_study':origin,'status':'proposed','updated_at':datetime.now(timezone.utc).isoformat()},'claim':'A user investigation plan is not a function annotation, measured result or completed method.'}
        raw=_canonical(packet)
        if len(raw)>512*1024:raise AnalysisError('Function investigation exceeds the snapshot limit')
        sha=hashlib.sha256(raw).hexdigest();obj=self.studies.analysis.objects_root/sha[:2]/sha
        obj.parent.mkdir(exist_ok=True)
        if obj.exists():self.studies._object(sha)
        else:
            with obj.open('xb') as stream:stream.write(raw)
        self.studies._attach(s,s['subjects'][0],filename='protein-function-investigation-'+sha[:12]+'.json',digest=sha,size=len(raw),origin='protein_function_investigation')
        s['question']=question;s['protein_function']={'schema':INVESTIGATION,'key':current['identity']['key'],'function_record_sha256':current['record_sha256'],'packet_sha256':sha}
        return self.studies._save(s)
