"""Optional hash-pinned local functional catalogue; no bundled research data."""
import hashlib,json,os,re
from pathlib import Path
from urllib.parse import urlsplit
from .store import AnalysisError
from .local_files import local_bytes
from .function_record_store import _text
SHA=re.compile(r'^[0-9a-f]{64}$')
ID=re.compile(r'^[a-z0-9_-]{1,120}$')
LIMIT=16*1024*1024

def digest(raw):return hashlib.sha256(raw).hexdigest()

def bounded(path,limit=LIMIT):
 p=Path(path).absolute()
 if p.resolve()!=p or not p.is_file() or p.stat().st_size>limit:raise AnalysisError('Use a bounded resident file without shortcuts')
 raw=local_bytes(p)
 if not raw or len(raw)>limit:raise AnalysisError('Local source size is invalid')
 return raw

def catalogue():
 name=os.environ.get('YAUVI_FUNCTION_RESEARCH_CONFIG')
 if not name:raise AnalysisError('Functional catalogue is not configured for this Workbench')
 config_raw=bounded(name,16384);config=json.loads(config_raw)
 if config.get('schema')!='function-research-config-v1' or not SHA.fullmatch(str(config.get('catalogue_sha256',''))):raise AnalysisError('Invalid functional catalogue configuration')
 raw=bounded(config['catalogue_path'])
 if digest(raw)!=config['catalogue_sha256']:raise AnalysisError('Functional catalogue changed; review and pin a successor before saving')
 data=json.loads(raw)
 if data.get('schema') not in {'research-decision-design-v1','functional-protein-catalog-v1'} or not SHA.fullmatch(str(data.get('dataset_id',''))):raise AnalysisError('Unsupported functional catalogue contract')
 if not isinstance(data.get('records'),list) or not 0<len(data['records'])<=100000:raise AnalysisError('Catalogue membership is invalid')
 if not SHA.fullmatch(str(data.get('source_ledger_sha256',''))):raise AnalysisError('Declare the original ledger identity')
 for name in ['groups','species']:
  registry=data.get(name)
  if not isinstance(registry,list) or not 0<len(registry)<=1000 or any(not isinstance(x,dict) or not isinstance(x.get('id'),str) or not ID.fullmatch(x['id']) or not isinstance(x.get('label'),str) or not x['label'].strip() or len(x['label'])>200 for x in registry):raise AnalysisError('Use bounded registry identifiers and labels')
  if len({x['id'] for x in registry})!=len(registry):raise AnalysisError('Duplicate catalogue registry identity')
 groups={g['id'] for g in data['groups']};species={s['id'] for s in data['species']};seen=set()
 for r in data['records']:
  key=r.get('key')
  if not isinstance(key,str) or len(key)>300 or key in seen or r.get('species') not in species or r.get('primary_group') not in groups or not SHA.fullmatch(str(r.get('sequence_sha256',''))):raise AnalysisError('Duplicate or unbound reference identity/category')
  if not all(isinstance(r.get(f),str) and r[f] for f in ['id','name','strain']):raise AnalysisError('Reference accession, annotation and strain are required')
  if type(r.get('length')) is not int or r['length']<1:raise AnalysisError('Reference length is invalid')
  e=data['evidence'].get(key)
  if not isinstance(e,dict) or not all(isinstance(e.get(f),dict) and isinstance(e[f].get('status'),str) for f in ['structure','location']):raise AnalysisError('Each member requires explicit form/location states')
  seen.add(key)
 for ref in data.get('source_references',[]):
  if urlsplit(ref.get('url','')).scheme not in {'http','https'}:raise AnalysisError('Unsafe reference URL')
 base=config.get('evidence_base_url','')
 if base:
  u=urlsplit(base)
  if u.scheme!='http' or u.hostname not in {'127.0.0.1','localhost','::1'} or u.username or u.password or u.path or u.query or u.fragment:raise AnalysisError('Evidence links must use a plain loopback origin')
 examples=config.get('example_sets',{})
 if not isinstance(examples,dict) or any(not ID.fullmatch(k) or not isinstance(v,dict) or not isinstance(v.get('label'),str) or not isinstance(v.get('keys'),list) or not 1<=len(v['keys'])<=3 or any(key not in seen for key in v['keys']) for k,v in examples.items()):raise AnalysisError('Invalid catalogue comparison sets')
 data.setdefault('profiles',{});data.setdefault('source_references',[]);data.setdefault('source_relationships',[])
 data['summary']={'total':len(data['records']),'forms':sum(e['structure'].get('available') is True for e in data['evidence'].values()),'location_gaps':sum(e['location'].get('zone')=='unplaced' for e in data['evidence'].values()),'source_cross_records':sum(e.get('a','').split('-')[0] in {s.split('-')[0] for s in species} and e.get('b','').split('-')[0] in {s.split('-')[0] for s in species} and e.get('a','').split('-')[0]!=e.get('b','').split('-')[0] for e in data['source_relationships'])}
 return data,{'example_sets':examples,'catalogue_sha256':digest(raw),'configuration_sha256':digest(config_raw),'label':_text(config.get('label','Functional research catalogue'),'Catalogue label',160),'evidence_base_url':base}
