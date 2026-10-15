from __future__ import annotations

import hashlib
import html
import json
import math
from pathlib import Path
import re
import shutil


class InputError(ValueError):
    pass


LIMIT = 'Early numerical/setup signals only. No equilibrium, convergence, force-field suitability or biological validity is established.'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    try:
        value=json.loads(path.read_text(), parse_constant=lambda s: (_ for _ in ()).throw(InputError('Non-finite JSON value '+s)))
    except (OSError, ValueError) as exc:
        raise InputError('Unreadable JSON: '+str(path)) from exc
    if not isinstance(value,dict):raise InputError('Expected a JSON object')
    return value


def bound(root: Path, name: str) -> Path:
    if not isinstance(name,str) or not name or Path(name).is_absolute() or '..' in Path(name).parts:
        raise InputError('File references must be relative within the supplied bundle')
    path=root/name
    if not path.is_file() or path.is_symlink():
        raise InputError('Missing or linked input: '+name)
    if not path.resolve().is_relative_to(root.resolve()):raise InputError('Input escapes bundle: '+name)
    # Reject linked parent directories as well as linked files.
    current=path
    while current!=root:
        if current.is_symlink():raise InputError('Linked input component: '+name)
        current=current.parent
    return path


def mdp_values(text: str) -> dict[str,str]:
    values={}
    for i,line in enumerate(text.splitlines(),1):
        line=line.split(';',1)[0].strip()
        if not line:continue
        if '=' not in line:raise InputError(f'MDP line {i} has no assignment')
        key,value=line.split('=',1);key=key.strip().lower().replace('_','-')
        if not key or key in values:raise InputError('Duplicate or empty MDP key: '+key)
        values[key]=value.strip()
    return values


def gro_check(path: Path) -> dict:
    lines=path.read_text().splitlines()
    try:
        count=int(lines[1]);assert count>0 and len(lines)==count+3
        for line in lines[2:2+count]:
            assert len(line)>=44 and line[5:10].strip() and line[10:15].strip()
            xyz=[float(line[i:i+8])for i in [20,28,36]]
            assert all(math.isfinite(x)for x in xyz)
            if line[44:].strip():
                assert len(line)>=68
                assert all(math.isfinite(float(line[i:i+8]))for i in [44,52,60])
        box=[float(x)for x in lines[-1].split()]
        assert len(box)in [3,9] and all(math.isfinite(x)for x in box)
        if len(box)==3:assert all(x>0 for x in box)
        else:
            a,b,c=box[:3];ay,az,bx,bz,cx,cy=box[3:]
            assert a*(b*c-bz*cy)-ay*(bx*c-bz*cx)+az*(bx*cy-b*cx)>0
    except (ValueError,IndexError,AssertionError) as exc:raise InputError('Invalid GRO count, atom record, finite coordinates/velocities or box')from exc
    return {'atom_count':count,'box_nm':box,'topology_atom_correspondence':'requires grompp; not checked by text audit'}


def inspect_bundle(root: Path, manifest: dict) -> dict:
    root=root.resolve()
    if manifest.get('schema')!='md-health-input-v1':raise InputError('Expected md-health-input-v1 schema')
    if manifest.get('stage')not in ['minimization','nvt','npt','production']:raise InputError('Declare the simulation stage')
    if not isinstance(manifest.get('force_field'),str)or not manifest['force_field'].strip():raise InputError('Declare the supplied force field; it is not chosen automatically')
    roles=manifest.get('files',{})
    required={'coordinates','topology','mdp'};optional={'restraint_coordinates','index','checkpoint','energy'}
    if not isinstance(roles,dict)or not required.issubset(roles)or set(roles)-required-optional:raise InputError('Declare coordinates, topology and mdp; only supported optional input roles are allowed')
    files={name:bound(root,path)for name,path in roles.items()}
    if files['coordinates'].suffix.lower()!='.gro':raise InputError('v1 requires GRO coordinates')
    gro=gro_check(files['coordinates'])
    if 'restraint_coordinates'in files:gro_check(files['restraint_coordinates'])
    settings=mdp_values(files['mdp'].read_text())
    try:
        steps=int(settings['nsteps']);assert steps>0
        integ=settings['integrator'].lower()
        if manifest['stage']=='minimization':assert integ in ['steep','cg','l-bfgs']
        else:assert integ in ['md','md-vv','md-vv-avek'];assert math.isfinite(float(settings['dt']))and float(settings['dt'])>0
    except (KeyError,ValueError,AssertionError)as exc:raise InputError('Declare finite positive dt and positive nsteps with a stage-compatible integrator')from exc
    if settings.get('include','').strip():raise InputError('External MDP include search paths are not supported; supply self-contained relative topology includes')
    refs=set(roles.values());visited=set();conditional=False
    def scan(path):
        nonlocal conditional
        name=path.relative_to(root).as_posix()
        if name in visited:return
        visited.add(name)
        for line in path.read_text().splitlines():
            text=line.split(';',1)[0].strip()
            if re.match(r'#\s*(ifdef|ifndef|else|endif|define|undef)\b',text):conditional=True
            if re.match(r'#\s*include\b',text):
                match=re.fullmatch(r'#\s*include\s*["<]([^">]+)[">]\s*',text)
                if not match:raise InputError('Unresolved topology include: '+text)
                inc=match[1]
                if Path(inc).is_absolute()or '..'in Path(inc).parts:raise InputError('Topology include escapes its declared bundle')
                candidates=[path.parent/inc,root/inc];target=next((p for p in candidates if p.is_file()),None)
                if target is None:raise InputError('Missing topology dependency (including conditional branches): '+inc)
                safe=bound(root,target.relative_to(root).as_posix());refs.add(safe.relative_to(root).as_posix());scan(safe)
    scan(files['topology'])
    parents={name:sha(bound(root,name))for name in sorted(refs)}
    return {'status':'text_audit_only','files':roles,'source_sha256':parents,'settings':settings,
            'stage':manifest['stage'],'force_field':manifest['force_field'],'gro':gro,
            'conditional_topology':conditional,'grompp_required':True,'physics_unchanged':True,'claim_ceiling':LIMIT}


def prepare(root: Path, manifest: dict, out: Path, pilot_steps: int) -> dict:
    if not isinstance(pilot_steps,int)or isinstance(pilot_steps,bool)or not 1<=pilot_steps<=1000000:raise InputError('Pilot steps must be 1..1000000 (execution cap, not a recommended scientific duration)')
    root=root.resolve();out=out.resolve()
    if out.exists()or out.is_relative_to(root)or root.is_relative_to(out):raise InputError('Use a new output directory outside the parent bundle')
    result=inspect_bundle(root,manifest)
    out.mkdir(parents=True)
    for name,value in result['source_sha256'].items():
        source=bound(root,name)
        if sha(source)!=value:raise InputError('Source changed during preparation')
        child=out/'inputs'/name;child.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,child)
        if sha(child)!=value:raise InputError('Copy identity mismatch')
    parent_mdp=out/'inputs'/result['files']['mdp']
    text=parent_mdp.read_text()
    child=re.sub(r'(?im)^\s*nsteps\s*=.*$',f'nsteps = {pilot_steps} ; MD Health explicit bounded pilot override',text)
    if mdp_values(child)['nsteps']!=str(pilot_steps):raise InputError('Could not derive the bounded pilot MDP')
    (out/'pilot.mdp').write_text(child)
    record={**result,'schema':'md-health-preparation-v1','input_manifest':manifest,
            'pilot_steps':pilot_steps,'changes':[{'field':'nsteps','parent':result['settings']['nsteps'],'child':str(pilot_steps)}],
            'derived_mdp_sha256':sha(out/'pilot.mdp'),'source_preserved':all(sha(bound(root,n))==v for n,v in result['source_sha256'].items())}
    (out/'PREPARATION.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


def read_xvg(path: Path) -> dict:
    labels={};rows=[];units='unknown'
    for number,line in enumerate(path.read_text().splitlines(),1):
        s=line.strip()
        if not s or s.startswith('#'):continue
        if s.startswith('@'):
            m=re.match(r'@\s*s(\d+)\s+legend\s+"([^"]+)"',s)
            if m:labels[int(m[1])]=m[2]
            if re.fullmatch(r'@\s*xaxis\s+label\s+"Time \(ps\)"',s):units='ps'
            continue
        if s=='&':raise InputError('Multiple XVG datasets need separate declared inputs')
        try:values=[float(x)for x in s.split()]
        except ValueError as exc:raise InputError(f'Unreadable XVG row {number}')from exc
        if len(values)<2 or (rows and len(values)!=len(rows[0]['values'])):raise InputError('Inconsistent XVG column count')
        rows.append({'line':number,'values':values})
    if labels and rows and (set(labels)!=set(range(len(rows[0]['values'])-1))or len(set(labels.values()))!=len(labels)):raise InputError('Incomplete/duplicate XVG metric labels')
    return {'labels':labels,'rows':rows,'time_unit':units}


def assess(log: Path, series: Path|None=None, policy: dict|None=None, run_receipt: Path|None=None) -> dict:
    policy=policy or {};findings=[];missing=[];metrics={};inputs={'log':{'name':log.name,'sha256':sha(log)}}
    def flag(code,severity,message,locator):findings.append({'code':code,'severity':severity,'message':message,'locator':locator})
    text=log.read_text(errors='replace');lines=text.splitlines()
    if sha(log)!=inputs['log']['sha256']:raise InputError('Log changed during reading; assess a captured copy')
    for i,line in enumerate(lines,1):
        if re.search(r'Fatal error:|Segmentation fault|cannot be settled',line,re.I):flag('fatal_log','failure',line.strip(),{'file':log.name,'line':i})
        elif re.search(r'LINCS WARNING|SHAKE.*(?:warning|error)|SETTLE.*(?:warning|error)|constraint.*(?:warning|error)',line,re.I):flag('constraint_warning','review',line.strip(),{'file':log.name,'line':i})
        elif re.match(r'\s*WARNING\b',line):flag('engine_warning','review',line.strip(),{'file':log.name,'line':i})
    if not lines:missing.append('nonempty_log')
    if not re.search(r'GROMACS|GROMACS version',text,re.I):missing.append('GROMACS_log_identity')
    finished=bool(re.search(r'^\s*Finished mdrun\b',text,re.M))
    if not finished:missing.append('finished_mdrun_marker (may still be running or stopped)')
    if run_receipt:
        receipt=load(run_receipt);inputs['run_receipt']={'name':run_receipt.name,'sha256':sha(run_receipt)}
        if receipt.get('exit_code')!=0:flag('execution_exit','failure','Recorded execution did not exit successfully',{'file':run_receipt.name,'pointer':'/exit_code'})
        if receipt.get('log_sha256')!=sha(log):raise InputError('Run receipt is not bound to this log')
        if receipt.get('kind')!='gromacs_mdrun_receipt_v1':raise InputError('Unsupported run receipt')
    else:missing.append('log-bound execution exit receipt')
    if series:
        inputs['series']={'name':series.name,'sha256':sha(series)}
        data=read_xvg(series);rows=data['rows']
        if sha(series)!=inputs['series']['sha256']:raise InputError('Series changed during reading; assess a captured copy')
        if len(rows)<2:missing.append('at least two metric samples')
        for row in rows:
            if any(not math.isfinite(v)for v in row['values']):flag('nonfinite_series','failure','Non-finite value in recorded series',{'file':series.name,'line':row['line']})
        times=[r['values'][0]for r in rows]
        if any(b<=a for a,b in zip(times,times[1:])):flag('nonincreasing_time','review','Times repeat or reverse; inspect concatenation/restart history',{'file':series.name})
        if data['time_unit']!='ps':missing.append('declared XVG time units in ps')
        end=policy.get('expected_end_ps');tol=policy.get('time_tolerance_ps',0)
        if end is None:missing.append('expected_end_ps')
        elif not isinstance(end,(int,float))or isinstance(end,bool)or not math.isfinite(end)or end<0:raise InputError('expected_end_ps must be finite and nonnegative')
        elif not isinstance(tol,(int,float))or isinstance(tol,bool)or not math.isfinite(tol)or tol<0:raise InputError('time_tolerance_ps must be finite and nonnegative')
        elif not times or not math.isfinite(times[-1])or times[-1]+tol<end:flag('short_series','review','Recorded series stops before the declared endpoint; capture may be incomplete',{'file':series.name})
        required=policy.get('required_metrics',[]);bands=policy.get('bands',{})
        if not isinstance(bands,dict):raise InputError('bands must be an object')
        if not isinstance(required,list)or not required or any(not isinstance(x,str)for x in required):missing.append('nonempty required_metrics declaration')
        if not isinstance(required,list)or any(not isinstance(x,str)for x in required):raise InputError('required_metrics must be a list of names')
        for label in required:
            if label not in data['labels'].values():missing.append('metric:'+label)
            if label not in bands:missing.append('declared band:'+label)
        for col,label in data['labels'].items():
            finite=[r['values'][col+1]for r in rows if math.isfinite(r['values'][col+1])]
            metrics[label]={'samples':len(rows),'finite_samples':len(finite),'min':min(finite)if finite else None,'max':max(finite)if finite else None,'unit':'as declared in input; verify GROMACS units'}
            if label in bands:
                band=bands[label]
                if not isinstance(band,list)or len(band)!=2 or any(not isinstance(x,(int,float))or isinstance(x,bool)or not math.isfinite(x)for x in band)or band[0]>band[1]:raise InputError('Declare ordered finite bounds for '+label)
                for row in rows:
                    value=row['values'][col+1]
                    if math.isfinite(value)and not band[0]<=value<=band[1]:flag('declared_band_departure','review',f'{label}: {value} outside user-declared {band}',{'file':series.name,'line':row['line'],'metric':label})
    else:missing.append('named metric series')
    severity={f['severity']for f in findings}
    status='failure_signal_observed'if 'failure'in severity else 'review_required'if 'review'in severity else 'insufficient_evidence'if missing else 'no_failure_signal_observed'
    return {'schema':'md-health-report-v1','status':status,'findings':findings,'missing_evidence':missing,
            'metrics':metrics,'inputs':inputs,'policy':policy,'finished_marker':finished,'claim_ceiling':LIMIT,
            'trajectory_integrity':'not checked; run gmx check with recorded outputs separately','scientific_qualification':'not_established'}


def write_report(report: dict, out: Path):
    if out.exists():raise InputError('Report output already exists; preserve it and choose a new directory')
    out.mkdir(parents=True)
    (out/'REPORT.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    e=html.escape
    items=''.join('<li><strong>'+e(f['code'])+'</strong>: '+e(f['message'])+'<br>Source: '+e(json.dumps(f['locator']))+'</li>'for f in report['findings'])
    gaps=''.join('<li>'+e(x)+'</li>'for x in report['missing_evidence'])
    cards=''.join('<tr><td>'+e(name)+'</td><td>'+str(value['samples'])+'</td><td>'+str(value['min'])+'</td><td>'+str(value['max'])+'</td></tr>'for name,value in report['metrics'].items())
    page='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MD Health · early signals</title><style>body{max-width:960px;margin:40px auto;padding:0 20px;font:17px/1.6 system-ui;color:#20322e;background:#f5f7f5}section{background:white;padding:24px;border:1px solid #d7e1dc;border-radius:12px;margin:20px 0}h1{line-height:1.2}table{width:100%;border-collapse:collapse}th,td{text-align:left;padding:10px;border-bottom:1px solid #d7e1dc}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px}.status{font-size:22px;font-weight:bold}.boundary{border-left:4px solid #286b65;padding-left:16px}</style><p>YAUVI · GROMACS MD HEALTH · PROTOTYPE</p><h1>Check the run before interpreting the motion.</h1>'''
    page+='<section><p class="status">'+e(report['status'].replace('_',' '))+'</p><p class="boundary">'+e(LIMIT)+'</p></section><section><h2>Observed signals</h2><ul>'+items+'</ul>'+('<p>No configured signal was detected in the supplied evidence.</p>'if not items else '')+'</section><section><h2>Evidence still needed</h2><ul>'+gaps+'</ul><p>'+e(report['trajectory_integrity'])+'</p></section><section><h2>Recorded series</h2><p>Bands are explicitly user-declared. Instantaneous pressure fluctuations and RMSD alone are not treated as universal failure criteria.</p><table><tr><th>Metric</th><th>Samples</th><th>Minimum</th><th>Maximum</th></tr>'+cards+'</table></section><details><summary>Exact report, input fingerprints and policy</summary><pre>'+e(json.dumps(report,indent=2,allow_nan=False))+'</pre></details></html>'
    (out/'REPORT.html').write_text(page)
    (out/'CHECKSUMS.json').write_text(json.dumps({p.name:sha(p)for p in [out/'REPORT.json',out/'REPORT.html']},indent=2)+'\n')
