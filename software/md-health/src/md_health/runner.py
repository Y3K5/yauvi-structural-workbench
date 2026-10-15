"""Explicit strict grompp plus single-process bounded pilot. No automatic repairs."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess

from .core import InputError, bound, load, sha


def pilot(prepared: Path, out: Path, *, gmx: str='gmx', timeout_seconds: int=60, threads: int=1) -> dict:
    prepared=prepared.resolve();out=out.resolve()
    if not 1<=timeout_seconds<=3600 or not 1<=threads<=64:raise InputError('Declare timeout 1..3600 seconds and threads 1..64')
    if out.exists()or out.is_relative_to(prepared)or prepared.is_relative_to(out):raise InputError('Use a new execution directory outside prepared inputs')
    record=load(prepared/'PREPARATION.json')
    if record.get('schema')!='md-health-preparation-v1':raise InputError('Unsupported preparation record')
    source=prepared/'inputs'
    for name,value in record['source_sha256'].items():
        if sha(bound(source,name))!=value:raise InputError('Prepared parent changed: '+name)
    if sha(prepared/'pilot.mdp')!=record['derived_mdp_sha256']:raise InputError('Derived MDP changed')
    values=record['pilot_steps']
    if not isinstance(values,int)or isinstance(values,bool)or not 1<=values<=1000000:raise InputError('Invalid pilot steps')
    executable=shutil.which(gmx)
    if executable is None:raise InputError('GROMACS runtime unavailable; no engine execution occurred')
    # Warning-suppression environment must not bypass the strict handoff.
    env=os.environ.copy()
    if 'GMX_MAXCONSTRWARN'in env:raise InputError('Remove GMX_MAXCONSTRWARN; this runner does not suppress constraint warnings')
    out.mkdir(parents=True)
    for name,value in record['source_sha256'].items():
        child=out/'inputs'/name;child.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(bound(source,name),child)
        if sha(child)!=value:raise InputError('Input changed during execution copy')
    shutil.copyfile(prepared/'pilot.mdp',out/'pilot.mdp')
    if sha(out/'pilot.mdp')!=record['derived_mdp_sha256']:raise InputError('Pilot MDP changed during copy')
    tpr=out/'pilot.tpr';roles=record['files']
    grompp=[executable,'grompp','-f',str(out/'pilot.mdp'),'-c',roles['coordinates'],'-p',roles['topology'],
            '-o',str(tpr),'-pp',str(out/'processed.top'),'-po',str(out/'mdout.mdp'),'-maxwarn','0']
    for role,flag in [('restraint_coordinates','-r'),('index','-n'),('checkpoint','-t'),('energy','-e')]:
        if role in roles:grompp.extend([flag,roles[role]])
    mdrun=[executable,'mdrun','-s',str(tpr),'-deffnm',str(out/'pilot'),'-nt',str(threads),'-nsteps',str(values)]
    registration={'schema':'md-health-pilot-registration-v1','utc':datetime.now(timezone.utc).isoformat(),
                  'preparation_sha256':sha(prepared/'PREPARATION.json'),'gmx_executable':str(Path(executable).resolve()),
                  'gmx_executable_sha256':sha(Path(executable)),'timeout_seconds_per_command':timeout_seconds,
                  'commands':[grompp,mdrun],'source_sha256':record['source_sha256'],'derived_mdp_sha256':sha(out/'pilot.mdp'),
                  'module_sources':{p.name:sha(p)for p in Path(__file__).parent.glob('*.py')}}
    (out/'REGISTRATION.json').write_text(json.dumps(registration,indent=2)+'\n')
    observations=[]
    for name,command in [('version',[executable,'--version']),('grompp',grompp),('mdrun',mdrun)]:
        try:
            value=subprocess.run(command,cwd=out/'inputs',env=env,capture_output=True,text=True,timeout=timeout_seconds)
            stdout,stderr,code=value.stdout,value.stderr,value.returncode
        except subprocess.TimeoutExpired as exc:
            stdout=exc.stdout or '';stderr=exc.stderr or ''
            if isinstance(stdout,bytes):stdout=stdout.decode(errors='replace')
            if isinstance(stderr,bytes):stderr=stderr.decode(errors='replace')
            code=124
        except OSError as exc:stdout='';stderr=str(exc);code=127
        (out/(name+'.stdout')).write_text(stdout);(out/(name+'.stderr')).write_text(stderr)
        observations.append({'step':name,'exit_code':code,'stdout_sha256':sha(out/(name+'.stdout')),'stderr_sha256':sha(out/(name+'.stderr'))})
        if code!=0:break
        if name=='grompp'and (not tpr.is_file()or not tpr.stat().st_size):
            observations.append({'step':'required_tpr','exit_code':2});break
    all_commands= len(observations)==3 and all(x['exit_code']==0 for x in observations)
    log=out/'pilot.log'
    result={'schema':'md-health-pilot-result-v1','status':'engine_exited_without_error'if all_commands and log.is_file()else 'failed_or_incomplete',
            'observations':observations,'registration_sha256':sha(out/'REGISTRATION.json'),
            'source_preserved':all(sha(bound(source,n))==v for n,v in record['source_sha256'].items()),
            'claim_ceiling':'Engine exit only; assess logs, metric series and trajectory integrity before interpretation.'}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    if log.is_file():
        receipt={'kind':'gromacs_mdrun_receipt_v1','exit_code':0 if all_commands else next((x['exit_code']for x in observations if x['exit_code']!=0),2),
                 'log_sha256':sha(log),'registration_sha256':sha(out/'REGISTRATION.json')}
        (out/'MDRUN_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return result


def check_trajectory(trajectory: Path, run_input: Path, out: Path, *, gmx: str='gmx', timeout_seconds: int=60) -> dict:
    """Retain gmx check output; exit zero is not an automatic integrity verdict."""
    if not 1<=timeout_seconds<=3600:raise InputError('Declare timeout 1..3600 seconds')
    if not trajectory.is_file()or not run_input.is_file()or trajectory.suffix.lower()not in ['.xtc','.trr','.tng']or run_input.suffix.lower()!='.tpr':raise InputError('Supply an existing XTC/TRR/TNG and its declared TPR')
    if out.exists():raise InputError('Trajectory-check output already exists')
    executable=shutil.which(gmx)
    if executable is None:raise InputError('GROMACS runtime unavailable; trajectory integrity is not checked')
    parents={'trajectory':sha(trajectory),'run_input':sha(run_input)}
    out.mkdir(parents=True)
    a=out/('trajectory'+trajectory.suffix.lower());b=out/'declared.tpr'
    shutil.copyfile(trajectory,a);shutil.copyfile(run_input,b)
    if sha(a)!=parents['trajectory']or sha(b)!=parents['run_input']:raise InputError('Trajectory/TPR changed during capture')
    command=[executable,'check','-f',str(a.resolve()),'-s1',str(b.resolve())]
    registration={'command':command,'input_sha256':parents,'timeout_seconds':timeout_seconds,'gmx_executable_sha256':sha(Path(executable))}
    (out/'REGISTRATION.json').write_text(json.dumps(registration,indent=2)+'\n')
    try:
        r=subprocess.run(command,cwd=out,capture_output=True,text=True,timeout=timeout_seconds)
        stdout,stderr,code=r.stdout,r.stderr,r.returncode
    except subprocess.TimeoutExpired as exc:
        stdout,stderr,code=exc.stdout or '',exc.stderr or '',124
        if isinstance(stdout,bytes):stdout=stdout.decode(errors='replace')
        if isinstance(stderr,bytes):stderr=stderr.decode(errors='replace')
    except OSError as exc:stdout,stderr,code='',str(exc),127
    (out/'check.stdout').write_text(stdout);(out/'check.stderr').write_text(stderr)
    result={'schema':'md-health-trajectory-check-v1','exit_code':code,
            'status':'review_engine_diagnostics'if code==0 else 'failed_or_incomplete',
            'inputs':parents,'registration_sha256':sha(out/'REGISTRATION.json'),
            'stdout_sha256':sha(out/'check.stdout'),'stderr_sha256':sha(out/'check.stderr'),
            'claim_ceiling':'Read the retained frame/bond diagnostic output; exit zero alone does not establish trajectory integrity or scientific convergence.'}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
