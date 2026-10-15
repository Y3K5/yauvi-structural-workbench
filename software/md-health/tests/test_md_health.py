import json
from pathlib import Path
import subprocess
import sys

import pytest

from md_health.core import InputError, assess, inspect_bundle, load, prepare, sha, write_report
from md_health.runner import pilot, check_trajectory


@pytest.fixture
def bundle(tmp_path):
    root=tmp_path/'parent';root.mkdir()
    (root/'system.gro').write_text('Invented single-atom software fixture; not a physical MD system\n1\n'+f'{1:5d}{"TST":<5}{"X":>5}{1:5d}{0.1:8.3f}{0.2:8.3f}{0.3:8.3f}'+'\n1.0 1.0 1.0\n')
    (root/'system.top').write_text('#include "toy.itp"\n[ system ]\nInvented software fixture\n')
    (root/'toy.itp').write_text('; Invented include; no force-field parameters or runnable MD\n')
    (root/'run.mdp').write_text('integrator = md\ndt = 0.002\nnsteps = 1000\ntcoupl = no\npcoupl = no\n')
    manifest={'schema':'md-health-input-v1','stage':'nvt','force_field':'invented-parser-fixture-only',
              'files':{'coordinates':'system.gro','topology':'system.top','mdp':'run.mdp'}}
    return root,manifest


def evidence(tmp_path, *, log='GROMACS version fixture\nFinished mdrun on test\n', values='0 300\n1 300\n'):
    p=tmp_path/'run.log';p.write_text(log)
    x=tmp_path/'energy.xvg';x.write_text('@ xaxis label "Time (ps)"\n@ s0 legend "Temperature"\n'+values)
    r=tmp_path/'receipt.json';r.write_text(json.dumps({'kind':'gromacs_mdrun_receipt_v1','exit_code':0,'log_sha256':sha(p)}))
    policy={'expected_end_ps':1,'required_metrics':['Temperature'],'bands':{'Temperature':[290,310]}}
    return p,x,r,policy


def test_derivation_preserves_parent_and_physics(bundle,tmp_path):
    root,manifest=bundle;parent={p.name:sha(p)for p in root.iterdir()}
    record=prepare(root,manifest,tmp_path/'child',40)
    assert record['changes']==[{'field':'nsteps','parent':'1000','child':'40'}]
    assert record['source_preserved']
    assert {p.name:sha(p)for p in root.iterdir()}==parent
    assert 'dt = 0.002' in (tmp_path/'child/pilot.mdp').read_text()
    assert sha(tmp_path/'child/inputs/run.mdp')==parent['run.mdp']
    assert set(record['source_sha256'])==set(parent)


@pytest.mark.parametrize('name',['missing','escape','symlink','macro','external_mdp','duplicate_mdp','nan_gro','wrong_count','negative_box','unbounded_steps','bad_dt','wrong_stage'])
def test_input_refusals(bundle,tmp_path,name):
    root,manifest=bundle
    if name=='missing':(root/'toy.itp').unlink()
    elif name=='escape':(root/'system.top').write_text('#include "../outside.itp"\n')
    elif name=='symlink':(root/'toy.itp').unlink();(root/'toy.itp').symlink_to(root/'run.mdp')
    elif name=='macro':(root/'system.top').write_text('#include FILE\n')
    elif name=='external_mdp':(root/'run.mdp').write_text((root/'run.mdp').read_text()+'include = -I/private\n')
    elif name=='duplicate_mdp':(root/'run.mdp').write_text((root/'run.mdp').read_text()+'nsteps = 20\n')
    elif name=='nan_gro':(root/'system.gro').write_text((root/'system.gro').read_text().replace('   0.100','     nan'))
    elif name=='wrong_count':(root/'system.gro').write_text((root/'system.gro').read_text().replace('\n1\n','\n2\n'))
    elif name=='negative_box':(root/'system.gro').write_text((root/'system.gro').read_text().replace('1.0 1.0 1.0','-1.0 1.0 1.0'))
    elif name=='unbounded_steps':(root/'run.mdp').write_text((root/'run.mdp').read_text().replace('1000','-1'))
    elif name=='bad_dt':(root/'run.mdp').write_text((root/'run.mdp').read_text().replace('0.002','nan'))
    elif name=='wrong_stage':manifest['stage']='minimization'
    with pytest.raises(InputError):inspect_bundle(root,manifest)


def test_conditional_dependencies_are_declared_and_not_silently_skipped(bundle):
    root,manifest=bundle;(root/'system.top').write_text('#ifdef POSRES\n#include "toy.itp"\n#endif\n')
    r=inspect_bundle(root,manifest);assert r['conditional_topology'] and r['grompp_required']
    (root/'toy.itp').unlink()
    with pytest.raises(InputError):inspect_bundle(root,manifest)


@pytest.mark.parametrize('nested',['parent','ancestor','existing'])
def test_output_never_overwrites_or_nests_parent(bundle,tmp_path,nested):
    root,manifest=bundle
    out=root/'child'if nested=='parent'else tmp_path if nested=='ancestor'else root
    with pytest.raises(InputError):prepare(root,manifest,out,10)


@pytest.mark.parametrize('log,code,severity',[
    ('GROMACS\nFatal error:\nToo many constraint warnings','fatal_log','failure'),
    ('GROMACS\nStep 20, LINCS WARNING','constraint_warning','review'),
    ('GROMACS\nSHAKE error','constraint_warning','review'),
    ('GROMACS\nWater molecule cannot be settled','fatal_log','failure'),
    ('GROMACS\nWARNING 1: incomplete topology','engine_warning','review'),
])
def test_concrete_log_signals(tmp_path,log,code,severity):
    p,x,r,policy=evidence(tmp_path,log=log);report=assess(p,x,policy,r)
    assert any(f['code']==code and f['severity']==severity and f['locator']['line']>0 for f in report['findings'])


def test_lincs_settings_are_not_a_failure_warning(tmp_path):
    p,x,r,policy=evidence(tmp_path,log='GROMACS\nLINCS: order 4\nFinished mdrun on test\n')
    assert assess(p,x,policy,r)['status']=='no_failure_signal_observed'


@pytest.mark.parametrize('value',['nan','inf','-inf'])
def test_nonfinite_metrics_fail_and_report_is_strict_json(tmp_path,value):
    p,x,r,policy=evidence(tmp_path,values=f'0 300\n1 {value}\n')
    report=assess(p,x,policy,r);assert report['status']=='failure_signal_observed'
    write_report(report,tmp_path/'report');json.loads((tmp_path/'report/REPORT.json').read_text())
    assert 'NaN'not in (tmp_path/'report/REPORT.json').read_text()


def test_user_band_departure_is_review_not_universal_failure(tmp_path):
    p,x,r,policy=evidence(tmp_path,values='0 300\n1 340\n')
    report=assess(p,x,policy,r)
    assert report['status']=='review_required'
    assert report['findings'][0]['code']=='declared_band_departure'


@pytest.mark.parametrize('what',['empty','receipt','finish','metric','policy','sample','units','short','reverse'])
def test_missing_or_partial_evidence_cannot_pass(tmp_path,what):
    p,x,r,policy=evidence(tmp_path)
    if what=='empty':p.write_text('');r.write_text(json.dumps({'kind':'gromacs_mdrun_receipt_v1','exit_code':0,'log_sha256':sha(p)}))
    if what=='receipt':r=None
    if what=='finish':p.write_text('GROMACS\nStep 0\n');r.write_text(json.dumps({'kind':'gromacs_mdrun_receipt_v1','exit_code':0,'log_sha256':sha(p)}))
    if what=='metric':policy['required_metrics'].append('Density')
    if what=='policy':policy={}
    if what=='sample':x.write_text('@ xaxis label "Time (ps)"\n@ s0 legend "Temperature"\n1 300\n')
    if what=='units':x.write_text(x.read_text().replace('Time (ps)','Time (ns)'))
    if what=='short':policy['expected_end_ps']=5
    if what=='reverse':x.write_text(x.read_text().replace('0 300\n1 300','1 300\n0 300'))
    assert assess(p,x,policy,r)['status']in ['insufficient_evidence','review_required']


@pytest.mark.parametrize('bad',['ragged','labels','datasets'])
def test_corrupt_series_refused(tmp_path,bad):
    p,x,r,policy=evidence(tmp_path)
    if bad=='ragged':x.write_text(x.read_text()+'2 300 2\n')
    if bad=='labels':x.write_text(x.read_text().replace('s0','s1'))
    if bad=='datasets':x.write_text(x.read_text()+'&\n')
    with pytest.raises(InputError):assess(p,x,policy,r)


def test_mismatched_exit_receipt_and_failed_exit(tmp_path):
    p,x,r,policy=evidence(tmp_path);d=json.loads(r.read_text());d['log_sha256']='wrong';r.write_text(json.dumps(d))
    with pytest.raises(InputError):assess(p,x,policy,r)
    d['log_sha256']=sha(p);d['exit_code']=1;r.write_text(json.dumps(d))
    assert assess(p,x,policy,r)['status']=='failure_signal_observed'


def test_no_runtime_is_explicit_and_creates_no_execution(bundle,tmp_path,monkeypatch):
    root,m=bundle;prepare(root,m,tmp_path/'child',10)
    monkeypatch.setattr('md_health.runner.shutil.which',lambda _:None)
    with pytest.raises(InputError,match='no engine execution'):pilot(tmp_path/'child',tmp_path/'run')
    assert not (tmp_path/'run').exists()


@pytest.mark.parametrize('failure',['grompp','false_tpr','mdrun','timeout','changed_input','suppression'])
def test_runner_failure_is_retained_and_never_promoted(bundle,tmp_path,monkeypatch,failure):
    root,m=bundle;prepare(root,m,tmp_path/'child',10)
    monkeypatch.setattr('md_health.runner.shutil.which',lambda _:sys.executable)
    calls=[]
    def fake(command,**kwargs):
        calls.append(command)
        if failure=='timeout':raise subprocess.TimeoutExpired(command,kwargs['timeout'],output=b'partial')
        if 'grompp'in command:
            assert command[command.index('-maxwarn')+1]=='0'
            if failure=='grompp':return subprocess.CompletedProcess(command,1,'','WARNING 1')
            if failure!='false_tpr':Path(command[command.index('-o')+1]).write_bytes(b'mock-tpr')
        if 'mdrun'in command:
            assert command[command.index('-nsteps')+1]=='10'
            Path(str(tmp_path/'run/pilot')+'.log').write_text('GROMACS mock\n')
            return subprocess.CompletedProcess(command,1,'','mock fatal')
        return subprocess.CompletedProcess(command,0,'GROMACS mocked version','')
    monkeypatch.setattr('md_health.runner.subprocess.run',fake)
    if failure=='changed_input':(tmp_path/'child/inputs/toy.itp').write_text('changed')
    if failure=='suppression':monkeypatch.setenv('GMX_MAXCONSTRWARN','-1')
    if failure in ['changed_input','suppression']:
        with pytest.raises(InputError):pilot(tmp_path/'child',tmp_path/'run')
        assert not calls
    else:
        result=pilot(tmp_path/'child',tmp_path/'run')
        assert result['status']=='failed_or_incomplete' and result['source_preserved']
        assert (tmp_path/'run/RESULT.json').is_file()
        if failure in ['grompp','false_tpr']:assert not any('mdrun'in c for c in calls)


def test_no_report_overwrite_and_html_escapes_supplied_text(tmp_path):
    p,x,r,policy=evidence(tmp_path,log='GROMACS\nFatal error: <script>bad</script>\n')
    report=assess(p,x,policy,r);write_report(report,tmp_path/'report')
    assert '<script>bad</script>'not in (tmp_path/'report/REPORT.html').read_text()
    with pytest.raises(InputError):write_report(report,tmp_path/'report')


def test_json_policy_nonfinite_is_rejected(tmp_path):
    p=tmp_path/'bad.json';p.write_text('{"expected_end_ps":NaN}')
    with pytest.raises(InputError):load(p)


@pytest.mark.parametrize('code',[0,1,124])
def test_trajectory_engine_exit_is_not_promoted_to_qualification(tmp_path,monkeypatch,code):
    a=tmp_path/'test.xtc';a.write_bytes(b'invented binary placeholder, not a trajectory')
    b=tmp_path/'test.tpr';b.write_bytes(b'invented placeholder')
    monkeypatch.setattr('md_health.runner.shutil.which',lambda _:sys.executable)
    def fake(command,**kwargs):
        assert command[1]=='check'and '-s1'in command
        if code==124:raise subprocess.TimeoutExpired(command,kwargs['timeout'],output=b'partial')
        return subprocess.CompletedProcess(command,code,'mock diagnostic text','')
    monkeypatch.setattr('md_health.runner.subprocess.run',fake)
    result=check_trajectory(a,b,tmp_path/'check')
    assert result['exit_code']==code
    assert result['status']==('review_engine_diagnostics'if code==0 else 'failed_or_incomplete')
    assert (tmp_path/'check/check.stdout').is_file()
    assert sha(a)==result['inputs']['trajectory']


def test_execution_copy_uses_only_declared_members(bundle,tmp_path,monkeypatch):
    root,m=bundle;prepare(root,m,tmp_path/'child',10)
    (tmp_path/'child/inputs/unrelated.private').write_text('must not be copied')
    monkeypatch.setattr('md_health.runner.shutil.which',lambda _:sys.executable)
    monkeypatch.setattr('md_health.runner.subprocess.run',lambda command,**kwargs:subprocess.CompletedProcess(command,1,'mock',''))
    result=pilot(tmp_path/'child',tmp_path/'run')
    assert not (tmp_path/'run/inputs/unrelated.private').exists()
    assert result['status']=='failed_or_incomplete'
