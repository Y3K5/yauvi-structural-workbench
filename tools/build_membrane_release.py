#!/usr/bin/env python3
"""Build local code/metadata-only candidates and inspect every archive member.

No installation, network, Git, publication-checkout sync or scientific execution.
The manifest is an explicit review list, never permission to release its files.
"""
import argparse
import getpass
import fnmatch
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
try:import tomllib
except ImportError:import tomli as tomllib


def sha(raw):return hashlib.sha256(raw).hexdigest()


def snapshot(source,target,standalone=False):
    config=tomllib.loads((source/'pyproject.toml').read_text());setting=config['tool']['setuptools']
    paths={'pyproject.toml','MANIFEST.in'}
    paths.add(config['project']['readme'])
    paths.update(setting.get('license-files',[]))
    if standalone:paths.add('LICENSE')
    else:paths.update(['NOTICE.md','tools/build_membrane_release.py','tools/verify_structural_workbench_wheel.py','tools/test_workbench_wheel.py'])
    find=setting['packages']['find'];patterns=find.get('include',['memorient*'])
    for directory in find['where']:
        root=source/directory
        for path in root.rglob('*.py'):
            relative=path.relative_to(root);module='.'.join(relative.parts[:-1])
            if any(part.startswith(('_archive','__pycache__')) for part in relative.parts):continue
            if any(fnmatch.fnmatch(module,p) for p in patterns):paths.add(path.relative_to(source).as_posix())
        for package,assets in setting.get('package-data',{}).items():
            package_root=root/Path(*package.split('.'))
            if not package_root.is_dir():continue
            for pattern in assets:
                for path in package_root.glob(pattern):
                    if path.is_file():paths.add(path.relative_to(source).as_posix())
    # New regression controls and public guide; no biological test fixtures.
    extra=(['tests/test_run_bundle.py','tests/test_md_geometry.py'] if standalone else [
        'software/platform/tests/test_membrane_workflow.py','software/platform/tests/test_membrane_process_control.py','software/membrane-orientor/memorient/tests/test_run_bundle.py',
        'software/membrane-orientor/memorient/tests/test_md_geometry.py','software/structural-workbench/tests/test_membrane_server.py',
        'software/structural-workbench/tests/conftest.py','tools/test_membrane_release.py',
        'software/membrane-orientor/memorient/pyproject.toml','software/membrane-orientor/memorient/MANIFEST.in',
        'software/membrane-orientor/memorient/README.md','software/membrane-orientor/memorient/LICENSE',
        '.github/workflows/structural-workbench.yml'])
    paths.update(extra)
    guide='docs/OPEN_SOURCE.md' if standalone else 'software/membrane-orientor/memorient/docs/OPEN_SOURCE.md'
    if (source/guide).is_file():paths.add(guide)
    runtime='requirements-md-prep.txt' if standalone else 'software/membrane-orientor/memorient/requirements-md-prep.txt'
    if (source/runtime).is_file():paths.add(runtime)
    selected=[]
    for name in sorted(paths):
        path=source/name
        if path.is_symlink() or not path.is_file():raise ValueError('Candidate input is missing or a shortcut: '+name)
        raw=path.read_bytes();out=target/name;out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(raw)
        selected.append({'path':name,'bytes':len(raw),'sha256':sha(raw)})
    return selected


RULES={
    'personal_absolute_path':re.compile(rb'(?:/(?:Users|home)/[A-Za-z0-9._-]+/[^\s"\x27<>]*|\$HOME/[A-Za-z0-9._-][^\s"\x27<>]*)',re.I),
    'private_key':re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
    'credential_shape':re.compile(rb'(?:gh[pousr]_[A-Za-z0-9]{20,}|sk-(?:proj-)?[A-Za-z0-9_-]{25,}|AKIA[A-Z0-9]{16})'),
    'credential_assignment':re.compile(rb'(?:api_key|password|token|secret)\s*[=:]\s*["\x27][A-Za-z0-9_\-/+]{20,}["\x27]',re.I),
    'local_account_name':re.compile(re.escape(getpass.getuser().encode()),re.I),
}


def inspect_archive(path):
    rows=[];findings=[]
    if path.suffix=='.whl':
        with zipfile.ZipFile(path) as archive:members=[(n,archive.read(n)) for n in archive.namelist() if not n.endswith('/')]
    else:
        with tarfile.open(path) as archive:members=[(r.name,archive.extractfile(r).read()) for r in archive.getmembers() if r.isfile()]
    for name,raw in members:
        rows.append({'path':name,'bytes':len(raw),'sha256':sha(raw)})
        if any(marker in name for marker in ['/ecoli_om/','/catalogue/sources/','/md-prep/','/model-v2/','/minimized-view','catalogue/starter.json']):
            findings.append({'member':name,'rule':'excluded_biological_payload','count':1})
        for rule,pattern in RULES.items():
            count=len(pattern.findall(raw))
            if count:findings.append({'member':name,'rule':rule,'count':count})
    return {'artifact':path.name,'sha256':sha(path.read_bytes()),'bytes':path.stat().st_size,'members':rows,'findings':findings}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',required=True)
    parser.add_argument('--private-terms',help='local-only text file of additional private terms; never copied to candidates')
    args=parser.parse_args()
    if args.private_terms:
        for index,term in enumerate(Path(args.private_terms).read_text().splitlines()):
            if term.strip():RULES['private_term_'+str(index+1)]=re.compile(re.escape(term.strip().encode()),re.I)
    home=Path(__file__).resolve().parents[1];output=Path(args.out).resolve();output.mkdir(parents=True,exist_ok=False)
    cases=[('workbench',home,False),('memorient',home/'software/membrane-orientor/memorient',True)]
    report={'schema_version':'memorient.release_candidate.v1','publication':'not_authorized','source_snapshots':{},'archives':[]}
    # Reuse setuptools' already bundled wheel if the host has no standalone
    # wheel distribution. No package download or installation occurs.
    env=os.environ.copy()
    import setuptools
    vendor=Path(setuptools.__file__).parent/'_vendor'
    if (vendor/'wheel').is_dir():
        env['PYTHONPATH']=str(vendor)+os.pathsep+env.get('PYTHONPATH','')
    for name,source,standalone in cases:
        target=output/name;target.mkdir()
        report['source_snapshots'][name]=snapshot(source,target,standalone)
        log=output/(name+'-build.log')
        with log.open('w') as stream:
            subprocess.run([sys.executable,'-m','build','--no-isolation','--outdir',str(output/'dist'),str(target)],env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
    for path in sorted((output/'dist').iterdir()):report['archives'].append(inspect_archive(path))
    report['status']='HOLD' if any(r['findings'] for r in report['archives']) else 'READY_FOR_CONTENT_REVIEW'
    (output/'RELEASE_MANIFEST.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'archives':len(report['archives']),
                      'findings':sum(len(r['findings']) for r in report['archives'])},indent=2))


if __name__=='__main__':main()
