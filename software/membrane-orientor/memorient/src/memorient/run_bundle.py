"""Portable, local result bundles. No execution, fetching or curvature remapping.

Artifact paths are relative to the manifest and every artifact is checksum-bound.
Original scientific receipts remain byte-identical inside an imported bundle.
"""
from pathlib import Path, PurePosixPath
import hashlib
import json
import math
import shutil

SCHEMA='memorient.run.v1'
STAGES={'synthetic','unrelaxed','prepared','minimized','failed','held','interrupted'}


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def artifact_path(root, relative):
    if not isinstance(relative,str) or '\\' in relative:
        raise ValueError('Invalid artifact path')
    if any(ord(c)<32 or ord(c)==127 or c=='\"' for c in relative):
        raise ValueError('Unsafe artifact path characters')
    p=PurePosixPath(relative)
    if p.is_absolute() or '..' in p.parts or ':' in relative or not p.parts:
        raise ValueError('Artifact path must stay inside its bundle')
    root=Path(root).resolve();target=root/relative
    if not target.resolve().is_relative_to(root) or target.is_symlink():
        raise ValueError('Artifact path escapes bundle')
    return target


def load_manifest(path):
    path=Path(path);doc=json.loads(path.read_text())
    if doc.get('schema_version')!=SCHEMA or doc.get('stage') not in STAGES:
        raise ValueError('Unknown membrane run schema or stage')
    if doc.get('units') not in {'angstrom','nm'} or not isinstance(doc.get('artifacts'),dict):
        raise ValueError('Run needs units and checksum-bound artifacts')
    if doc.get('equilibration')!='not_started':
        raise ValueError('Equilibration is outside this workflow')
    for name,record in doc['artifacts'].items():
        target=artifact_path(path.parent,record['path'])
        if not target.is_file() or digest(target)!=record['sha256']:
            raise ValueError('Run artifact checksum mismatch: '+name)
    if doc['stage']=='minimized':
        validate_minimum(doc)
        for name in ('coordinates','topology','minimization_log'):
            if name not in doc['artifacts']:raise ValueError('Minimized result needs '+name)
        from .md_preparation import minimization_result
        log=artifact_path(path.parent,doc['artifacts']['minimization_log']['path']).read_text()
        if minimization_result(log)!=doc['minimization']:raise ValueError('Minimization log does not match the accepted result')
        checks=doc['geometry_checks'].get('original',doc['geometry_checks'])
        if (checks.get('coordinate_sha256')!=doc['artifacts']['coordinates']['sha256']
            or checks.get('topology_sha256')!=doc['artifacts']['topology']['sha256']):
            raise ValueError('Geometry checks are not bound to these coordinates and topology')
    return doc


def validate_minimum(doc):
    m=doc.get('minimization',{})
    force=m.get('maximum_force_kj_mol_nm');energy=m.get('potential_energy_kj_mol');charge=doc.get('charge_e')
    if (m.get('converged') is not True or m.get('constraint_or_fatal_messages')
        or type(force) not in (int,float) or type(energy) not in (int,float)
        or not math.isfinite(force) or not math.isfinite(energy) or not 0<=force<1000
        or type(charge) not in (int,float) or not math.isfinite(charge) or abs(charge)>1e-4
        or doc.get('geometry_checks',{}).get('passed') is not True):
        raise ValueError('Minimized result has not passed acceptance checks')


def write_manifest(root, *,stage,artifacts,units='angstrom',**metadata):
    root=Path(root)
    doc={**metadata,'schema_version':SCHEMA,'stage':stage,'units':units,
         'equilibration':'not_started','artifacts':{}}
    for name,relative in artifacts.items():
        target=artifact_path(root,relative)
        doc['artifacts'][name]={'path':relative,'sha256':digest(target),'bytes':target.stat().st_size}
    if stage=='minimized':validate_minimum(doc)
    out=root/'MEMBRANE_RUN.json'
    out.write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n')
    load_manifest(out)
    return out


def scene(path):
    doc=load_manifest(path)
    record=doc['artifacts'].get('scene')
    if record is None:raise ValueError('This run has no display artifact')
    result=json.loads(artifact_path(Path(path).parent,record['path']).read_text())
    return {**result,'run_metadata':{k:v for k,v in doc.items() if k!='artifacts'}}


def import_legacy(run,display,destination):
    """Copy an accepted legacy result without modifying its original receipts."""
    from .md_scene import minimized_scene
    from .md_preparation import read_minimization
    run=Path(run).resolve();display=Path(display).resolve()
    minimized_scene(run=run,display=display)
    receipt_path=run/'MINIMIZATION_CONTINUATION.json'
    if not receipt_path.is_file():receipt_path=run/'PREPARATION_RECEIPT.json'
    receipt=json.loads(receipt_path.read_text())
    checks=json.loads((run/'FINAL_GEOMETRY_CHECK.json').read_text())
    # A legacy check must explicitly declare acceptance, never merely exist.
    passed=checks.get('passed') is True or checks.get('state')=='passed'
    if not passed:raise ValueError('Legacy geometry checks are not accepted')
    if checks.get('coordinate_sha256')!=digest(run/'em.gro') or checks.get('topology_sha256')!=digest(run/'system.top'):
        raise ValueError('Legacy geometry checks belong to different coordinates')
    root=Path(destination);root.mkdir(parents=True,exist_ok=False)
    files={'coordinates':run/'em.gro','topology':run/'system.top','scene':display/'scene.json',
           'atom_map':display/'ATOM_MAP.tsv','receipt':receipt_path,'geometry':run/'FINAL_GEOMETRY_CHECK.json',
           'minimization_log':run/'em.log'}
    artifacts={}
    for name,source in files.items():
        shutil.copyfile(source,root/source.name);artifacts[name]=source.name
    return write_manifest(root,stage='minimized',artifacts=artifacts,units='nm',
        profile_id=receipt['id'],organism=receipt.get('organism'),strain=receipt.get('strain'),
        membrane_geometry={'kind':'plane','periodic':True},sources=receipt['sources'],
        assumptions=receipt['assumptions'],charge_e=receipt['charge_e'],
        minimization=read_minimization(run),geometry_checks={'passed':True,'original':checks},
        limits=['Restrained minimum only; not equilibrated or biologically validated.'])
