"""Local, append-only review receipts and allowlisted region replay bundles.

Bundles contain data and frozen method copies, never commands to execute. Replay
uses installed, identity-checked methods; supplied code is never imported.
"""
from __future__ import annotations
import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import shutil
from pathlib import Path

from .biological_case import BiologicalCase
from .regions import canonical, digest, view_binding

ADJUDICATIONS = {'supported_discrepancy', 'reference_applicability_limitation', 'unresolved'}
METHOD_FILES = ('biological_case.py', 'regions.py', 'evidence_tools.py')
PACKAGES = ('gemmi', 'numpy', 'scipy', 'biopython')


def runtime():
    return {'python': '.'.join(platform.python_version().split('.')[:2]),
            'system': platform.system(), 'machine': platform.machine(),
            'dependencies': {name: importlib.metadata.version(name) for name in PACKAGES}}


def method_locks():
    return {name: hashlib.sha256(path.read_bytes()).hexdigest() for name,path in method_paths().items()}


def method_paths():
    paths = {name:Path(__file__).parent/name for name in METHOD_FILES}
    # Assembly rotation checks and optional engine/placement adapters are part
    # of replay, so lock their executable Python sources as well.
    for package in ('structqc','memorient','bio_orient'):
        spec=importlib.util.find_spec(package)
        if spec is None: raise ValueError('Required replay method package absent: '+package)
        root=Path(next(iter(spec.submodule_search_locations)))
        paths.update({package+'/'+p.relative_to(root).as_posix():p for p in sorted(root.rglob('*.py')) if '__pycache__' not in p.parts})
    return paths


def safe_file(root, relative):
    root = Path(root).resolve()
    p = Path(relative)
    if p.is_absolute() or '..' in p.parts or not p.parts:
        raise ValueError('Bundle paths must be relative and stay within the bundle')
    path = root/p
    if any((root/Path(*p.parts[:i])).is_symlink() for i in range(1, len(p.parts)+1)):
        raise ValueError('Symlink bundle inputs are unsupported')
    if not path.is_file(): raise ValueError('Required bundle input is missing: '+relative)
    return path


def selectors(case):
    return [{'structure_id':s['id'], 'model_id':s['models'][0], 'assembly_id':'asu', 'conformer':'auto'}
            for s in case.summary()['structures']]


def selected_view(case, selector):
    if set(selector) != {'structure_id','model_id','assembly_id','conformer'}:
        raise ValueError('Complete explicit view selection is required')
    return case.view(selector['structure_id'], selector['model_id'], selector['assembly_id'], selector['conformer'])


def review(case, selector, record, destination):
    """Retain a user adjudication beside, without altering, the original finding."""
    view = selected_view(case, selector)
    if record.get('binding') != view_binding(view): raise ValueError('Review view identity mismatch')
    if record.get('adjudication') not in ADJUDICATIONS: raise ValueError('Unsupported review classification')
    for key in ('reviewer','rationale','target_id','target_sha256'):
        if not isinstance(record.get(key),str) or not record[key].strip(): raise ValueError('Review requires '+key)
    if not record.get('source_ids') or any(s not in case.sources for s in record['source_ids']):
        raise ValueError('Review requires existing checksum-locked source evidence')
    targets = {r['id']:r for r in view['regions']}
    # Original engine observations can be reviewed through their bound evidence record.
    targets.update({r['id']:r for r in view['evidence']})
    target = targets.get(record['target_id'])
    original_parent_sha256 = digest(target)
    pointer = record.get('target_pointer','')
    if pointer:
        if not isinstance(pointer,str) or not pointer.startswith('/'): raise ValueError('Invalid finding JSON pointer')
        try:
            for part in pointer[1:].split('/'):
                key=part.replace('~1','/').replace('~0','~')
                target=target[int(key)] if isinstance(target,list) and key.isdecimal() else target[key]
        except (TypeError,KeyError,IndexError,ValueError) as exc: raise ValueError('Review finding pointer is unavailable') from exc
    if target is None or digest(target) != record['target_sha256']:
        raise ValueError('Review target identity mismatch')
    receipt = {'schema_version':'1.0','case_manifest_sha256':hashlib.sha256((case.directory/'case.json').read_bytes()).hexdigest(),
        'selection':selector,'record':record,'original_target':target,'original_parent_sha256':original_parent_sha256,
        'sources':[{k:v for k,v in case.sources[s].items() if k != 'path'} for s in record['source_ids']],
        'claim_limit':'Recorded human review; original findings remain unchanged. This is not independent chemical or biological validation.'}
    raw = canonical(receipt); receipt_id = hashlib.sha256(raw).hexdigest()
    destination = Path(destination)
    # Never write reviews into an immutable case/source directory.
    if destination.resolve() == case.directory or case.directory in destination.resolve().parents:
        raise ValueError('Review receipts must be stored outside the original case')
    destination.mkdir(parents=True,exist_ok=True)
    path = destination/(receipt_id+'.json')
    if path.exists():
        if path.read_bytes() != raw: raise ValueError('Existing receipt differs')
    else:
        with path.open('xb') as f: f.write(raw)
    return {'receipt_id':receipt_id,'file':str(path),'adjudication':record['adjudication']}


def bundle(case, destination, selections=None, *, atol=1e-6):
    if not isinstance(atol,(int,float)) or not 0 < atol <= 1e-4: raise ValueError('Comparison tolerance must be in (0, 0.0001]')
    destination = Path(destination).resolve()
    if destination.exists() and any(destination.iterdir()): raise ValueError('Bundle destination must be empty')
    if destination == case.directory or case.directory in destination.parents: raise ValueError('Bundle must be outside original case')
    selections = selections or selectors(case)
    expected = [selected_view(case,s) for s in selections]
    destination.mkdir(parents=True,exist_ok=True)
    inputs = destination/'case'; inputs.mkdir()
    (inputs/'case.json').write_bytes((case.directory/'case.json').read_bytes())
    for source in case.sources.values():
        target = inputs/source['path']; target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(case.payloads[source['id']])
    (destination/'methods').mkdir()
    locks = method_locks()
    for name,path in method_paths().items():
        target=destination/'methods'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
    (destination/'EXPECTED_VIEWS.json').write_bytes(canonical(expected))
    (destination/'README.txt').write_text(
        'Offline region replay. Requires this exact installed Workbench and the locked runtime below.\n'
        'Run: yauvi evidence replay --bundle /path/to/extracted/bundle\n'
        'Region mapping, component contacts and selections are recalculated. Other engine evidence is checksum-verified only.\n'
        'Frozen method files are supplied for inspection; replay never executes bundled code.\n'
        'Same-machine replay does not establish independent reproduction.\n')
    files = {p.relative_to(destination).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted(destination.rglob('*')) if p.is_file()}
    manifest = {'schema_version':'1.0','recipe':'coordinate_bound_regions_v1','selections':selections,
        'runtime':runtime(),'method_locks':locks,'files':files,'comparison':{'absolute_tolerance':atol,'relative_tolerance':0},
        'recalculation':['reference annotation mapping','selected conformer and chain-copy atom identities','5 Å heavy-atom contacts'],
        'checksum_verification_only':[e['id'] for e in case.document.get('evidence',[])],
        'scope_limit':'Static observations only. Independent researcher and second-machine reproduction are not established.'}
    (destination/'REPLAY_MANIFEST.json').write_bytes(canonical(manifest))
    return manifest


def differences(expected, actual, atol, path='$'):
    if type(expected) is not type(actual): return [path+': type differs']
    if isinstance(expected,dict):
        if expected.keys() != actual.keys(): return [path+': fields differ']
        return [d for k in expected for d in differences(expected[k],actual[k],atol,path+'.'+k)]
    if isinstance(expected,list):
        if len(expected) != len(actual): return [path+': length differs']
        return [d for i,(a,b) in enumerate(zip(expected,actual)) for d in differences(a,b,atol,path+'['+str(i)+']')]
    if isinstance(expected,float):
        import math
        return [] if math.isfinite(expected) and math.isfinite(actual) and abs(expected-actual) <= atol else [path+': numeric value differs']
    return [] if expected == actual else [path+': value differs']


def replay(directory):
    directory = Path(directory).resolve()
    manifest = json.loads(safe_file(directory,'REPLAY_MANIFEST.json').read_bytes())
    if manifest.get('schema_version') != '1.0' or manifest.get('recipe') != 'coordinate_bound_regions_v1':
        raise ValueError('Unsupported replay recipe')
    if manifest.get('method_locks') != method_locks(): raise ValueError('Installed method identity differs from the frozen bundle')
    try: current = runtime()
    except importlib.metadata.PackageNotFoundError as exc: raise ValueError('Required replay dependency is absent: '+str(exc)) from exc
    if manifest.get('runtime') != current: raise ValueError('Required replay runtime/dependency identity differs')
    for relative,sha in manifest['files'].items():
        if hashlib.sha256(safe_file(directory,relative).read_bytes()).hexdigest() != sha:
            raise ValueError('Bundle checksum differs: '+relative)
    expected = json.loads(safe_file(directory,'EXPECTED_VIEWS.json').read_bytes())
    case = BiologicalCase(directory/'case')
    actual = [selected_view(case,s) for s in manifest['selections']]
    atol = manifest['comparison']['absolute_tolerance']
    if manifest['comparison'].get('relative_tolerance') != 0 or not 0 < atol <= 1e-4: raise ValueError('Unsupported comparison tolerance')
    changed = differences(expected,actual,atol)
    return {'schema_version':'1.0','state':'passed' if not changed else 'different',
        'replay_manifest_sha256':hashlib.sha256((directory/'REPLAY_MANIFEST.json').read_bytes()).hexdigest(),
        'runtime':current,'method_locks':method_locks(),'view_count':len(actual),'differences':changed[:100],
        'recalculated':manifest['recalculation'],'checksum_verified_only':manifest['checksum_verification_only'],
        'comparison':manifest['comparison'],'independent_reproduction':'not_established'}
