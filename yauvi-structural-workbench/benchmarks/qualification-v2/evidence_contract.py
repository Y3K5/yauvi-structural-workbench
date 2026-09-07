"""Portable identities for newly executed qualification evidence, never retrofitted."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()

def file_digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def execution_identity():
    files = {}
    for directory in ('structural-workbench', 'platform/src', 'sources/src', 'structqc/src',
                      'Membrane Orientor/memorient/src', 'state-atlas/src', 'site-context/src',
                      'activity-state/src', 'assembly-context/src', 'sf-csa/src', 'tools'):
        for path in sorted((ROOT / directory).rglob('*.py')):
            if not {'__pycache__', 'build', '.venv'}.intersection(path.parts):
                files[path.relative_to(ROOT).as_posix()] = file_digest(path)
    for path in sorted(HERE.glob('*.py')):
        files[path.relative_to(ROOT).as_posix()] = file_digest(path)
    files['pyproject.toml'] = file_digest(ROOT / 'pyproject.toml')
    return {'code_sha256': digest(files), 'manifest_sha256': file_digest(HERE / 'PANEL_MANIFEST.json'),
            'source_lock_sha256': file_digest(HERE / 'SOURCE_LOCK.json'),
            'protocols_sha256': digest({p.name: file_digest(p) for p in sorted(HERE.glob('ADOPTION_DRAFT*.json'))}),
            'execution_panels': {json.loads(p.read_text()).get('workflow', 'structure_qc'): file_digest(p)
                                 for p in sorted(HERE.glob('ADOPTION_DRAFT*.json'))}}

def record_verdicts(status):
    """Keep exact IDs and required check verdicts; omit paths and numeric drift."""
    return {kind: [{ 'record_id': r['record_id'], 'passed': r['passed'],
                    'checks': [{k: c[k] for k in ('check', 'required', 'passed')} for c in r['checks'] if c.get('required', True)]}
                   for r in status.get(kind, [])] for kind in ('cases', 'controls')}
