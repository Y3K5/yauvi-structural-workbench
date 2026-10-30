"""Read-only capability cards over existing definitions, sources and saved checks.

No executable probe, provider import, input validation or job is started here.
File presence is deliberately not reported as runtime or scientific readiness.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import stat
import sys

from yauvi_platform.structural_workbench import store as store_module
from yauvi_platform.structural_workbench.local_files import LocalSourceUnavailable

METHODS = {
    'structure_qc': ('Structure QC', 'StructQC and existing Workbench custody',
                     ['Coordinate/provenance quality report', 'Exact input and method records']),
    'subtractive_collection_screen': ('Collection research dossier', 'Existing local collection provider',
                                     ['Per-protein dossiers', 'All-input accounting, comparisons and holds']),
    'antigenbench_synthetic_check': ('AntigenBench synthetic check', 'Separate AntigenBench, through the existing controlled example',
                                   ['Synthetic software contract report', 'Exact fixture and source identities']),
}
SHA = re.compile(r'^[0-9a-f]{64}$')
MAX_BYTES = 2 * 1024 * 1024


def _bytes(path):
    """Bounded, no-follow read, refusing offloaded or changing files."""
    path = Path(path).absolute()
    if path.resolve() != path:
        raise ValueError('Source path uses a shortcut')
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, 'rb') as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode) or before.st_size > MAX_BYTES:
            raise ValueError('Source is not a bounded regular file')
        if sys.platform == 'darwin' and getattr(before, 'st_flags', 0) & 0x40000000:
            raise LocalSourceUnavailable('Source is not stored locally')
        raw = stream.read(MAX_BYTES + 1)
        after = os.fstat(stream.fileno())
        fields = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns')
        if len(raw) > MAX_BYTES or any(getattr(before, f) != getattr(after, f) for f in fields):
            raise ValueError('Source changed while reading')
        return raw


def _json(path):
    value = json.loads(_bytes(path))
    if not isinstance(value, dict):
        raise ValueError('Expected a saved object record')
    return value


def _digest(value):
    return hashlib.sha256(store_module._canonical(value)).hexdigest()


class CapabilityView:
    def __init__(self, store, *, reference_fetch_enabled=False):
        self.store = store
        self.reference_fetch_enabled = bool(reference_fetch_enabled)

    def _sources(self, method):
        # The existing execution configuration is the authority; never search for a
        # substitute checkout or run a provider to make a card look available.
        adapter_root = Path(store_module.__file__).resolve().parent
        files = {'workbench-adapter': adapter_root / {
            'structure_qc': 'store.py', 'subtractive_collection_screen': 'collection_screen.py',
            'antigenbench_synthetic_check': 'module_examples.py'}[method]}
        configuration = {'method':method, 'python':sys.executable}
        state = 'configured'
        runtime = 'not_checked'
        if method == 'structure_qc':
            if not self.store.installed_mode:
                root = self.store.source_root / 'structqc' / 'src' / 'structqc'
            else:
                spec = importlib.util.find_spec('structqc')
                root = Path(spec.origin).parent if spec and spec.origin else None
            if root:
                files.update({'structqc-core': root/'core.py', 'structqc-cli': root/'cli.py'})
                configuration['source_root'] = str(root)
            else:
                state = 'source_missing'
            if not Path(sys.executable).is_file():
                runtime, state = 'missing_interpreter', 'runtime_missing'
        elif method == 'subtractive_collection_screen':
            root_value = os.environ.get('YAUVI_ATLAS_ROOT')
            python_value = os.environ.get('YAUVI_ATLAS_PYTHON')
            configuration.update(root=root_value, python=python_value,
                                 provider=os.environ.get('YAUVI_COLLECTION_PROVIDER'))
            if not root_value or not python_value:
                state = 'not_configured'
            else:
                root = Path(root_value).expanduser().absolute()
                python = Path(python_value).expanduser().absolute()
                if not python.is_file():
                    runtime, state = 'missing_interpreter', 'runtime_missing'
                # Reuse the executor's allowlist and fallback order; this imports
                # only the already registered adapter, never the private provider.
                from yauvi_platform.structural_workbench.collection_screen import provider_module
                try:
                    module = provider_module(root)
                    files.update({'collection-worker': root/module/'worker.py',
                                  'sequence-worker': root/'sequence_screen_v1/worker.py',
                                  'atlas-worker': root/'atlas_bridge_v1/worker.py'})
                except ValueError:
                    if state == 'configured':
                        state = 'source_missing'
        else:
            value = os.environ.get('YAUVI_ANTIGENBENCH_ROOT')
            configuration['root'] = value
            if not value:
                state = 'not_configured'
            else:
                root = Path(value).expanduser().absolute()
                files.update({'antigenbench-package': root/'antigenbench/__init__.py',
                              'synthetic-manifest': root/'corpus/synthetic/manifest.yaml'})
            if not Path(sys.executable).is_file():
                runtime, state = 'missing_interpreter', 'runtime_missing'
        pins = {}
        for name, path in files.items():
            try:
                pins[name] = {'state': 'present', 'sha256': hashlib.sha256(_bytes(path)).hexdigest()}
            except LocalSourceUnavailable:
                pins[name] = {'state': 'offloaded'}
            except (OSError, ValueError, RuntimeError):
                pins[name] = {'state': 'unavailable'}
        if state == 'configured' and any(p['state'] != 'present' for p in pins.values()):
            state = 'source_unavailable'
        return {'configuration': state, 'configuration_sha256':_digest(configuration),
                'runtime': runtime, 'files': pins,
                'scope': 'Named adapter/source files only. Full dependencies, references and executable behavior are not checked by inspection.'}

    def _saved_check(self, manifest):
        record_path = self.store._case_dir(manifest['analysis_id'])/'PREFLIGHT.json'
        if not record_path.exists():
            return {'state': 'not_recorded', 'label': 'No saved input check for this case'}
        try:
            raw = _bytes(record_path)
            record = json.loads(raw)
            if not isinstance(record, dict):
                raise ValueError('Expected a saved check object')
            revision = record.get('revision_sha256')
            if (record.get('contract_id') != 'analysis_preflight'
                    or record.get('analysis_id') != manifest['analysis_id']
                    or record.get('analysis_type') != manifest['analysis_type']
                    or not isinstance(revision, str) or not SHA.fullmatch(revision)
                    or type(record.get('valid')) is not bool):
                raise ValueError('Saved check attribution is invalid')
            old = _json(self.store._case_dir(manifest['analysis_id'])/'revisions'/(revision+'.json'))
            if (_digest({k:v for k,v in old.items() if k != 'revision_sha256'}) != revision
                    or old.get('analysis_id') != manifest['analysis_id']):
                raise ValueError('Saved revision identity differs')
            fields = ('analysis_type','inputs','parameters','question','subject_id','source_adoptions')
            current = all(old.get(k) == manifest.get(k) for k in fields)
            state = ('recorded_pass' if record['valid'] else 'recorded_blocked') if current else 'inputs_changed'
            label = {'recorded_pass':'A pass was recorded for these declarations',
                     'recorded_blocked':'A refusal was recorded for these declarations',
                     'inputs_changed':'Inputs or declarations changed after the saved check'}[state]
            return {'state': state, 'label': label, 'revision_sha256': revision,
                    'record_sha256': hashlib.sha256(raw).hexdigest(),
                    'limit':'Historical record only. Input bytes, method sources and runtime have not been rechecked now.'}
        except (OSError, ValueError, TypeError, KeyError, RuntimeError):
            return {'state':'unavailable', 'label':'Saved check cannot be attributed safely; use the existing input check'}

    def _inputs(self, method, definition, analysis_id):
        if analysis_id is None:
            return {'state':'not_selected', 'label':'Choose a saved input case to inspect its recorded check'}
        if not isinstance(analysis_id, str):
            raise ValueError('Choose one saved input case')
        manifest = _json(self.store._case_dir(analysis_id)/'ANALYSIS_CASE.json')
        if manifest.get('analysis_id') != analysis_id or manifest.get('analysis_type') != method:
            raise ValueError('Selected case does not belong to this method')
        if _digest({k:v for k,v in manifest.items() if k != 'revision_sha256'}) != manifest.get('revision_sha256'):
            raise ValueError('Selected case identity changed')
        counts = {}
        for item in manifest.get('inputs', []):
            counts[item['role']] = counts.get(item['role'], 0) + 1
        missing = [r['label'] for r in definition['inputs']
                   if counts.get(r['role'], 0) < int(r.get('minimum_files', 1 if r['required'] else 0))]
        return {'state':'missing_inputs' if missing else 'present_unchecked',
                'label':'Required inputs missing' if missing else 'Input entries present; contents not rechecked',
                'analysis_id': analysis_id, 'revision_sha256':manifest['revision_sha256'],
                'missing':missing, 'recorded_check':self._saved_check(manifest),
                'open_inputs_url':'/?case='+analysis_id}

    def detail(self, method, *, expected_sha256=None, analysis_id=None):
        if method not in METHODS:
            raise KeyError('Unknown capability')
        if expected_sha256 is not None and (not isinstance(expected_sha256,str) or not SHA.fullmatch(expected_sha256)):
            raise ValueError('Refresh the capability view before inspecting a method')
        definition = self.store._definitions.get(method)
        title, owner, outputs = METHODS[method]
        if definition is None:
            return {'id':method,'title':title,'registration':'unavailable','source_state':'not_registered',
                    'next_action':'Refresh after checking the current Workbench build','scientific_state':'not_assessed'}
        source = self._sources(method)
        identity = _digest({'definition':definition,'source':source})
        changed = expected_sha256 is not None and identity != expected_sha256
        inputs = self._inputs(method, definition, analysis_id)
        return {'id':method,'title':title,'owner':owner,'registration':'registered',
                'question':definition['question'],'source':source,'identity_sha256':identity,
                'source_state':'changed_since_view' if changed else source['configuration'],
                'synthetic_only':method=='antigenbench_synthetic_check',
                'inputs':[{'role':r['role'],'label':r['label'],'required':r['required'],
                           'formats':r.get('extensions',[])} for r in definition['inputs']],
                'outputs':definition.get('receives') or outputs,'input_state':inputs,
                'scientific_state':'synthetic_software_only' if method=='antigenbench_synthetic_check' else 'not_qualified_by_this_view',
                'claim_ceiling':definition['claim_ceiling'],
                'next_action':'Refresh source status before using the existing input check' if changed
                    else 'Use the existing Check inputs action before deliberately executing a method',
                'inspection_only':True}

    def catalogue(self):
        return {'schema_version':'1.0','capabilities':[self.detail(key) for key in METHODS],
                'privacy':{'inspection':'Local metadata only; no tool probe, job or public retrieval',
                           'public_reference_fetch_enabled':self.reference_fetch_enabled,
                           'network_qualification':'Tool/subprocess network behavior is not established by this view'},
                'scope':'Source/configuration, input records, execution and scientific evidence are separate states.'}
