#!/usr/bin/env python3
"""Verify complete, passing, identity-bound case evidence across environments.

Exit 0: every expected blocking scope passed and reproduced; 1: observed failure
or disagreement; 2: missing, malformed, incompatible, or insufficient evidence.
Both nonzero outcomes block release. Machine agreement never grants independent
scientific approval. Legacy summaries remain readable but cannot pass this gate.
"""
from __future__ import annotations
import argparse
import json
import sys
import tarfile
from pathlib import Path, PurePosixPath
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
QUALIFICATION = ROOT / 'yauvi-structural-workbench' / 'benchmarks' / 'qualification-v2'
sys.path.insert(0, str(QUALIFICATION))
from evidence_contract import execution_identity
MAX_SUMMARY_BYTES = 16 * 1024 * 1024
OUTCOME_FIELDS = ('stratum_state', 'cases_passed', 'cases_failed', 'cases_total', 'controls_passed', 'controls_total')

def environment_key(runtime):
    return f"{runtime.get('platform', '?')}/{runtime.get('machine', '?')}"

def load_summaries(evidence_dir):
    found = []
    for path in sorted(evidence_dir.rglob('EXECUTION_SUMMARY.json')):
        source = path.relative_to(evidence_dir).as_posix()
        try:
            if path.is_symlink() or path.stat().st_size > MAX_SUMMARY_BYTES:
                raise ValueError('unsafe or oversized summary')
            found.append({'source': source, 'summary': json.loads(path.read_bytes())})
        except (OSError, ValueError) as exc:
            found.append({'source': source, 'error': str(exc)})
    for archive in sorted(evidence_dir.rglob('*.tar.gz')):
        source = archive.relative_to(evidence_dir).as_posix()
        try:
            with tarfile.open(archive, 'r:gz') as tar:
                for member in tar:
                    if not member.name.endswith('EXECUTION_SUMMARY.json'):
                        continue
                    name = PurePosixPath(member.name)
                    if not member.isfile() or name.is_absolute() or '..' in name.parts or member.size > MAX_SUMMARY_BYTES:
                        raise ValueError('unsafe or oversized summary member')
                    # Read bytes only: never extract untrusted archives, on any Python version.
                    handle = tar.extractfile(member)
                    if handle is None:
                        raise ValueError('summary has no file content')
                    found.append({'source': source + ':' + member.name,
                                  'summary': json.loads(handle.read(MAX_SUMMARY_BYTES + 1))})
        except (OSError, ValueError, tarfile.TarError) as exc:
            found.append({'source': source, 'error': str(exc)})
    return found

def _records(panel, definition):
    """Validate counts against exact IDs and conjunctions of required checks."""
    verdicts = panel['record_verdicts']
    normalized = {}
    for kind in ('cases', 'controls'):
        rows = verdicts[kind]
        ids = [r['record_id'] for r in rows]
        expected = [r['record_id'] for r in definition.get('records' if kind == 'cases' else 'controls', [])]
        if len(ids) != len(set(ids)) or set(ids) != set(expected):
            raise ValueError(f'{kind}: duplicate, missing, or unexpected record IDs')
        for row in rows:
            checks = row['checks']
            names = [c['check'] for c in checks]
            if not checks or len(names) != len(set(names)):
                raise ValueError('missing or duplicate checks')
            if any(c.get('required') is not True or type(c.get('passed')) is not bool for c in checks):
                raise ValueError('invalid required-check verdict')
            if type(row.get('passed')) is not bool or row['passed'] != all(c['passed'] for c in checks):
                raise ValueError('record verdict contradicts required checks')
        counts = panel[kind]
        passed = sum(r['passed'] for r in rows)
        if any(type(counts.get(k)) is not int for k in ('total', 'passed')):
            raise ValueError('counts must be integers')
        if counts.get('total') != len(rows) or counts.get('passed') != passed:
            raise ValueError('summary count contradicts records')
        if kind == 'cases' and counts.get('failed') != len(rows) - passed:
            raise ValueError('failed count contradicts records')
        normalized[kind] = sorted((r['record_id'], r['passed'], sorted((c['check'], c['passed']) for c in r['checks'])) for r in rows)
    return normalized

def build_report(loaded, min_environments=2, *, manifest=None, identity=None):
    if min_environments < 2:
        raise ValueError('at least two environments are required')
    manifest = manifest if manifest is not None else json.loads((QUALIFICATION / 'PANEL_MANIFEST.json').read_text())
    identity = identity if identity is not None else execution_identity()
    definitions = {p['workflow']: p for p in manifest['panels']}
    expected = {s.split(':', 1)[0] for s in manifest['release_blocking_scopes']}
    runners, errors, seen = [], [], set()
    for entry in loaded:
        source = entry['source']
        try:
            if entry.get('error'):
                raise ValueError(entry['error'])
            summary = entry['summary']
            runtimes = summary['recorded_on']
            if len(runtimes) != 1 or any(not runtimes[0].get(k) for k in ('platform', 'machine', 'python')):
                raise ValueError('a single fully identified runtime is required')
            runtime = runtimes[0]
            env = environment_key(runtime)
            label = env + ' py' + runtime['python']
            if label in seen:
                raise ValueError('duplicate runner runtime')
            seen.add(label)
            panels = {}
            for panel in summary['panels']:
                workflow = panel['workflow']
                if workflow not in definitions or workflow in panels:
                    raise ValueError('unknown or duplicate panel')
                if panel.get('evidence_identity') != identity:
                    raise ValueError('missing or incompatible code, manifest, source lock, or protocol identity')
                expected_digest = identity.get('execution_panels', {}).get(workflow)
                if not expected_digest or panel.get('execution_panel_sha256') != expected_digest:
                    raise ValueError('missing or incompatible execution-panel digest')
                records = _records(panel, definitions[workflow])
                cases, controls = panel['cases'], panel['controls']
                passed = (panel['stratum_state'] == 'passed' and cases['total'] > 0
                          and cases['failed'] == 0 and cases['passed'] == cases['total']
                          and controls['passed'] == controls['total'] and not panel.get('coverage', {}).get('unmet'))
                panels[workflow] = {'passed': passed, 'panel_sha256': panel['execution_panel_sha256'],
                                    'signature': json.dumps({'state': panel['stratum_state'], 'records': records,
                                                            'coverage': panel.get('coverage', {})}, sort_keys=True)}
            runners.append({'source': source, 'environment': env, 'label': label, 'panels': panels})
        except (KeyError, TypeError, ValueError) as exc:
            errors.append({'source': source, 'why': str(exc)})
    panels = []
    for workflow in sorted(expected | {w for r in runners for w in r['panels']}):
        observations = [(r, r['panels'][workflow]) for r in runners if workflow in r['panels']]
        environments = sorted({r['environment'] for r, _ in observations})
        signatures = {p['signature'] for _, p in observations}
        digests = {p['panel_sha256'] for _, p in observations}
        agreed = len(signatures) == 1
        passed = bool(observations) and all(p['passed'] for _, p in observations)
        internally = sorted(e for e in environments if len({p['signature'] for r,p in observations if r['environment'] == e}) > 1)
        definition = definitions[workflow]
        requirements = [r for r in definition['requirements'] if f"{workflow}:{r['stratum']}" in manifest['release_blocking_scopes']]
        coverage_complete = bool(requirements) and all(
            r['count'] > 0 and sum(row['stratum'] == r['stratum'] and row.get('split') == r.get('split')
                                  for row in definition.get('records', [])) == r['count']
            for r in requirements)
        complete = bool(observations) and len(observations) == len(runners) and (workflow not in expected or coverage_complete)
        panels.append({'workflow': workflow, 'release_blocking': workflow in expected,
                       'environment_count': len(environments), 'environments': environments,
                       'agreed_across_environments': agreed, 'passed_on_every_environment': passed,
                       'internally_inconsistent_environments': internally,
                       'complete': complete, 'protocols_match': len(digests) == 1,
                       'reproduced': bool(complete and agreed and passed and len(digests) == 1
                                          and len(environments) >= min_environments and not errors)})
    blocking = [p for p in panels if p['release_blocking']]
    failed = [p['workflow'] for p in blocking if p['environment_count'] and not p['passed_on_every_environment']]
    disagreed = [p['workflow'] for p in blocking if p['environment_count'] and not p['agreed_across_environments']]
    incomplete = [p['workflow'] for p in blocking if not p['complete'] or p['environment_count'] < min_environments or not p['protocols_match']]
    success = bool(blocking) and all(p['reproduced'] for p in blocking) and not errors
    return {'schema_version': '2.0', 'min_environments': min_environments,
            'evidence_identity': identity, 'runners_read': len(loaded), 'runners_usable': len(runners),
            'environments_observed': sorted({r['environment'] for r in runners}),
            'environment_count': len({r['environment'] for r in runners}), 'panels': panels,
            'unusable_runners': errors, 'release_blocking_panels_failed': failed,
            'release_blocking_panels_disagreed': disagreed,
            'release_blocking_panels_under_minimum': [p['workflow'] for p in blocking if p['environment_count'] < min_environments],
            'release_blocking_panels_incomplete': incomplete,
            'release_blocking_panels_reproduced': [p['workflow'] for p in blocking if p['reproduced']],
            'every_release_blocking_panel_reproduced': success,
            'exit_code': 2 if errors or incomplete or not blocking else 1 if failed or disagreed else 0 if success else 2,
            'independence_note': 'Environment agreement does not establish independent scientific approval.'}

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--evidence-dir', type=Path, required=True)
    parser.add_argument('--json-out', type=Path, default=QUALIFICATION / 'results/CROSS_MACHINE_REPRODUCTION.json')
    parser.add_argument('--min-environments', type=int, default=2)
    args = parser.parse_args(argv)
    if not args.evidence_dir.is_dir() or args.min_environments < 2:
        print('Missing evidence directory or fewer than two requested environments.', file=sys.stderr)
        return 2
    report = build_report(load_summaries(args.evidence_dir), args.min_environments)
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: report[k] for k in ('exit_code', 'every_release_blocking_panel_reproduced',
          'release_blocking_panels_failed', 'release_blocking_panels_disagreed', 'release_blocking_panels_incomplete', 'unusable_runners')}, indent=2))
    return report['exit_code']

if __name__ == '__main__':
    raise SystemExit(main())
