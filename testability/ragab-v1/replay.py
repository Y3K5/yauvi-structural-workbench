"""Standalone RagAB replay. No download, installation or upload during a run."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import importlib.metadata as metadata
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
ENGINE_NAMES = ('structqc', 'state-atlas', 'assembly-context')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def verify_payload():
    records = read(ROOT / 'PAYLOAD.json')['files']
    for relative, expected in records.items():
        if sha(ROOT / relative) != expected:
            raise RuntimeError('Kit member changed: ' + relative)
    return records


def verify_sources(source_root):
    records = read(ROOT / 'SOURCES.json')['sources']
    checked = {}
    for record in records:
        items = [(record['artifact'], record['sha256'])]
        if record['validation']['identity_revision_match']:
            items.append((record['validation']['xml_artifact'], record['validation']['xml_sha256']))
        for relative, expected in items:
            # The caller supplies a directory containing the sources/ tree.
            if sha(source_root / relative) != expected:
                raise RuntimeError('Source missing or checksum mismatch: ' + relative)
            checked[relative] = expected
    return checked


def runtime(freesasa):
    profile = read(ROOT / 'RUNTIME_PROFILE.json')
    if f'{sys.version_info.major}.{sys.version_info.minor}' != profile['python_major_minor']:
        raise RuntimeError('Python differs from the declared runtime profile')
    versions = {name: metadata.version(name) for name in profile['dependencies']}
    if versions != profile['dependencies']:
        raise RuntimeError('Dependencies differ from the declared runtime profile')
    binary = Path(freesasa).resolve(strict=True)
    check = subprocess.run([str(binary), '--version'], capture_output=True, text=True, timeout=30, check=True)
    if not check.stdout.startswith('FreeSASA ' + profile['freesasa_version'] + '\n'):
        raise RuntimeError('Unsupported FreeSASA version')
    return {
        'python_version': sys.version, 'python_executable': str(Path(sys.executable).resolve()),
        'python_executable_sha256': sha(Path(sys.executable).resolve()),
        'platform': platform.platform(), 'machine': platform.machine(), 'packages': versions,
        'freesasa_executable': str(binary), 'freesasa_sha256': sha(binary),
        'freesasa_version_output': check.stdout,
        'exact_original_runtime_replay': False,
    }


def checked_worker(argv, cwd, log, required, env=None):
    """Nonzero exits, missing outputs and nonfinite results are never passes."""
    with Path(log).open('x') as stream:
        result = subprocess.run(argv, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT,
                                env=env, timeout=300)
    if result.returncode != 0:
        raise RuntimeError('Worker failed; inspect ' + str(Path(log).name))
    for path in required:
        if not Path(path).is_file():
            raise RuntimeError('Required output missing: ' + str(Path(path).name))
        if Path(path).suffix == '.json':
            json.loads(Path(path).read_text(), parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError('Nonfinite JSON value: ' + value)))


def differences(expected, actual, policy, pointer=''):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(expected) != set(actual):
            return [pointer + ': field membership differs']
        result = []
        for key in expected:
            child = (pointer + '.' + key).lstrip('.')
            if child not in policy['omit_raw_fields']:
                result.extend(differences(expected[key], actual[key], policy, child))
        return result
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            return [pointer + ': list membership differs']
        return [item for i, value in enumerate(expected)
                for item in differences(value, actual[i], policy, pointer + '[' + str(i) + ']')]
    if isinstance(expected, float):
        tolerance = policy['surface_area_absolute_tolerance_A2'] if 'A2' in pointer else policy['other_float_absolute_tolerance']
        if isinstance(actual, bool) or not isinstance(actual, (int, float)) or not math.isfinite(actual):
            return [pointer + ': nonfinite or nonnumeric result']
        return [] if math.isclose(expected, actual, rel_tol=0, abs_tol=tolerance) else [pointer + ': numerical difference']
    return [] if type(expected) is type(actual) and expected == actual else [pointer + ': exact field differs']


def compare_run(destination):
    expected = read(ROOT / 'EXPECTED.json')
    measurements = destination / 'runs' / 'replay-001'
    comparisons = []
    for query, record in expected['state'].items():
        relative = query + '/state/summary.json'
        comparisons.append({'path': relative, 'differences': differences(record, read(measurements / relative), expected['policy'])})
    for key, record in expected['burial'].items():
        relative = key + '/summary.json'
        comparisons.append({'path': relative, 'differences': differences(record, read(measurements / relative), expected['policy'])})
    for relative, record in expected['raw_assembly'].items():
        comparisons.append({'path': relative, 'differences': differences(record, read(measurements / relative), expected['policy'])})
    return {'comparison_groups': len(comparisons), 'passed': sum(not row['differences'] for row in comparisons),
            'comparisons': comparisons, 'independent_operator_verified': False, 'independent_second_host_verified': False,
            'biological_validation': False, 'scope': 'Numerical replay agreement under the declared profile only'}


def execute(destination, source_root, freesasa):
    if destination.exists() or destination.is_symlink():
        raise RuntimeError('Destination already exists; previous attempts are preserved')
    payload = verify_payload()
    source_files = verify_sources(source_root)
    observed = runtime(freesasa)
    # Runtime and exact kit/source identities are recorded before calculation.
    destination.mkdir(parents=True, exist_ok=False)
    write(destination / 'REGISTRATION.json', {'utc': utc(), 'payload_sha256': sha(ROOT / 'PAYLOAD.json'),
        'payload_files': payload, 'sources': source_files, 'runtime': observed,
        'policy': read(ROOT / 'EXPECTED.json')['policy'],
        'scope': 'Reprepare fixed coordinates, controls, ten states and twelve burial legs; primary stopping rule retained'})
    try:
        for relative in source_files:
            target = destination / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_root / relative, target)
        for name in ('CONFIG.json', 'SOURCES.json'):
            shutil.copyfile(ROOT / name, destination / name)
        shutil.copytree(ROOT / 'frozen-adapter', destination / 'scripts')
        shutil.copytree(ROOT / 'engines', destination / 'engines', ignore=shutil.ignore_patterns('*.egg-info', '__pycache__'))
        write(destination / 'environment/ENVIRONMENT.json', {'freesasa': {'sha256': observed['freesasa_sha256']}})
        locked = {str(path.relative_to(destination)): sha(path) for path in destination.rglob('*') if path.is_file()}
        write(destination / 'PROTOCOL_LOCK.json', {'created_utc': utc(), 'files': locked,
            'note': 'New replay lock; original protocol identity is preserved in ORIGIN.json'})
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        env['PYTHONPATH'] = os.pathsep.join(str(destination / 'engines' / name / 'src') for name in ENGINE_NAMES)
        env['PATH'] = str(Path(freesasa).resolve().parent) + os.pathsep + os.defpath
        study = destination / 'scripts/study.py'
        for action, required in (
            ('prepare', [destination / 'prepared/PREPARATION.json']),
            ('controls', [destination / 'controls/RESULT.json']),
            ('freeze', [destination / 'EXECUTION_LOCK.json']),
        ):
            checked_worker([sys.executable, str(study), action], destination, destination / (action + '.log'), required, env)
        if read(destination / 'prepared/PREPARATION.json')['status'] != 'ready':
            raise RuntimeError('Primary preparation blocked')
        if read(destination / 'controls/RESULT.json')['exit_code'] != 0:
            raise RuntimeError('Scientific controls failed')
        for relative, expected_hash in read(ROOT / 'EXPECTED.json')['prepared_coordinates'].items():
            if sha(destination / relative) != expected_hash:
                raise RuntimeError('Prepared coordinate identity differs: ' + relative)
        base = destination / 'runs/replay-001'
        base.mkdir(parents=True)
        for i in range(1, 11):
            query = 'Q' + str(i)
            if i == 3:
                primary = [read(base / q / 'state/summary.json')['call'] for q in ('Q1', 'Q2')]
                if primary != ['open_like', 'closed_like']:
                    raise RuntimeError('Primary stopping rule: secondary queries withheld')
            for kind in ('state', 'burial', 'common') if i <= 2 else ('state', 'burial'):
                output = base / query / kind
                output.parent.mkdir(parents=True, exist_ok=True)
                required = [output / 'summary.json']
                if kind == 'state':
                    required.extend(output / 'engine' / name for name in
                                    ('STATE_ENSEMBLE.json', 'STATE_LAYER.json', 'RUN_MANIFEST.json',
                                     'FRAME_METRICS.tsv', 'CLUSTERS.tsv'))
                else:
                    required.extend(output / ('engine_' + ch) / 'ASSEMBLY_CONTEXT.json' for ch in ('A', 'B'))
                checked_worker([sys.executable, str(study), 'leg', kind, query, str(output)], destination,
                               output.parent / (kind + '.log'), required, env)
        result = compare_run(destination)
        write(destination / 'COMPARISON.json', result)
        if result['passed'] != result['comparison_groups']:
            raise RuntimeError('Measurements differ; retained comparisons do not pass')
        write(destination / 'RESULT.json', {'status': 'pass', 'utc': utc(), 'comparison_groups': result['passed'],
                                           'registration_sha256': sha(destination / 'REGISTRATION.json')})
        print(json.dumps({'status': 'pass', 'comparison_groups': result['passed'], 'output': str(destination)}))
    except Exception as error:
        write(destination / 'RESULT.json', {'status': 'failed', 'utc': utc(), 'error': str(error),
                                           'registration_sha256': sha(destination / 'REGISTRATION.json')})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('verify-kit')
    check = commands.add_parser('check')
    check.add_argument('--source-root', type=Path, default=ROOT)
    check.add_argument('--freesasa', required=True)
    run = commands.add_parser('run')
    run.add_argument('--source-root', type=Path, default=ROOT)
    run.add_argument('--freesasa', required=True)
    run.add_argument('--output', type=Path, required=True)
    compare = commands.add_parser('compare')
    compare.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify_payload()
    if args.command == 'verify-kit':
        print('Kit identities pass; no calculation or runtime check performed')
    elif args.command == 'check':
        verify_sources(args.source_root.resolve())
        runtime(args.freesasa)
        print('Source and runtime preflight pass; no calculation performed')
    elif args.command == 'run':
        execute(args.output.resolve(), args.source_root.resolve(), args.freesasa)
    else:
        result = compare_run(args.output.resolve())
        print(json.dumps(result, indent=2))
        if result['passed'] != result['comparison_groups']:
            raise SystemExit(1)


if __name__ == '__main__':
    main()
