"""Local developer entry point. No third-party dependency needed for preview."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'source'

def environment():
    env = os.environ.copy()
    roots = sorted((SOURCE / 'software').glob('*/src'))
    roots += [SOURCE / 'software/Membrane Orientor/memorient/src']
    env['PYTHONPATH'] = os.pathsep.join(str(p) for p in roots)
    env['PATH'] = str(Path(sys.executable).parent) + os.pathsep + env.get('PATH', '')
    env['YAUVI_SHOWCASE_DIR'] = str(ROOT / 'preview/showcase/workbench')
    return env

def smoke():
    with tempfile.TemporaryDirectory(prefix='yauvi-dev-smoke-') as temp:
        base = [sys.executable, '-m', 'yauvi_structural_workbench.cli', '--workspace', temp]
        for name, extra in [('complete', []), ('missing', ['--without-validation'])]:
            steps = [(['example', '--analysis', name] + extra, 0),
                     (['analysis', 'validate', '--analysis', name], 0),
                     (['analysis', 'run', '--analysis', name], 0 if name == 'complete' else 1),
                     (['analysis', 'export', '--analysis', name, '--out', str(Path(temp)/(name+'-report'))], 0)]
            for args, expected in steps:
                result = subprocess.run(base + args, env=environment(), text=True, capture_output=True)
                if result.returncode != expected:
                    raise RuntimeError(f'{name} {args[0:2]} returned {result.returncode}; expected {expected}: {result.stderr}')
                if args[:2] == ['analysis', 'run']:
                    record = json.loads(result.stdout)
                    expected_status = 'completed' if name == 'complete' else 'scientifically_incomplete'
                    if record['status'] != expected_status:
                        raise RuntimeError('Unexpected evidence status: '+record['status'])
            assert (Path(temp)/(name+'-report')/'REPORT.html').is_file()
            print(name + ': expected execution status and report export verified')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['serve', 'preview', 'test', 'smoke'])
    args = parser.parse_args()
    if args.action == 'smoke':
        smoke(); return 0
    if args.action == 'preview':
        command = [sys.executable, '-m', 'http.server', '8964', '--bind', '127.0.0.1', '--directory', str(ROOT/'preview')]
        print('Open http://127.0.0.1:8964/START_HERE.html', flush=True)
    elif args.action == 'test':
        command = [sys.executable, '-m', 'pytest', str(SOURCE/'tests'), '-q']
    else:
        command = [sys.executable, '-m', 'yauvi_structural_workbench.cli', '--workspace', str(ROOT/'local-workspace'), 'workbench', 'serve', '--port', '8963', '--label', 'YAUVI Developer Preview']
    return subprocess.call(command, cwd=SOURCE, env=environment())

if __name__ == '__main__':
    raise SystemExit(main())
