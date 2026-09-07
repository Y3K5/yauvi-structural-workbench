#!/usr/bin/env python3
"""Create a portable, inspectable publication projection of execution evidence.

Originals are never changed. Only summary/status documents are selected; pipeline
scratch and logs are not implicitly published. Every replaced output_dir is
accounted for by JSON pointer, original-value digest, and replacement. Unexpected
local paths or secret-shaped values block export for human review.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import re
from pathlib import Path

LOCAL = re.compile(r'/Users/[^/\s"\\]+|/home/[^/\s"\\]+|[A-Za-z]:\\\\Users\\\\|\$HOME', re.I)
SECRET = re.compile(r'(?i)(?:gh[pousr]_[A-Za-z0-9]{20,}|sk-[A-Za-z0-9_-]{24,}|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----)')

def sha(content):
    return hashlib.sha256(content).hexdigest()

def project(document):
    result = copy.deepcopy(document)
    changes = []
    def visit(value, pointer=''):
        if isinstance(value, dict):
            for key, item in list(value.items()):
                position = pointer + '/' + key.replace('~','~0').replace('/','~1')
                if key == 'output_dir' and isinstance(item, str):
                    record_id = value.get('record_id')
                    if not isinstance(record_id, str) or not re.fullmatch(r'[A-Za-z0-9_.-]+', record_id) or record_id in ('.','..'):
                        raise ValueError('output directory has no safe record identifier')
                    if item != record_id:
                        changes.append({'pointer': position, 'original_value_sha256': sha(item.encode()), 'replacement': record_id})
                        value[key] = record_id
                else:
                    visit(item, position)
        elif isinstance(value, list):
            for i, item in enumerate(value): visit(item, pointer + '/' + str(i))
    visit(result)
    encoded = json.dumps(result, sort_keys=True)
    local_username = Path.home().name
    if LOCAL.search(encoded) or SECRET.search(encoded) or (local_username and local_username.casefold() in encoded.casefold()):
        raise ValueError('unresolved local-path or secret-shaped content; review original locally')
    return result, changes

def export(source, target):
    source, target = Path(source).resolve(), Path(target).resolve()
    if source == target or source in target.parents or target in source.parents:
        raise ValueError('projection must be outside the source tree')
    if target.exists(): raise ValueError('projection destination must not already exist')
    paths = sorted(p for p in source.rglob('*.json') if p.name in {'EXECUTION_STATUS.json','EXECUTION_SUMMARY.json'})
    if not paths: raise ValueError('no execution evidence found')
    prepared, ledger = [], []
    for path in paths:
        if path.is_symlink(): raise ValueError('symlink evidence is not allowed')
        original = path.read_bytes(); document, changes = project(json.loads(original))
        content = (json.dumps(document, indent=2, sort_keys=True)+'\n').encode()
        relative = path.relative_to(source)
        prepared.append((relative, content))
        ledger.append({'file': relative.as_posix(), 'original_sha256': sha(original),
                       'projected_sha256': sha(content), 'replacements': changes})
    # All selected documents must pass before any output is created.
    target.mkdir(parents=True)
    for relative, content in prepared:
        p=target/relative; p.parent.mkdir(parents=True, exist_ok=True); p.write_bytes(content)
    report={'schema_version':'1.0','publication_authorized':False,
            'selection':'EXECUTION_SUMMARY.json and EXECUTION_STATUS.json only; logs, raw files and scratch excluded',
            'files':ledger,'unresolved_findings':[]}
    (target/'PUBLICATION_PROJECTION.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    return report

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        report=export(args.source,args.out)
    except (OSError, ValueError) as exc:
        print(f'Export blocked: {exc}');return 2
    print(json.dumps({'files':len(report['files']), 'replacements':sum(len(r['replacements']) for r in report['files']),
                      'publication_authorized':False}));return 0

if __name__=='__main__':raise SystemExit(main())
