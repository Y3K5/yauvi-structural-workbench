#!/usr/bin/env python3
"""Record the exact arm64 regression fixture without changing its acceptance bar."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import platform
import sys
from pathlib import Path
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / 'Membrane Orientor/memorient'
sys.path.insert(0, str(ENGINE / 'tests'))
from synthetic import make_barrel
from memorient.contexts import get_context
from memorient.geometry import canonical_rotation
from memorient.orientor import rotation_validate

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('use a new output file; previous investigations are preserved')
    parameters = dict(n_strands=12, strand_len=10, ec_loop_len=10, peri_loop_len=2, seed=0)
    structure = make_barrel(**parameters)
    result = rotation_validate(structure, get_context('gram_negative_om'), n_points=160)
    _, _, frame = canonical_rotation(structure.ca)
    passed = result['mean_jaccard'] >= .90 and min(result['jaccards']) >= .80
    report = dict(schema_version='1.0', fixture=parameters, sasa_points=160,
        thresholds=dict(mean_jaccard_min=.90, individual_jaccard_min=.80),
        test_passed=bool(passed), rotation=result, frame=frame,
        runtime=dict(python=platform.python_version(), system=platform.system(), machine=platform.machine(),
                     numpy=np.__version__, scipy=scipy.__version__,
                     openblas_threads=os.getenv('OPENBLAS_NUM_THREADS', 'unset'),
                     numerical_build=np.__config__.show(mode='dicts')),
        source_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted((ENGINE / 'src').rglob('*.py'))},
        interpretation='Rotation self-consistency only. Passing this fixture does not resolve the historical CI failure or establish agreement with OPM.')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True, default=str)+'\n')
    print(json.dumps(dict(test_passed=bool(passed), mean_jaccard=result['mean_jaccard'], jaccards=result['jaccards'])))
    return 0 if passed else 1

if __name__ == '__main__':
    raise SystemExit(main())
