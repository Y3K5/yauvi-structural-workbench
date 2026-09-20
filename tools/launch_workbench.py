#!/usr/bin/env python3
"""Launch this checkout's browser and engines using an installed dependency environment.

Run with a Python environment containing the project's declared dependencies.
All source roots are resolved from this file; no private path is shipped.
Arguments are forwarded to yauvi, e.g. --workspace ./analyses workbench open.
"""
from __future__ import annotations

import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    if not (3, 10) <= sys.version_info[:2] <= (3, 12):
        raise SystemExit('Use a Python 3.10–3.12 environment with the workbench dependencies installed.')
    try:
        import tomllib
    except ModuleNotFoundError:
        import tomli as tomllib
    metadata = tomllib.loads((ROOT / 'pyproject.toml').read_text())
    paths = [str(ROOT / p) for p in metadata['tool']['setuptools']['packages']['find']['where']]
    sys.path[:0] = paths
    # Scientific child processes must resolve the same source and environment.
    os.environ['PYTHONPATH'] = os.pathsep.join(paths)
    os.environ['PATH'] = str(Path(sys.executable).parent) + os.pathsep + os.environ.get('PATH', '')
    sys.dont_write_bytecode = True
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    from yauvi_structural_workbench.cli import main as cli
    return cli()


if __name__ == '__main__':
    raise SystemExit(main())
