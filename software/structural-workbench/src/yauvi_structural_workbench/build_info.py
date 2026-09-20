"""Content identities for the local application, independent of install paths."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from importlib.resources import files

from . import __version__

ASSETS = ('index.html', 'app.js', 'style.css', 'vendor/3Dmol-min.js',
          'vendor/membrane-bilayer.js', 'vendor/3DMOL-LICENSE.txt')
PACKAGES = ('yauvi_structural_workbench', 'yauvi_platform.structural_workbench',
            'yauvi_sources', 'structqc', 'memorient', 'state_atlas', 'site_context',
            'actstate', 'assembly_context', 'sf_csa')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def ui_assets():
    root = files('yauvi_structural_workbench').joinpath('ui')
    return {name: root.joinpath(*name.split('/')).read_bytes() for name in ASSETS}


def application_build(assets=None):
    """Identify shipped application files, not scientific qualification or dependencies."""
    assets = ui_assets() if assets is None else assets
    ui_hashes = {name: hashlib.sha256(data).hexdigest() for name, data in assets.items()}
    hashes = {}
    for package in PACKAGES:
        spec = importlib.util.find_spec(package)
        if spec is None or not spec.submodule_search_locations:
            raise RuntimeError(f'Missing application package: {package}')
        for location in spec.submodule_search_locations:
            root = Path(location)
            for path in sorted(root.rglob('*')):
                if not path.is_file() or '__pycache__' in path.parts or path.suffix == '.pyc':
                    continue
                relative = path.relative_to(root).as_posix()
                if package == 'yauvi_structural_workbench' and relative.startswith('ui/'):
                    continue
                hashes[package + '/' + relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    hashes.update({'yauvi_structural_workbench/ui/' + name: value for name, value in ui_hashes.items()})
    return {'application_id': 'yauvi-structural-workbench', 'name': 'YAUVI Structural Workbench',
            'version': __version__, 'channel': 'Development preview',
            'build_id': digest(hashes), 'interface_id': digest(ui_hashes),
            'identity_scope': 'Application files and browser assets; dependency environment and scientific qualification are separate.',
            'distribution': 'Source checkout' if (Path(__file__).resolve().parents[4] / 'pyproject.toml').is_file() else 'Installed package'}


def workspace_identity(workspace):
    return hashlib.sha256(str(Path(workspace).resolve()).encode()).hexdigest()
