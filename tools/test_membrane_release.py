"""Distribution screen controls use invented paths and secrets only."""
import importlib.util
from pathlib import Path
import zipfile

spec=importlib.util.spec_from_file_location('membrane_release',Path(__file__).with_name('build_membrane_release.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_real_paths_and_credentials_are_flagged_but_pattern_syntax_is_not(tmp_path):
    wheel=tmp_path/'invented.whl'
    with zipfile.ZipFile(wheel,'w') as archive:
        archive.writestr('memorient/receipt.json','{"path":"'+str(Path('/').joinpath('Users','invented-user','private','file'))+'"}')
        archive.writestr('memorient/secret.txt','ghp_'+'x'*32)
        archive.writestr('memorient/catalogue/sources/data.bin',b'invented')
        archive.writestr('memorient/rules.py',r"rb'/(?:Users|home)/[A-Za-z0-9._-]+/'")
    report=module.inspect_archive(wheel)
    assert {r['rule'] for r in report['findings']}=={'personal_absolute_path','credential_shape','excluded_biological_payload'}
    assert all(r['member']!='memorient/rules.py' for r in report['findings'])


def test_workbench_and_standalone_ship_the_same_memorient_assets():
    try:import tomllib
    except ImportError:import tomli as tomllib
    home=Path(__file__).resolve().parents[1]
    configs=[tomllib.loads((home/'pyproject.toml').read_text()),
             tomllib.loads((home/'software/membrane-orientor/memorient/pyproject.toml').read_text())]
    lists=[set(c['tool']['setuptools']['package-data']['memorient']) for c in configs]
    assert lists[0]==lists[1]
    assert {'fold_bridge.js','viewer_assets/workbench.js','viewer_assets/3DMOL-LICENSE.txt','catalogue/metadata.json'}<=lists[0]
    assert all('sources/' not in p and 'model-v2' not in p and 'starter.json' not in p for p in lists[0])
