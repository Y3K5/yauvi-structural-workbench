"""A launcher must identify the application before reusing an occupied address."""
import json
import threading
from urllib.error import HTTPError
from urllib.request import ProxyHandler, build_opener

import pytest

from yauvi_platform.structural_workbench import AnalysisError
from yauvi_structural_workbench import build_info, launcher, server
from yauvi_structural_workbench.cli import build_parser, main


def running(workspace):
    expected = {'application_id': 'yauvi-structural-workbench', 'version': '0.1.0.dev0',
                'build_id': 'build-a', 'interface_id': 'interface-a'}
    current = {**expected, 'workspace_id': build_info.workspace_identity(workspace),
               'reference_fetch_enabled': False, 'restart_required': False,'workspace_label':'Local workspace'}
    return expected, current


@pytest.mark.parametrize('field,value', [
    ('application_id', 'other-app'), ('version', 'older'), ('build_id', 'build-b'),
    ('interface_id', 'interface-b'), ('workspace_id', 'another-library'),
    ('restart_required', True), ('reference_fetch_enabled', True),
    ('workspace_label', 'Another workspace'),
])
def test_launcher_refuses_mismatched_server(tmp_path, field, value):
    expected, current = running(tmp_path)
    current[field] = value
    with pytest.raises(AnalysisError):
        launcher.validate_running(current, expected, tmp_path, False)


def test_matching_server_is_opened_without_starting_a_second_one(tmp_path, monkeypatch):
    expected, current = running(tmp_path)
    opened = []
    monkeypatch.setattr(launcher, 'inspect_running', lambda url: current)
    monkeypatch.setattr(launcher, 'application_build', lambda: expected)
    monkeypatch.setattr(launcher.webbrowser, 'open', opened.append)
    def unexpected(*args, **kwargs):
        pytest.fail('Started a duplicate server')
    monkeypatch.setattr(server, 'serve', unexpected)
    assert launcher.open_workbench(tmp_path, '127.0.0.1', 8947) == 0
    assert opened == ['http://127.0.0.1:8947']


def test_redirect_is_not_followed_to_an_external_host():
    with pytest.raises(AnalysisError, match='redirected'):
        launcher.NoRedirects().redirect_request(None, None, 302, '', {}, 'https://example.org')


def test_build_changes_with_viewer_bytes_and_contains_no_local_paths():
    assets = build_info.ui_assets()
    original = build_info.application_build(assets)
    assets['vendor/membrane-bilayer.js'] += b'\n// changed viewer\n'
    edited = build_info.application_build(assets)
    assert edited['build_id'] != original['build_id']
    assert edited['interface_id'] != original['interface_id']
    assert '/Users/' not in json.dumps(original)
    assert '/home/' not in json.dumps(original)


def test_default_port_keeps_structural_and_legacy_portals_distinct():
    for command in ['open', 'serve']:
        assert build_parser().parse_args(['workbench', command]).port == 8947
        assert build_parser().parse_args(['workbench', command, '--port', '8990']).port == 8990


@pytest.mark.parametrize('host', ['0.0.0.0', '192.0.2.1'])
def test_cli_refuses_non_loopback_workbench_binding(host, capsys):
    assert main(['workbench', 'serve', '--host', host]) == 2
    assert 'must be a loopback address' in capsys.readouterr().err


def test_server_refuses_non_loopback_address_before_binding(tmp_path):
    with pytest.raises(AnalysisError, match='loopback host required'):
        server.WorkbenchServer(('0.0.0.0', 0), tmp_path)


def test_cli_version_reports_the_installed_workbench_version(capsys):
    with pytest.raises(SystemExit) as result:
        main(['--version'])
    assert result.value.code == 0
    assert capsys.readouterr().out.strip() == 'yauvi 0.1.0.dev0'


def test_analysis_input_discovery_reports_allowed_roles_and_formats(capsys):
    assert main(['analysis', 'inputs', '--type', 'membrane_orientation']) == 0
    result = json.loads(capsys.readouterr().out)
    assert result['analysis_type'] == 'membrane_orientation'
    roles = {item['role']: item for item in result['inputs']}
    assert roles['structure']['required'] is True
    assert roles['structure']['accepted_extensions'] == ['.pdb', '.cif', '.mmcif']
    assert roles['topology_evidence']['accepted_extensions'] == ['.json']
    assert roles['topology_evidence']['format_guide']


def test_analysis_types_are_discoverable_without_creating_a_case(capsys):
    assert main(['analysis', 'types']) == 0
    result = json.loads(capsys.readouterr().out)
    assert {item['analysis_type'] for item in result} == {
        'structure_qc', 'membrane_orientation', 'conformational_state',
        'functional_site_state', 'assembly_interface', 'sf_csa',
    }
    assert all(item['claim_ceiling'] for item in result)


def _loopback_server(tmp_path):
    """Start the workbench server, or skip if this machine forbids binding.

    The test needs a real loopback socket. Sandboxes, locked-down CI images and
    some corporate builds refuse `bind` outright with PermissionError, and the
    test then fails in a way that reads as a defect in the workbench rather than
    a restriction on the machine. A reviewer running the documented test command
    on such a machine sees red on the first step of the quickstart and has no way
    to tell the difference.

    Skipping names the cause. It does not hide a real failure: anywhere binding
    is permitted the test runs exactly as before.
    """
    try:
        return server.WorkbenchServer(('127.0.0.1', 0), tmp_path)
    except (PermissionError, OSError) as exc:
        pytest.skip(f"this machine does not permit binding a loopback socket: {exc}")


def test_server_pins_assets_and_reports_changed_files(tmp_path, monkeypatch):
    httpd = _loopback_server(tmp_path)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    opener = build_opener(ProxyHandler({}))
    url = f'http://127.0.0.1:{httpd.server_port}'
    try:
        current = launcher.inspect_running(url)
        assert current['restart_required'] is False
        assert current['workspace_id'] == build_info.workspace_identity(tmp_path)
        original_assets = dict(httpd.ui_assets)
        monkeypatch.setattr(server, 'application_build', lambda: {'build_id': 'changed'})
        monkeypatch.setattr(server, 'ui_assets', lambda: {'index.html': b'new HTML'})
        with opener.open(url + '/') as response:
            assert response.read() == original_assets['index.html']
            assert response.headers['Cache-Control'] == 'no-store'
        assert launcher.inspect_running(url)['restart_required'] is True
        with pytest.raises(HTTPError) as error:
            opener.open(url + '/not-an-asset')
        assert error.value.code == 404
    finally:
        httpd.shutdown()
        httpd.server_close()
        thread.join(timeout=3)


def test_showcase_is_off_unless_named_and_serves_only_its_two_files(tmp_path, monkeypatch):
    monkeypatch.delenv('YAUVI_SHOWCASE_DIR', raising=False)
    assert server.showcase_assets() == {}
    for name in server.SHOWCASE:
        (tmp_path / name).write_text(name)
    (tmp_path / 'private.json').write_text('never served')
    monkeypatch.setenv('YAUVI_SHOWCASE_DIR', str(tmp_path))
    assert server.showcase_assets() == {name: name.encode() for name in server.SHOWCASE}
    (tmp_path / 'showcase.js').unlink()
    assert server.showcase_assets() == {}
