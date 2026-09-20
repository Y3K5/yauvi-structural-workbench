"""Open a matching local workbench; never silently reuse another build or library."""
from __future__ import annotations

import errno
import json
import webbrowser
from urllib.error import HTTPError, URLError
from urllib.request import HTTPRedirectHandler, ProxyHandler, build_opener

from yauvi_platform.structural_workbench import AnalysisError
from .build_info import application_build, workspace_identity


class NoRedirects(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise AnalysisError('The local application redirected its identity check. It has not been opened.')


def inspect_running(url):
    try:
        # A system proxy must never receive a request intended for localhost.
        with build_opener(ProxyHandler({}), NoRedirects()).open(url + '/api/build', timeout=3) as response:
            if response.headers.get_content_type() != 'application/json':
                raise AnalysisError('Another application is using this address. Choose another port.')
            raw = response.read(65537)
            if len(raw) > 65536:
                raise AnalysisError('The application at this address returned an invalid identity.')
            value = json.loads(raw)
            if not isinstance(value, dict):
                raise ValueError('identity must be an object')
            return value
    except HTTPError as exc:
        raise AnalysisError('An older workbench or another application is using this address. Close its server before opening this build.') from exc
    except URLError as exc:
        if isinstance(exc.reason, ConnectionRefusedError) or getattr(exc.reason, 'errno', None) == errno.ECONNREFUSED:
            return None
        raise AnalysisError('The local address could not be checked. The existing server has not been changed.') from exc
    except (ValueError, TimeoutError) as exc:
        raise AnalysisError('The application at this address did not return a valid workbench identity.') from exc


def validate_running(current, expected, workspace, allow_reference_fetch,*,label='Local workspace'):
    for field in ('application_id', 'version', 'build_id', 'interface_id'):
        if current.get(field) != expected.get(field):
            raise AnalysisError('A different build is already running here. Finish any active jobs, close that server, then launch again.')
    if current.get('workspace_id') != workspace_identity(workspace):
        raise AnalysisError('This address is serving another analysis library. Use its launcher or choose another port.')
    if current.get('restart_required'):
        raise AnalysisError('Application files changed after this server started. Finish any active jobs and restart the server.')
    if current.get('reference_fetch_enabled') != allow_reference_fetch:
        raise AnalysisError('This server has different reference-retrieval settings. Restart it with the intended settings.')
    if current.get('workspace_label') != label.strip():
        raise AnalysisError('This server has another workspace label. Use its launcher or restart with the intended label.')


def open_workbench(workspace, host, port, allow_reference_fetch=False,*,label='Local workspace'):
    from .server import serve
    hostname = '[::1]' if host == '::1' else host
    url = f'http://{hostname}:{port}'
    current = inspect_running(url)
    if current is not None:
        validate_running(current, application_build(), workspace, allow_reference_fetch,label=label)
        print(f'Opening the matching YAUVI Structural Workbench: {url}', flush=True)
        webbrowser.open(url)
        return 0
    return serve(workspace, host, port, allow_reference_fetch, open_browser=True,label=label)
