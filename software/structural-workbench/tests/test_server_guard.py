"""The local server's request guard.

The workbench serves a browser page on localhost with a user's private case
evidence behind it. These checks are the only thing stopping another page open
in the same browser from reading that evidence or submitting work. Each one is
tested on its own so a future edit cannot drop one quietly.

The guard is exercised directly on a stand-in request rather than over a socket:
it reads only the request headers, the bound port and the session token, so a
real connection would add nothing but flakiness.
"""
from __future__ import annotations

import pytest

from yauvi_structural_workbench.server import Handler

PORT = 8947
TOKEN = "session-token-for-tests"


class StubServer:
    server_port = PORT
    token = TOKEN


class StubRequest:
    """Enough of a request for the guard: headers, the port, and the token."""

    def __init__(self, **headers):
        self.headers = {k.replace("_", "-"): v for k, v in headers.items()}
        self.server = StubServer()


def guard(mutating=False, **headers):
    return Handler._guard(StubRequest(**headers), mutating)


@pytest.mark.parametrize("host", [f"127.0.0.1:{PORT}", f"localhost:{PORT}", f"[::1]:{PORT}"])
def test_the_three_local_host_spellings_are_accepted(host):
    guard(Host=host)


def test_a_request_without_a_host_header_is_refused():
    with pytest.raises(PermissionError, match="invalid local host"):
        guard()


def test_a_request_for_another_address_is_refused():
    """A DNS name pointing at this machine must not reach the server."""
    with pytest.raises(PermissionError, match="invalid local host"):
        guard(Host="workbench.example.com")


def test_a_request_on_another_port_is_refused():
    with pytest.raises(PermissionError, match="invalid local host"):
        guard(Host=f"127.0.0.1:{PORT + 1}")


def test_a_matching_origin_is_accepted():
    guard(Host=f"127.0.0.1:{PORT}", Origin=f"http://127.0.0.1:{PORT}")


def test_another_sites_origin_is_refused():
    """Any other page in the browser is a different origin."""
    with pytest.raises(PermissionError, match="same-origin"):
        guard(Host=f"127.0.0.1:{PORT}", Origin="http://evil.example.com")


def test_an_https_origin_on_the_same_host_is_refused():
    with pytest.raises(PermissionError, match="same-origin"):
        guard(Host=f"127.0.0.1:{PORT}", Origin=f"https://127.0.0.1:{PORT}")


def test_a_cross_site_fetch_is_refused_even_without_an_origin():
    with pytest.raises(PermissionError, match="cross-site"):
        guard(Host=f"127.0.0.1:{PORT}", Sec_Fetch_Site="cross-site")


def test_a_same_origin_fetch_marker_is_accepted():
    guard(Host=f"127.0.0.1:{PORT}", Sec_Fetch_Site="same-origin")


def test_a_change_without_a_session_token_is_refused():
    """Reading is guarded by origin; changing anything needs the session token."""
    with pytest.raises(PermissionError, match="valid local session"):
        guard(mutating=True, Host=f"127.0.0.1:{PORT}")


def test_a_change_with_the_wrong_token_is_refused():
    with pytest.raises(PermissionError, match="valid local session"):
        guard(mutating=True, Host=f"127.0.0.1:{PORT}", X_Yauvi_Token="not-the-token")


def test_a_change_with_a_truncated_token_is_refused():
    with pytest.raises(PermissionError, match="valid local session"):
        guard(mutating=True, Host=f"127.0.0.1:{PORT}", X_Yauvi_Token=TOKEN[:-1])


def test_a_change_with_the_session_token_is_accepted():
    guard(mutating=True, Host=f"127.0.0.1:{PORT}", X_Yauvi_Token=TOKEN)


def test_reading_does_not_require_a_token():
    guard(mutating=False, Host=f"127.0.0.1:{PORT}")
