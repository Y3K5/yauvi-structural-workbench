"""Screen a resolved payload before it leaves the machine; report rule names and counts only.

Shared by the tester and developer kit builders. Public source carries only generic
rules. Names that are private to one maintainer live in a local denylist that is
never committed, never packaged and never echoed into a findings report:

  1. ``--denylist PATH`` on the builder's command line, else
  2. the file named by ``YAUVI_PRIVATE_DENYLIST``, else
  3. ``~/.config/yauvi/private-denylist.json``.

The file is a JSON list of non-empty strings. A builder refuses to run when none is
found unless the caller opts out explicitly, so the private-name check cannot vanish
silently. The builder's own home directory and account name are added automatically
at run time, so they never need to be written down anywhere.

Pattern screening narrows risk; it does not prove that nothing sensitive remains.
"""
from __future__ import annotations

import getpass
import importlib.util
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DENYLIST = Path("~/.config/yauvi/private-denylist.json")
ENV_DENYLIST = "YAUVI_PRIVATE_DENYLIST"
MAX_EXPANDED_BYTES = 256 * 1024 * 1024
MAX_DEPTH = 3

# One definition of a home path, a secret shape and a placeholder account: the
# repository-wide verifier's. Two scanners that disagree are worse than one.
_spec = importlib.util.spec_from_file_location("_verify_paths", Path(__file__).with_name("verify_no_local_paths_in_evidence.py"))
_verify = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_verify)
HOME_PATH, SECRET, EXEMPT_ACCOUNTS = _verify.HOME_PATH, _verify.SECRET, _verify.EXEMPT_ACCOUNTS

WINDOWS_HOME = re.compile(r"[A-Za-z]:\\+Users\\+([A-Za-z0-9._-]+)", re.IGNORECASE)
SYSTEM_TEMP = re.compile(r"/private/(?:tmp|var)/")
# A separator in a private term matches any run of space, underscore, hyphen, dot or
# URL-encoded space, so "Project 7", "project_7" and "PROJECT-7" are one term.
_SEPARATOR = r"(?:[\s_.\-]|%20)*"


class DenylistError(ValueError):
    pass


def _term_pattern(term: str) -> re.Pattern:
    parts = [re.escape(p) for p in re.split(r"[\s_.\-]+", term.strip()) if p]
    if not parts:
        raise DenylistError("denylist terms must contain a letter or digit")
    return re.compile(_SEPARATOR.join(parts), re.IGNORECASE)


def load_denylist(path: str | os.PathLike | None = None, *, allow_missing: bool = False) -> tuple[list[str], str]:
    """Return (terms, source description). Never returns the path of a file inside the repository."""
    candidates = [("--denylist", path)] if path else []
    if os.environ.get(ENV_DENYLIST):
        candidates.append((ENV_DENYLIST, os.environ[ENV_DENYLIST]))
    candidates.append(("default location", DEFAULT_DENYLIST))
    for label, candidate in candidates:
        file = Path(candidate).expanduser()
        if not file.is_file():
            if label != "default location":
                raise DenylistError(f"private denylist named by {label} does not exist")
            continue
        resolved = file.resolve()
        if resolved == ROOT or ROOT in resolved.parents:
            raise DenylistError("the private denylist must live outside the repository tree")
        try:
            terms = json.loads(resolved.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise DenylistError(f"private denylist from {label} is not readable JSON") from exc
        if not isinstance(terms, list) or not terms or any(not isinstance(t, str) or not t.strip() for t in terms):
            raise DenylistError("private denylist must be a non-empty JSON list of non-empty strings")
        for term in terms:
            _term_pattern(term)
        return terms, f"{label} ({len(terms)} terms)"
    if allow_missing:
        return [], "not used (explicit opt-out)"
    raise DenylistError(
        "no private denylist found; create ~/.config/yauvi/private-denylist.json, set "
        f"{ENV_DENYLIST}, pass --denylist, or pass --no-private-denylist to build without one")


def runtime_terms() -> list[str]:
    """This machine's own home path and account name, unless they are placeholder accounts."""
    terms = []
    home = Path.home()
    if home.name and home.name not in EXEMPT_ACCOUNTS:
        terms.append(str(home))
    try:
        user = getpass.getuser()
    except (KeyError, OSError):
        user = ""
    # Short or common account names ("sam", "admin") would flag ordinary words; the
    # home path above still catches them where it matters.
    if len(user) >= 6 and user not in EXEMPT_ACCOUNTS:
        terms.append(user)
    return terms


def _text_findings(name: str, text: str, private: list[re.Pattern]) -> list[dict]:
    found = []

    def add(category, count):
        if count:
            found.append({"file": name, "category": category, "occurrences": count})

    add("absolute_local_path",
        sum(1 for m in HOME_PATH.finditer(text) if m.group(1) not in EXEMPT_ACCOUNTS)
        + sum(1 for m in WINDOWS_HOME.finditer(text) if m.group(1) not in EXEMPT_ACCOUNTS)
        + len(SYSTEM_TEMP.findall(text)))
    add("credential_shape", len(SECRET.findall(text)))
    for index, pattern in enumerate(private, 1):
        add(f"private_term_{index}", len(pattern.findall(text)))
    return found


def _decode(data: bytes) -> str:
    # UTF-16 files carry a BOM; everything else is read as UTF-8 with replacement so a
    # binary blob cannot raise, and its ASCII runs are still screened.
    if data[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return data.decode("utf-16", errors="replace")
    return data.decode("utf-8", errors="replace")


def _screen(name: str, data: bytes, private: list[re.Pattern], depth: int, budget: list[int]) -> list[dict]:
    budget[0] -= len(data)
    budget[1] += 1
    if budget[0] < 0 or depth > MAX_DEPTH:
        return [{"file": name, "category": "inspection_limit_reached", "occurrences": 1}]
    # The member name is screened too: a path can leak through a file name alone.
    found = _text_findings(name, name + "\n" + _decode(data), private)
    if zipfile.is_zipfile(io.BytesIO(data)):
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            for item in archive.infolist():
                if item.is_dir():
                    continue
                child = f"{name}!{item.filename}"
                member = PurePosixPath(item.filename)
                if (member.is_absolute() or ".." in member.parts or stat.S_ISLNK(item.external_attr >> 16)
                        or item.file_size > budget[0]):
                    found.append({"file": child, "category": "unsafe_archive_member", "occurrences": 1})
                    continue
                found.extend(_screen(child, archive.read(item), private, depth + 1, budget))
    elif name.lower().endswith((".gz", ".tgz", ".tar", ".bz2", ".xz", ".7z", ".rar")):
        found.append({"file": name, "category": "unsupported_archive_requires_review", "occurrences": 1})
    return found


def screen(files: dict[str, bytes], denylist: list[str], denylist_source: str = "") -> dict:
    """Screen every file, recursing into ZIP archives (wheels included). Terms are never reported."""
    private = [_term_pattern(t) for t in list(dict.fromkeys(denylist + runtime_terms()))]
    budget = [MAX_EXPANDED_BYTES, 0]  # bytes left, entries screened
    findings = []
    for name, data in files.items():
        findings.extend(_screen(name, data, private, 0, budget))
    return {
        "files_scanned": len(files),
        "entries_screened_including_archive_members": budget[1],
        "findings": findings,
        "private_denylist": denylist_source or ("in use" if denylist else "not used"),
        "rules": ["absolute_local_path", "credential_shape", "private_term_N",
                  "unsafe_archive_member", "unsupported_archive_requires_review", "inspection_limit_reached"],
        "scope": ("Resolved payload, file names and every ZIP member to depth 3; private terms matched "
                  "case-insensitively across space, underscore, hyphen and dot. Pattern screening is not "
                  "proof that all sensitive content is absent."),
        "external_release_performed": False,
    }
