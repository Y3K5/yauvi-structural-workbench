#!/usr/bin/env python3
"""Refuse any tracked file that publishes a local home directory.

The generator was fixed and the affected evidence was remediated, but neither
of those stops the next one. `run_execution.py` writes panel-relative paths now
(`panel_relative`); a different writer, or a hand-assembled record, can still put
an absolute path into a document that goes to a public repository -- which is
exactly how 51 occurrences of one home directory sat in public history across
three EXECUTION_STATUS.json files for eleven days, in a class already fixed once.

So this checks the resolved set of files actually going out, not a diff. Reading
your own change never catches this, because the leaking field is generated.

Scope: every file `git ls-files` reports, minus binaries. Both a home path and a
secret-shaped string fail the check. A hosted CI runner's own home is exempt --
`/home/runner` names a disposable virtual machine, not a person -- and that
exemption is spelled as an exact path so a real account called `runner` is still
caught.

Usage:  python tools/verify_no_local_paths_in_evidence.py [--root DIR]
Exit:   0 nothing found, 1 something was found, 2 the check could not run.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

#: An account name inside an absolute home path. `$HOME` and `%USERPROFILE%` are
#: deliberately not matched: they are variable references, and a shell step or a
#: scanner's own pattern containing the literal five characters discloses nothing.
#: What discloses is the expansion -- a real name on a real machine.
HOME_PATH = re.compile(r'(?:/Users/|/home/)([A-Za-z0-9._-]+)')
SECRET = re.compile(
    r'(?:gh[pousr]_[A-Za-z0-9]{20,}'
    r'|sk-[A-Za-z0-9_-]{24,}'
    r'|-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----'
    r'|AKIA[0-9A-Z]{16})'
)

#: Account names that name nobody. Two kinds, kept in one list so the exemption is
#: one declared thing rather than scattered special cases:
#:
#:   * hosted-CI accounts, which name a disposable virtual machine -- a runner's
#:     own home legitimately appears in recorded environment metadata;
#:   * documented placeholders in this repository's own tests and scanner
#:     patterns, which exist precisely to prove the check fires.
#:
#: Exempt hits are counted and printed rather than dropped, so a name that joins
#: this list silently is visible in the output. A real person whose account is
#: called `runner` is still disclosed, and this check would miss it -- which is
#: why the list stays short and is reviewed rather than grown by reflex.
EXEMPT_ACCOUNTS = frozenset({
    "runner", "runneradmin",              # GitHub-hosted runners
    "researcher", "example-account",      # platform/tests/test_product_hardening.py
    "example", "user", "someone",         # documentation placeholders
})


def tracked_files(root: Path) -> list[Path]:
    out = subprocess.run(["git", "ls-files", "-z"], cwd=root,
                         capture_output=True, text=True, check=True).stdout
    return [root / name for name in out.split("\0") if name]


def findings(path: Path) -> tuple[list[str], list[str]]:
    """Return (disclosures, exempt hits) for one file."""
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return [], []  # binary or unreadable; nothing textual to disclose

    found: list[str] = []
    exempt: list[str] = []
    for match in HOME_PATH.finditer(text):
        (exempt if match.group(1).casefold() in EXEMPT_ACCOUNTS else found).append(match.group(0))
    for match in SECRET.finditer(text):
        # Never echo the value itself.
        found.append(f"<secret-shaped {match.group(0)[:4]}...>")
    return found, exempt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1],
                        help="Repository root to check. Defaults to this checkout.")
    args = parser.parse_args(argv)

    try:
        files = tracked_files(args.root)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"Cannot enumerate tracked files: {exc}", file=sys.stderr)
        return 2

    total = 0
    exempt_total = 0
    exempt_names: set[str] = set()
    for path in files:
        hits, exempt = findings(path)
        exempt_total += len(exempt)
        exempt_names.update(exempt)
        if not hits:
            continue
        total += len(hits)
        relative = path.relative_to(args.root)
        unique = sorted(set(hits))
        print(f"{len(hits):>4}  {relative}")
        for value in unique[:5]:
            print(f"        {value}")
        if len(unique) > 5:
            print(f"        ... and {len(unique) - 5} more distinct values")

    if exempt_total:
        print(f"\nexempt: {exempt_total} occurrences of "
              f"{', '.join(sorted(exempt_names))} (declared in EXEMPT_ACCOUNTS)")

    if total:
        print(f"\n{total} local-path or secret-shaped occurrences in tracked files.",
              file=sys.stderr)
        print("Fix the writer that produced them, then remediate the recorded output.",
              file=sys.stderr)
        return 1

    print(f"{len(files)} tracked files carry no undeclared local home path "
          f"or secret-shaped value")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
