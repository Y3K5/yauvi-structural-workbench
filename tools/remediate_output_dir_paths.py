#!/usr/bin/env python3
"""Retroactively apply the panel_relative() fix to evidence written before it existed.

Three EXECUTION_STATUS.json documents -- structqc, sitecontext and assembly --
were written by a `run_execution.py` that recorded each case's `output_dir` as an
absolute path on the recording machine. They therefore publish a home directory
to a public repository, 51 times between them.

`run_execution.py` has since recorded panel-relative paths (`panel_relative`,
line 99), and every panel executed after that change -- abl and membrane --
carries no absolute path at all. These three were simply never re-run.

This is deliberately not a scrubber. It performs exactly one transformation and
refuses everything else:

  * only `output_dir` is eligible, and only when it begins with the panel root;
  * the remainder must equal `results/execution-<panel>/<record_id>` exactly --
    the string `panel_relative()` would have emitted -- or the file is refused;
  * every other key, in the document and in each case, must be byte-identical
    afterwards, or the file is refused;
  * nothing is written until every selected file has passed every check.

So the result is not an edited record. It is the record the fixed generator
would have written, and the ledger states that claim in a form a reviewer can
recompute: before and after digests, and the replacement made for each record.

The digests of these files are recorded in `RELEASE_STATUS.json`, so applying
this invalidates them by design. Re-record them with

    python tools/verify_release_status_digests.py --update

and review that diff as part of the same change.

Usage:
    python tools/remediate_output_dir_paths.py            # dry run; writes nothing
    python tools/remediate_output_dir_paths.py --apply    # rewrite and write the ledger
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

PANEL_ROOT = Path(__file__).resolve().parent.parent / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2"
PANELS = ("structqc", "sitecontext", "assembly", "abl", "membrane")
LEDGER = PANEL_ROOT / "results" / "OUTPUT_DIR_REMEDIATION_2026-09-09.json"

#: Any absolute home directory, on any platform, in any of the three shapes a
#: recorded path takes here. Kept independent of the account name that actually
#: leaked, so a record written on a different machine is caught just as well.
LOCAL = re.compile(r'/Users/[^/\s"\\]+|/home/[^/\s"\\]+|[A-Za-z]:\\+Users\\+|\$HOME', re.I)


class Refused(Exception):
    """A file did not satisfy every precondition, so nothing is written."""


def plan(panel: str) -> tuple[Path, str, str, list[dict]]:
    """Return (path, original text, remediated text, replacements) or raise Refused."""
    path = PANEL_ROOT / "results" / f"execution-{panel}" / "EXECUTION_STATUS.json"
    if not path.is_file():
        raise Refused(f"{path} does not exist")

    original = path.read_text(encoding="utf-8")
    before = json.loads(original)

    # The panel root as it appears inside this document. Derived from a recorded
    # value rather than from this machine's layout, because the document may have
    # been written somewhere this script is not running.
    prefix = ""
    for group in ("cases", "controls"):
        for case in before.get(group) or []:
            value = case.get("output_dir")
            if isinstance(value, str) and value.startswith("/"):
                marker = f"/results/execution-{panel}/"
                if marker not in value:
                    raise Refused(f"{path}: {value!r} is absolute but not inside this panel's results")
                prefix = value[: value.index(marker) + 1]
                break
        if prefix:
            break

    if not prefix:
        return path, original, original, []

    replacements = []
    for group in ("cases", "controls"):
        for case in before.get(group) or []:
            value = case.get("output_dir")
            if not isinstance(value, str) or not value.startswith("/"):
                continue
            if not value.startswith(prefix):
                raise Refused(f"{path}: {value!r} does not share the panel root {prefix!r}")
            record_id = case.get("record_id")
            if not isinstance(record_id, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", record_id):
                raise Refused(f"{path}: a case has no safe record identifier")
            relative = value[len(prefix):]
            expected = f"results/execution-{panel}/{record_id}"
            if relative != expected:
                raise Refused(f"{path}: {record_id} would become {relative!r}, not {expected!r}")
            replacements.append({
                "record_id": record_id,
                "original_value_sha256": hashlib.sha256(value.encode()).hexdigest(),
                "replacement": relative,
            })

    remediated = original.replace(prefix, "")
    after = json.loads(remediated)

    for key in before:
        if key in ("cases", "controls"):
            continue
        if json.dumps(before[key], sort_keys=True) != json.dumps(after[key], sort_keys=True):
            raise Refused(f"{path}: top-level key {key!r} changed")
    for group in ("cases", "controls"):
        old, new = before.get(group) or [], after.get(group) or []
        if len(old) != len(new):
            raise Refused(f"{path}: {group} changed length")
        for a, b in zip(old, new):
            for key in a:
                if key == "output_dir":
                    continue
                if json.dumps(a[key], sort_keys=True) != json.dumps(b[key], sort_keys=True):
                    raise Refused(f"{path}: {a.get('record_id')} field {key!r} changed")

    if LOCAL.search(remediated):
        found = sorted(set(LOCAL.findall(remediated)))
        raise Refused(f"{path}: still carries local paths after the rewrite: {found}")

    return path, original, remediated, replacements


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--apply", action="store_true",
                        help="Rewrite the files and write the ledger. Without it, nothing is written.")
    args = parser.parse_args(argv)

    prepared, ledger_files = [], []
    try:
        for panel in PANELS:
            path, original, remediated, replacements = plan(panel)
            if not replacements:
                print(f"{panel:12} already panel-relative; nothing to do")
                continue
            prepared.append((path, remediated))
            ledger_files.append({
                "file": path.relative_to(PANEL_ROOT.parent.parent).as_posix(),
                "sha256_before": hashlib.sha256(original.encode()).hexdigest(),
                "sha256_after": hashlib.sha256(remediated.encode()).hexdigest(),
                "home_path_occurrences_removed": len(LOCAL.findall(original)),
                "replacements": replacements,
            })
            print(f"{panel:12} {len(replacements):>2} output_dir values -> panel-relative, "
                  f"{len(LOCAL.findall(original)):>2} home-path occurrences removed")
    except Refused as exc:
        print(f"\nRefused, nothing written: {exc}", file=sys.stderr)
        return 2

    total = sum(f["home_path_occurrences_removed"] for f in ledger_files)
    if not prepared:
        print("\nNo file needed remediation.")
        return 0

    if not args.apply:
        print(f"\nDry run. {len(prepared)} files would change, removing {total} home-path "
              f"occurrences. Re-run with --apply to write them.")
        return 0

    for path, remediated in prepared:
        path.write_text(remediated, encoding="utf-8")
    LEDGER.write_text(json.dumps({
        "schema_version": "1.0",
        "date": "2026-09-09",
        "what": "Retroactive application of run_execution.py's panel_relative() fix to evidence "
                "written before that fix existed.",
        "why": "These documents recorded output_dir as an absolute path on the recording machine, "
               "publishing a home directory to a public repository. Panels executed after the fix "
               "carry no absolute path; these were never re-run.",
        "method": "Only output_dir was rewritten, only where the remainder equalled "
                  "results/execution-<panel>/<record_id> exactly. No measurement, expectation, "
                  "verdict or checksum was touched. See tools/remediate_output_dir_paths.py.",
        "digests_to_rerecord": "yauvi-structural-workbench/RELEASE_STATUS.json",
        "home_path_occurrences_removed": total,
        "files": ledger_files,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"\nRewrote {len(prepared)} files, removed {total} home-path occurrences.")
    print(f"Ledger: {LEDGER}")
    print("Next: python tools/verify_release_status_digests.py --update")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
