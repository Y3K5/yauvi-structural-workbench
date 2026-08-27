#!/usr/bin/env python3
"""Fail if the reviewer wheel contains missing, archived, or out-of-scope code."""
from __future__ import annotations

from pathlib import Path
import sys
import zipfile


REQUIRED_ROOTS = {
    "actstate", "assembly_context", "memorient", "sf_csa", "site_context",
    "state_atlas", "structqc", "yauvi_platform", "yauvi_sources",
    "yauvi_structural_workbench",
}
FORBIDDEN_PARTS = {
    "_archive", "oral_atlas", "oral_context", "oral_ecosystem", "cell_fate",
    "dockprep", "pose_evidence", "structdesign", "vaxpipe", "subproteo",
}


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: verify_structural_workbench_wheel.py WHEEL", file=sys.stderr)
        return 2
    wheel = Path(argv[1])
    if not wheel.is_file():
        print(f"wheel does not exist: {wheel}", file=sys.stderr)
        return 2
    with zipfile.ZipFile(wheel) as archive:
        names = [name for name in archive.namelist() if not name.endswith("/")]
    roots = {name.split("/", 1)[0] for name in names if ".dist-info/" not in name}
    missing = sorted(REQUIRED_ROOTS - roots)
    unexpected = sorted(root for root in roots - REQUIRED_ROOTS if not root.endswith(".dist-info"))
    forbidden = sorted({part for name in names for part in FORBIDDEN_PARTS if part in name.lower()})
    platform_files = sorted(name for name in names if name.startswith("yauvi_platform/"))
    if missing or unexpected or forbidden:
        print({"missing": missing, "unexpected": unexpected, "forbidden": forbidden}, file=sys.stderr)
        return 1
    if not platform_files or any(not name.startswith("yauvi_platform/structural_workbench/") for name in platform_files):
        print("the reviewer wheel contains platform code outside structural_workbench", file=sys.stderr)
        return 1
    # setuptools >= 77 stores declared license files under
    # `<dist-info>/licenses/<declared path>`; older releases place a single
    # `<dist-info>/LICENSE`. Both are one complete project license, so accept
    # either rather than failing on the builder's version.
    license_files = [
        name for name in names
        if (".dist-info/licenses/" in name and name.endswith("/yauvi-structural-workbench/LICENSE"))
        or name.endswith(".dist-info/LICENSE")
    ]
    if len(license_files) != 1:
        print(
            "the reviewer wheel does not contain exactly one complete project LICENSE "
            f"(found {len(license_files)}: {license_files})",
            file=sys.stderr,
        )
        return 1
    print(f"verified {wheel.name}: {len(names)} files, canonical structural namespaces only")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
