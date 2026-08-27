#!/usr/bin/env python3
"""Keep every vendored copy of membrane_bilayer.js byte-identical to the canonical one.

    python3 sync_membrane_bilayer.py            # copy canonical -> consumers
    python3 sync_membrane_bilayer.py --check    # report drift, exit 1, change nothing

Why copies exist at all: each consumer emits a SELF-CONTAINED page that has to open over
file:// with no sibling fetches and no cross-project paths. A runtime link would couple
Triple Vax to Protein Platform (and back — memorient already ports numeric modules from
redvax), so the module is inlined at build time instead and the copies are enforced equal.

Drift between these viewers is not hypothetical — it is what this module was created to
fix, so --check is meant to be cheap enough to run on every build.
"""
import hashlib
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent                     # .../Membrane Orientor/memorient
CANONICAL = HERE / "src" / "memorient" / "membrane_bilayer.js"

# Resolved from the workspace root so this keeps working through the YAUVI compatibility
# symlinks (Protein Platform/yauvi -> projects/YAUVI-TDVax/portal/legacy).
WORKSPACE = HERE.parents[2]                                # .../Desktop/YYY
CONSUMERS = [
    WORKSPACE / "Protein Platform" / "yauvi" / "assets" / "membrane-bilayer.js",
    WORKSPACE / "Triple Vax" / "src" / "redvax" / "viz" / "membrane_bilayer.js",
]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    if not CANONICAL.exists():
        print(f"ERROR: canonical missing: {CANONICAL}", file=sys.stderr)
        return 1
    want = sha(CANONICAL)
    check = "--check" in sys.argv
    drift = 0

    print(f"canonical {want[:12]}  {CANONICAL.relative_to(WORKSPACE)}")
    for dest in CONSUMERS:
        rel = dest.relative_to(WORKSPACE)
        if not dest.parent.is_dir():
            print(f"  MISSING DIR  {rel}", file=sys.stderr)
            drift += 1
            continue
        have = sha(dest) if dest.exists() else None
        if have == want:
            print(f"  ok           {rel}")
            continue
        drift += 1
        if check:
            print(f"  DRIFT        {rel}  ({'absent' if have is None else have[:12]})",
                  file=sys.stderr)
        else:
            shutil.copyfile(CANONICAL, dest)
            print(f"  synced       {rel}")

    if check and drift:
        print(f"\n{drift} copy/copies differ — run without --check to sync", file=sys.stderr)
        return 1
    print(f"\n{'all copies match' if not drift else str(drift) + ' synced'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
