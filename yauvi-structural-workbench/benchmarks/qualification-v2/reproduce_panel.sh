#!/bin/sh
# Independent reproduction of a qualification-v2 panel — protocol rule 7.
#
# Run this on a machine that did NOT record the expectations. It executes the
# panel offline against the checksummed sources and compares every case to the
# values recorded on the first machine. Disagreement is the result worth having:
# the membrane stratum passed 16/16 where it was recorded and 9-10/16 elsewhere,
# and that gap is why this step exists.
#
# Usage:  sh reproduce_panel.sh [panel.json] [output-directory]
#         defaults to ADOPTION_DRAFT_ABL.json
# Needs:  python3.11+ with numpy, scipy, gemmi and biopython. No network.
#         The sf-csa panel additionally needs foldseek and diamond on PATH.

set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
DIST=$(CDPATH= cd -- "$HERE/../.." && pwd)
ROOT=$(CDPATH= cd -- "$DIST/.." && pwd)
PANEL=${1:-"$HERE/ADOPTION_DRAFT_ABL.json"}
case "$PANEL" in /*) : ;; *) PANEL="$HERE/$PANEL" ;; esac
OUT=${2:-"${TMPDIR:-/tmp}/panel-reproduction"}
VENV="${TMPDIR:-/tmp}/qualification-repro-venv"

say() { printf '\n== %s\n' "$1"; }
[ -f "$PANEL" ] || { echo "no such panel: $PANEL" >&2; exit 2; }

WORKFLOW=$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["workflow"])' "$PANEL")
say "panel $(basename "$PANEL")  workflow $WORKFLOW"

# Which engine packages this workflow needs. structqc is always needed: every
# chained case regenerates its upstream evidence from the same locked sources.
case "$WORKFLOW" in
  conformational_state) PKGS="structqc state-atlas"; CLIS="structqc state-atlas" ;;
  sf_csa)               PKGS="sf-csa";              CLIS="sf-csa" ;;
  *)                    PKGS="structqc";            CLIS="structqc" ;;
esac

say "1/6 interpreter and packages"
[ -x "$VENV/bin/python" ] || python3 -m venv --system-site-packages "$VENV"
for pkg in $PKGS; do
  "$VENV/bin/pip" install -q --no-deps --no-build-isolation -e "$ROOT/$pkg" 2>/dev/null || {
    echo "could not install $pkg from $ROOT" >&2; exit 2; }
done
"$VENV/bin/python" - <<'PY'
import sys
for mod in ("numpy", "scipy", "gemmi", "Bio"):
    try: __import__(mod)
    except Exception as exc: sys.exit(f"missing dependency {mod}: {exc}")
print("  shared dependencies ok")
PY
PATH="$VENV/bin:$PATH"; export PATH

say "2/6 engines on PATH"
for cli in $CLIS; do
  command -v "$cli" >/dev/null || { echo "$cli is not on PATH" >&2; exit 2; }
  printf '  %s\n' "$(command -v "$cli")"
done
if [ "$WORKFLOW" = "sf_csa" ]; then
  # These are pinned runtimes, and a version mismatch changes the numbers rather
  # than failing loudly. Report them; the panel's own manifest enforces them.
  for tool in foldseek diamond; do
    command -v "$tool" >/dev/null || { echo "$tool is not on PATH" >&2; exit 2; }
    printf '  %s %s\n' "$tool" "$("$tool" version 2>/dev/null | head -1)"
  done
fi

say "3/6 acquiring locked sources"
# The repository declares ships_public_records: false and carries no coordinates,
# so a fresh clone has the panel and none of the data it runs on. Every artifact
# the lock gives a URL for is fetched here and checksummed in the next step;
# anything without a URL is authored and must have come with the clone. Fetching
# is curation, not execution -- the panel itself still runs with no network,
# which is what execution_policy.network_access forbids.
"$VENV/bin/python" - "$HERE" "$PANEL" <<'ACQ'
import json, pathlib, subprocess, sys
here, panel = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
entries = json.loads(panel.read_text()).get("source_lock_additions", [])
need = [e for e in entries if not (here / e["artifact"]).is_file()]
if not need:
    print(f"  all {len(entries)} artifacts already present")
    raise SystemExit(0)
fetchable = [e for e in need if e.get("url")]
orphan = [e for e in need if not e.get("url")]
if orphan:
    print("  absent, and the lock gives no URL for them:", file=sys.stderr)
    for e in orphan:
        print("   ", e["artifact"], file=sys.stderr)
    sys.exit("they are authored files and must come with the repository")
print(f"  fetching {len(fetchable)} of {len(entries)} artifacts")
for e in fetchable:
    dest = here / e["artifact"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"    {e['artifact']}")
    r = subprocess.run(["curl", "-sS", "--http1.1", "--retry", "3", "--max-time", "2400",
                        "-o", str(dest), e["url"]])
    if r.returncode != 0:
        sys.exit(f"could not fetch {e['artifact']} from {e['url']}")
ACQ

say "4/6 source lock"
"$VENV/bin/python" - "$HERE" "$PANEL" <<'PY'
import hashlib, json, pathlib, sys
here, panel = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
entries = json.loads(panel.read_text()).get("source_lock_additions", [])
bad = []
for entry in entries:
    path = here / entry["artifact"]
    if not path.is_file():
        bad.append(f"missing {entry['artifact']}")
    elif hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
        bad.append(f"checksum drift {entry['artifact']}")
if bad:
    sys.exit("source lock failed:\n  " + "\n  ".join(bad))
print(f"  {len(entries)} artifacts verified")
PY

say "5/6 executing the panel"
"$VENV/bin/python" "$HERE/run_execution.py" --panel "$PANEL" --out "$OUT"

say "6/6 what to send back"
printf '  %s\n\n' "$OUT/EXECUTION_STATUS.json"
cat <<'TXT'
  It carries this machine's platform and python version alongside every case.
  Send the file itself, not a summary: a case that passed here and fails there
  is the finding, and the numbers behind it are what identify which case.
TXT
