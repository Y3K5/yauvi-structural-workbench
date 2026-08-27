#!/usr/bin/env bash
#
# Run the SF-CSA offline fixture end to end.
#
# Two scenarios, both driven entirely by stub executables on PATH -- no
# foldseek, no diamond, no network, no reference database:
#
#   main  a release that passes its own audit, exercising five decision paths
#   trap  a release that FAILS its audit, because a hit was promoted to
#         probable_same_function while its PDB title carries a trap substring
#
# The trap scenario is the point of having two: a fixture whose audit always
# returns clean does not demonstrate that the audit works.
#
# Usage:
#   ./run_fixture.sh                  # run both scenarios, compare against golden/
#   ./run_fixture.sh --update-golden  # rewrite golden/ from this run
#   PY=/path/to/python ./run_fixture.sh
#
# PY defaults to `python3`. Point it at an interpreter that has sf-csa
# importable (e.g. the repo venv).

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-python3}"
UPDATE_GOLDEN=0
[[ "${1:-}" == "--update-golden" ]] && UPDATE_GOLDEN=1

WORK="$(mktemp -d "${TMPDIR:-/tmp}/sfcsa-fixture.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

# The stubs must win over any real foldseek/diamond that happens to be
# installed, so they are prepended, not appended.
export PATH="$HERE/stub_bin:$PATH"

echo "== regenerating inputs =="
"$PY" "$HERE/build_inputs.py" >/dev/null
"$PY" "$HERE/build_inputs.py" --check

fail=0

run_scenario() {
    local scenario="$1" queries="$2" databases="$3" expect="$4"
    local out="$WORK/$scenario"

    echo
    echo "== scenario: $scenario (audit expected to $expect) =="
    SFCSA_FIXTURE_SCENARIO="$scenario" "$PY" -m sf_csa.cli run \
        --queries "$HERE/inputs/$queries" \
        --databases "$HERE/inputs/$databases" \
        --output "$out"

    set +e
    SFCSA_FIXTURE_SCENARIO="$scenario" "$PY" -m sf_csa.cli verify \
        --output "$out" --databases "$HERE/inputs/$databases"
    local rc=$?
    set -e

    if [[ "$expect" == "pass" && $rc -ne 0 ]]; then
        echo "FIXTURE BROKEN: $scenario audit was expected to pass, exited $rc"
        fail=1
    fi
    if [[ "$expect" == "fail" && $rc -eq 0 ]]; then
        echo "FIXTURE BROKEN: $scenario audit was expected to fail, exited 0"
        echo "  (the title trap did not catch the promoted hit -- the audit is not biting)"
        fail=1
    fi

    # Golden comparison happens on the canonical form, because a raw release
    # embeds its output directory name and the absolute proteome paths.
    local canon="$WORK/canon-$scenario"
    "$PY" "$HERE/canonicalise.py" --release "$out" --dest "$canon" --fixture-root "$HERE" >/dev/null

    if [[ $UPDATE_GOLDEN -eq 1 ]]; then
        rm -rf "$HERE/golden/$scenario"
        mkdir -p "$HERE/golden"
        cp -R "$canon" "$HERE/golden/$scenario"
        echo "golden/$scenario updated"
    elif [[ -d "$HERE/golden/$scenario" ]]; then
        if diff -r "$HERE/golden/$scenario" "$canon"; then
            echo "golden/$scenario: match"
        else
            echo "GOLDEN DIFF in $scenario (above)"
            fail=1
        fi
    else
        echo "no golden/$scenario yet; run with --update-golden"
        fail=1
    fi
}

run_scenario main query_manifest.json database_manifest.json pass
run_scenario trap query_manifest_trap.json database_manifest_trap.json fail

echo
if [[ $fail -eq 0 ]]; then
    echo "SF-CSA fixture OK: main release audits clean, trap release is caught, both match golden."
else
    echo "SF-CSA fixture FAILED (see above)."
fi
exit $fail
