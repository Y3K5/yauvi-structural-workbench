# Handoff validation — 2026-09-19

Ready for developer review and a guided usability pilot. This is not a production release or a claim of general scientific validity.

## Verified on the extracted handoff

- macOS / Apple silicon, Python 3.12, existing local scientific dependencies.
- `python dev.py test`: **81 passed**, including the real local-server check. An initial sandbox run skipped that one test; the run with loopback binding allowed passed all 81.
- `python dev.py smoke`: complete synthetic Structure QC and deliberately missing-validation variants returned their expected evidence states and exported HTML reports.
- The extracted source built a wheel successfully; an independent wheel-layout verifier found 95 files and only the expected structural namespaces.
- Snapshot rebuild reproduced the standalone tour and embedded welcome page from packaged assets without private source datasets or network access. The script supports updating the workbench tour as well.
- Builder verifies every archived file against MANIFEST.json and checks ZIP integrity. CONTENT_SCREEN.json covers source, preview and nested wheel members.
- Content screening found two deliberate generic path-prefix assertions in a test file. They contain no account identity or full personal path and were retained; no unresolved findings remain. Pattern screening is not proof that every possible sensitive string is absent.

## UI evidence and remaining review

Earlier browser checks covered workbench welcome/help, incomplete-evidence findings, the workbench-to-showcase link, showcase tab navigation, animation start and US spelling. These checks are recorded in preview/VALIDATION.md. The developer checklist starts untested for independent review.

The welcome button now targets its embedded tour, avoiding a sibling-file dependency. Its anchor and embedded document were checked structurally. Direct-file rendering was blocked by the browser tool's URL policy and has not been independently verified. Small-screen layout, keyboard-only full journeys, screen readers, color contrast, feedback downloads and browser compatibility need developer review.

Fresh dependency installation, Python versions other than 3.12, Windows and Linux execution were not tested here. Dependencies use version ranges, not a lockfile. Optional scientific runtimes and databases are not bundled; original scientific acquisition/qualification campaigns are outside this source handoff. All six guided panels are available, but that does not mean arbitrary new inputs are qualified across all six workflows.

No user outcomes, independent biological validation, external upload or deployment are claimed. The archive checksum identifies the handoff. Developer edits and local runs will change files or create local data; do not redistribute a working directory without a new content review.
