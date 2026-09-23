# YAUVI developer handoff

Start with `preview/START_HERE.html` for the product tour. This is a local development handoff for software and UX/UI evaluation, not a production or scientifically qualified release.

## Setup

Use Python 3.12 (supported minimum: 3.10). Extract the whole archive first.

macOS / Linux:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install './source[dev]'
python dev.py test
python dev.py smoke
python dev.py serve
```

Windows PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install './source[dev]'
python dev.py test
python dev.py smoke
python dev.py serve
```

If activation is restricted, invoke `.venv/Scripts/python.exe` directly; do not change system security policy. First installation downloads dependencies from PyPI. Python and dependencies are not bundled. Version ranges are not a lockfile; resolve and record your own test environment.

Open http://127.0.0.1:8963 after `serve` starts. It runs the editable source with a separate `local-workspace` library and public-reference fetching disabled. Stop with Ctrl+C. Restart after edits: the server deliberately pins its assets at startup. Do not run the installed `yauvi` command to check source edits; use `dev.py serve`.

For the welcome page and tour only, `python dev.py preview` serves the supplied preview on http://127.0.0.1:8964/START_HERE.html using only Python's standard library. Direct file opening is also supported by design, but remains a browser compatibility test.

## Where to work

| Location | Purpose |
|---|---|
| `source/software/structural-workbench/src/yauvi_structural_workbench/ui/` | Workbench HTML, CSS, JavaScript and licensed viewer |
| `source/software/platform/src/yauvi_platform/structural_workbench/` | Analysis contracts, evidence handling and orchestration |
| `source/software/*/src/` | Scientific engines and reference adapters |
| `source/showcase/` | Editable showcase template, CSS and JavaScript with frozen public example data |
| `source/tools/build_showcase_snapshot.py` | Rebuild the tour after UI edits, without private datasets or network |
| `source/tests/` | Offline regression tests shipped for this handoff |
| `preview/` | Runnable tester kit, installed-app wheel and validation record |
| `docs/UX_UI_ACCEPTANCE.md` | Tasks, expected behavior and review matrix |
| `docs/ISSUE_TEMPLATE.md` | Reproducible bug-report format |
| `HANDOFF_VALIDATION.md` | Checks and unresolved gaps for this exact handoff |
| `MANIFEST.json`, `CONTENT_SCREEN.json` | File fingerprints and resolved-payload screen |

After showcase edits, run `python source/tools/build_showcase_snapshot.py` from the handoff root. It updates both standalone and workbench tours and the embedded welcome tour. Then restart the workbench. This rebuild changes UI around frozen example data; it does not regenerate scientific measurements. Original scientific acquisition and full qualification campaigns are outside this handoff.

`dev.py smoke` creates two synthetic cases in a temporary folder and checks execution plus report export. Missing validation is expected to produce `scientifically_incomplete`, not a success claim. No public sequences are fetched.

## Review and return

Follow the acceptance checklist, record OS/browser/Python/dependency versions, and attach reproduction steps and sanitized screenshots. Return a patch or changed files with a test summary; do not send environments, caches, the local analysis library, or confidential protein data. Keep source changes separate from scientific data changes. No deployment, telemetry or external publication is part of this handoff.

Use US English: color, centered, license, practice. Preserve uncertainty, exact evidence sources, accessibility, keyboard access, and the distinction between a completed calculation and a biological conclusion.
