# YAUVI Structural Workbench

Local, evidence-bounded structural protein analysis with deterministic reports
and provenance. Ten installable Python packages covering six structural-analysis
workflows, taking protein coordinates to inspectable, checksum-bound evidence.

**This is a pre-public scientific build in open development.** It is not a
released, qualified, or peer-reviewed tool, and it is being prepared as a
candidate for eventual submission to the Journal of Open Source Software. Read
the status section below before relying on any output.

## What this is

Ten installable Python packages covering six structural-analysis workflows,
plus the documentation, paper, community files, benchmarks, and evidence
showcases that a JOSS reviewer would receive.

Nine of the ten carry a workflow or the case store. The tenth, `structprep`, is
**infrastructure**: it prepares coordinates for the others and is not a seventh
analysis scope. It answers no scientific question on its own and is not
release-blocking.

| Package | CLI | Workflow |
|---|---|---|
| `yauvi-structural-workbench` | `yauvi` | Analysis case store |
| `yauvi-structqc` | `structqc` | Coordinate trust |
| `memorient` | `memorient` | Membrane orientation |
| `yauvi-state-atlas` | `state-atlas` | Conformational state |
| `yauvi-site-context` | `site-context` | Functional site |
| `actstate` | `actstate` | Activity state |
| `yauvi-assembly-context` | `assembly-context` | Assembly interface |
| `sf-csa` | `sf-csa` | Structure/sequence function comparison |
| `yauvi-sources` | `yauvi-fetch` | Registered source acquisition |
| `structprep` | `structprep` | Structure preparation (infrastructure) |

## Current local verification

These records describe separate checks on the local candidate; their totals and
interpretations are not combined:

- The offline suite reports **684 passed, 15 skipped, 5 deselected, 0 failed,
  and 0 errors** on macOS arm64 / Python 3.12.0. See
  [`OFFLINE_TESTS.json`](evidence/preparation-2026-09-20/OFFLINE_TESTS.json).
- A separate public-case safety and input-validation run reports **5 passed**;
  it is not added to the full offline-suite count. See
  [`PUBLIC_CASE_SAFETY_TESTS.json`](evidence/preparation-2026-09-20/PUBLIC_CASE_SAFETY_TESTS.json).
- The fresh wheel was installed in an isolated environment. The installed
  distribution reproduced the public-data case twice with 12 compared files
  byte-identical. A browser check completed the synthetic StructQC flow with
  zero console errors. See
  [`CANDIDATE_VERIFICATION.json`](evidence/preparation-2026-09-20/CANDIDATE_VERIFICATION.json).
- An earlier README snapshot reported 590 passing tests and referred to a
  different build. That count is historical and superseded by the dated suite
  record above; it is not evidence about the present candidate.

Start with [`yauvi-structural-workbench/START_HERE.md`](yauvi-structural-workbench/START_HERE.md),
then [`docs/quickstart.md`](yauvi-structural-workbench/docs/quickstart.md) and
[`docs/cli-reference.md`](yauvi-structural-workbench/docs/cli-reference.md).

## Running the tests

```bash
python -m pip install -e ".[dev]"
python tools/run_structural_workbench_tests.py
```

**Use an interpreter whose architecture matches the machine.** On Apple Silicon
that means an arm64 Python 3.10–3.12:

```bash
python -c "import platform; print(platform.machine())"   # expect arm64, not x86_64
```

A mismatched interpreter — most easily an x86_64 conda running under Rosetta —
can acquire wheels whose compiled extensions will not load:

```
ImportError: dlopen(.../Bio/PDB/ccealign.cpython-312-darwin.so):
  (mach-o file, but is an incompatible architecture (have 'arm64', need 'x86_64'))
```

The runner's preflight catches declared-version drift and refuses rather than
reporting failures that describe the environment. An architecture mismatch is
not a version mismatch, so it slips past that guard and lands later as roughly
two dozen pytest collection errors that look like broken code and are not.

`ModuleNotFoundError: No module named 'yauvi_platform'` (or `yauvi_sources`) has
the same character: it means the editable install did not take, not that a module
is missing. Both failure modes are environment reports, not test results.

## Boundary

The current local development build packages both the command-line tools and a
standalone loopback browser. Start it with `yauvi --workspace ./my-analysis workbench open`, which starts the application and opens a browser on it; `workbench serve` starts the same application without opening one.
Both interfaces use the same structural analysis store and scientific engines.
The published commit recorded in the implementation baseline predates this extraction;
publication of these changes remains pending review.

Also excluded: private projects and campaign data, the 12 MB of downloaded
third-party benchmark coordinates (the source lock ships, the files do not), and
the out-of-scope module directories.

## Status — not submission-eligible

Assembling a working build does not clear the publication gates. What still
blocks submission:

- **The public-history gate remains open.** Public development began on
  2026-08-27; more than six months of active development would be reached only
  after 2027-02-27, conditional on continued substantive public activity.
- **The changed local candidate is not scientifically qualified.** Collection
  2.11 adopts 94/110 records. The retained single-machine summary reports
  67/110 across five executed panels and predates recorded SF-CSA execution.
  Separately, public CI run [35295451246](https://github.com/Y3K5/yauvi-structural-workbench/actions/runs/35295451246)
  passed its configured blocking gate on commit
  `7981148e70c5eaa6424e3608f09cd65bcfb35834`; this is historical CI job
  evidence, not qualification of the changed candidate. Membrane coverage is
  beta-barrel-only and research-only/nonblocking. These records do not form a
  combined accuracy score; see
  [`QUALIFICATION_RECONCILIATION.json`](evidence/preparation-2026-09-20/QUALIFICATION_RECONCILIATION.json).
- A fresh installed-wheel case reproduction and synthetic StructQC browser flow
  are recorded, but independent-human research-use interpretation remains
  pending. Software checks and reproducible output do not establish biological
  validity. Author review of scientific claims and journal-facing disclosures
  remains open.
- Both interpretation defects the pre-public audit recorded are now closed.
  SF-CSA computed reciprocal-best-hit after structural classification, so
  `probable_same_function` was unreachable end-to-end; fixed 2026-09-01
  ([`sf-csa/CHANGES.md`](software/sf-csa/CHANGES.md)). ActState reached
  `active_site_disrupted` from membership in a broad residue set, which is not a
  position-specific chemistry test; narrowed 2026-09-02 to require an expected
  residue from a validated reference
  ([`activity-state/CHANGES.md`](software/activity-state/CHANGES.md)). ActState's
  occupancy caveat is unchanged and still open: a non-solvent heteroatom is
  detected, but its identity is not proven against the declared cofactor.

Community channels are in place — see [`CONTRIBUTING.md`](CONTRIBUTING.md),
[`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md),
[`SECURITY.md`](SECURITY.md), and
[`SUPPORT.md`](SUPPORT.md). This is a
single-maintainer project with no conduct committee and no response-time
commitment, which those documents state directly rather than imply.

Passing local tests is not JOSS acceptance and not external scientific
validation. See [`CHANGES.md`](CHANGES.md) for how this distribution was assembled and what
was corrected in the process.
