# YAUVI Structural Biology Platform — Mark 1

Local, evidence-bounded structural protein analysis with deterministic reports
and provenance. Nine installable Python packages covering six structural-analysis
workflows, taking protein coordinates to inspectable, checksum-bound evidence.

**This is a pre-public scientific build in open development.** It is not a
released, qualified, or peer-reviewed tool, and it is being prepared as a
candidate for eventual submission to the Journal of Open Source Software. Read
the status section below before relying on any output.

## What this is

Nine installable Python packages covering six structural-analysis workflows,
plus the documentation, paper, community files, benchmarks, and evidence
showcases that a JOSS reviewer would receive.

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

## Verified working

Every claim below was executed in this folder, offline:

- `pip install -e ".[dev]"` succeeds; all nine console scripts land on `PATH`
- **526 tests pass, 0 fail** (6 network/adapter deselected, 1 skipped)
- A real StructQC analysis runs to completion and is **byte-identical across two runs**
- The fail-closed path exits `1` and names its missing evidence rather than scoring around it
- The wheel builds offline and contains only canonical structural namespaces
- Installed from that wheel into a fresh environment in an unrelated directory,
  every CLI works and a full analysis reproduces the documented input digest
  `a598a520…` with no absolute-path leakage
- Both evidence showcases rebuild and re-verify; all checksums intact

Start with [`yauvi-structural-workbench/START_HERE.md`](yauvi-structural-workbench/START_HERE.md),
then [`docs/quickstart.md`](yauvi-structural-workbench/docs/quickstart.md) and
[`docs/cli-reference.md`](yauvi-structural-workbench/docs/cli-reference.md).

## Boundary

This is the **command-line** distribution. The loopback browser workbench is
deliberately excluded: its controller imports private control-plane modules
(`yauvi_platform.oral_atlas`, `yauvi_platform.protein_case`) that sit outside the
published boundary. Only `yauvi_platform.structural_workbench` ships, and it has
no sibling imports.

Also excluded: private projects and campaign data, the 12 MB of downloaded
third-party benchmark coordinates (the source lock ships, the files do not), and
the out-of-scope module directories.

## Status — not submission-eligible

Assembling a working build does not clear the publication gates. What still
blocks submission:

- **Public development has only just begun** (first public commit 2026-08-27).
  JOSS expects sustained public history, tagged releases, and evidence of
  independent use. None of that exists yet.
- **No scope is scientifically qualified.** Qualification v2 requires 110 cases
  across six panels. Four panels are adopted and executed and **two are not** --
  ABL StateAtlas and SF-CSA, both release-blocking -- so the collection's
  composition state is `blocked_panel_incomplete` and no scope can be qualified
  whatever the executed panels report. Of the five release-blocking scopes,
  three are adopted and pass their predeclared gates offline on six OS/Python
  combinations in CI: StructQC 16/16 with 2 controls, site-context 16/16 with 1
  control, assembly-context 16/16. Membrane orientation executes at 5/16 against
  the accuracy gate collection 2.3 added, and collection 2.4 moved it to
  non-blocking, research-only, both strata: Mark 1 makes no accuracy claim for
  it. **No scope has completed the independent second-machine reproduction
  gate**, which is required on its own and is untouched by any of the above.
  Counts here come from
  [`EXECUTION_SUMMARY.json`](yauvi-structural-workbench/benchmarks/qualification-v2/results/EXECUTION_SUMMARY.json),
  derived from the executed evidence rather than typed. The historical v1
  collection passed four public cases and left two partial; those are named
  cases, not workflow-general accuracy evidence.
- **The paper still states its own ineligibility** and carries unresolved
  conflict-of-interest and funding statements.
- No independent installation or research-use evidence is recorded.
- Both interpretation defects the pre-public audit recorded are now closed.
  SF-CSA computed reciprocal-best-hit after structural classification, so
  `probable_same_function` was unreachable end-to-end; fixed 2026-09-01
  ([`sf-csa/CHANGES.md`](sf-csa/CHANGES.md)). ActState reached
  `active_site_disrupted` from membership in a broad residue set, which is not a
  position-specific chemistry test; narrowed 2026-09-02 to require an expected
  residue from a validated reference
  ([`activity-state/CHANGES.md`](activity-state/CHANGES.md)). ActState's
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
