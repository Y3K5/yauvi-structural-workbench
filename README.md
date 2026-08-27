# JOSS Publication — Structural Biology

A self-contained, working build of the **YAUVI Structural Biology Platform —
Mark 1** command-line distribution, assembled as the candidate for eventual JOSS
submission.

This folder is a **working copy**, not an authorized publication. Nothing here is
public, no repository or remote exists, and creating one requires explicit
approval.

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
- **495 tests pass, 0 fail** (6 network/adapter deselected, 1 skipped)
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

Assembling a working build does not clear the publication gates. The blocking
items are unchanged:

- **No version control.** The six-month public-history requirement has not begun.
- **Qualification v2 is unexecuted** — 0 of 114 required cases adopted across six
  panels. No scope is scientifically qualified for release.
- **No community channels.** The code of conduct, security, and support routes
  have no monitored contact.
- **The paper still states its own ineligibility** and carries unresolved
  conflict-of-interest and funding statements.
- No independent installation or research-use evidence is recorded.

Passing local tests is not JOSS acceptance and not external scientific
validation. See
[`docs/JOSS_READINESS_REVIEW_2026-08-26.md`](../docs/JOSS_READINESS_REVIEW_2026-08-26.md)
for the full audit and [`CHANGES.md`](CHANGES.md) for what was fixed here.
