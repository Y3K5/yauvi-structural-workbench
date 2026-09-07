# YAUVI Structural Workbench implementation

The local development preview now packages the browser and scientific engines,
with a shared case service, offline example, durable jobs, revision history and
checksum-verified run reuse. Publication and scientific completeness remain gated.

This is the maintained implementation record for the plan approved on 6 September
2026. The canonical source is this package; the separate GitHub checkout is its
reviewed publication destination. The initial 1,441-file comparison is in
the local baseline record (withheld from Git because it contains machine provenance). Later changes in the
publication checkout were observed; they must be compared before a future sync.
No commit, push, upload, public evidence replacement or history edit is part of
this implementation.

## Personal UX preview

The subsequent [UX preview](docs/UX_PREVIEW.md) adds a quieter four-section interface, real stage tracking, preserved run selection, clearer downloads and the template-download repair. Its separate [local build](artifacts/local-preview-20260906-ux/) supersedes r4 for the interface. Earlier evidence and artifacts below are preserved as historical records; consult the UX verification record for the new checks.

## Read next

- [Installation and offline browser](docs/install.md)
- [Six worked examples](docs/tutorials.md)
- [Maintained backlog](IMPLEMENTATION_BACKLOG.json)
- [Protocols awaiting scientific review](protocols/README.md)
- [Claim-to-evidence ledger](CLAIM_LEDGER.json)
- [Pilot and independent reproduction protocol](docs/PILOT_PROTOCOL.md)
- [Release review and approval boundaries](docs/RELEASE_REVIEW.md)
- [Fresh verification evidence](implementation-evidence/2026-09-06/)

## Separate gate states

The final local Python 3.12 suite passed 560 tests, with two skips and five deselections. The final installed Python 3.10 build also passed 560 tests, with the same two skips and five deselections. Browser/CLI scientific documents and residue tables matched for the bundled offline example. A missing-validation browser example correctly reported scientifically incomplete. Scientific
qualification of the changed code is incomplete; archived passing cases cannot
be relabeled as fresh evidence. Independent reproduction and scientific review
are not recorded. Publication authorization is false.

Membrane orientation remains experimental. The original Linux arm64 regression
is unresolved: the offline container passed the same .90/.80 thresholds and did
not reproduce the original .888 mean Jaccard. This is a useful environment
comparison, not evidence that the defect is fixed.

SF-CSA remains blocked on complete archived ten-proteome recovery. Four ABL
reference structures remain separate from ten held-out cases. New StructQC
modified-residue handling and the conservative ActState occupancy correction
require fresh scientific evaluation before any upgrade in qualification.

## Weekly checkpoint

Update the existing backlog instead of creating a second task list. Record what
works and its version-bound evidence, changes to scientific claims, remaining
failures, user misunderstandings, and the next bounded work package. Yuvraj owns
scientific/product decisions; the software engineering role owns execution
reliability. External people supply independent criticism and usability evidence.
AI review does not fulfill that requirement.

## Local artifacts ready for inspection

The final wheel and source archive are in
[local-preview-20260906-r4](artifacts/local-preview-20260906-r4/).
The preceding preview directories are retained as superseded local builds.
All use the existing development version identifier; compare their exact SHA-256
values rather than treating that shared version string as an identity.
The final archive inventory has 85 wheel files and 91 source-archive files.
[The screening report](artifacts/local-preview-20260906-r4/ARTIFACT_SCREEN.json)
found no matches in its declared path, username, private-project-name and
secret-shaped pattern checks. That is not an exhaustive privacy certification
or a completed license review. No artifact is approved for publication.

The revised paper is [available as source](paper/paper.md), with a separate
`paper-working-draft.pdf`; the original `paper.pdf` remains untouched.
The five-page draft was rendered offline and visually checked. The manuscript
and three reusable figures distinguish historical qualification from new software
verification. Author review, funding/conflicts and AI-version reconciliation
remain pending.

Remaining product acceptance includes broader input/security exercises and
browser/CLI equivalence across all six workflows. Scientific acceptance still
requires the reviewed protocols, recovered SF-CSA universe, relevant fresh panel
runs, independent criticism, and external pilot/reproduction evidence. These
open items prevent a version 1 or submission-ready declaration.

The proposed execution-evidence projection is retained locally for separate review before any replacement of public records. It is excluded from this Git preparation.

The membrane suite also retains one expected failure for its documented geometry limitation; this is separate from successful software tests and is not scientific qualification. The final standalone membrane metadata now labels beta-barrel orientation experimental and non-blocking, consistently with the workbench.
