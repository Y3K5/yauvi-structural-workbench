# JOSS preparation checklist

This checklist separates work completed locally from evidence that can exist
only after a public repository and independent research use exist.

## Implemented locally

- One root reviewer distribution with canonical structural package namespaces.
- Root and standalone wheels build offline; the root wheel passes an isolated
  smoke install using the current dependency set and can create an analysis case.
- Independently runnable StructQC, MembraneOrient, StateAtlas, SiteContext,
  ActState, AssemblyContext, and SF-CSA command lines.
- Six task-first builders with input purpose, formats, absence effects,
  templates, official source links, claim ceilings, and readiness disclosures.
- Explicit accession acquisition behind `--allow-reference-fetch`; acquisition
  never adopts a file into an analysis automatically.
- Content-addressed ingestion, checksum verification, preflight, registered
  command execution, deterministic report bundles, and print CSS.
- Offline software baseline recorded in [BASELINE.json](BASELINE.json).
- Checksum-locked public qualification runner and evidence report in
  [benchmarks/qualification-v1/](benchmarks/qualification-v1/README.md).
- Frozen Qualification v2 scope, stratum, split, evidence, and tolerance
  specification with a deterministic fail-closed panel audit in
  [benchmarks/qualification-v2/](benchmarks/qualification-v2/README.md).
- Shareable public-safe microsite with file-role guidance, six synthetic
  demonstrations, six separate public qualification narratives, biological
  context, raw evidence, reviewer quickstart, and a gated publication roadmap.
- Public cases passed for StructQC, functional-site mapping, AssemblyContext,
  and SF-CSA; the current Python 3.12.7 reviewer selection reports 494 passed
  and 6 network/adapter tests deselected.
- FreeSASA (version not captured by the runner), Foldseek 10.941cd33, and
  DIAMOND 2.1.11 were invoked in the
  local public qualification run.
- Apache-2.0 text, citation metadata, contribution/support/security/governance
  documents, paper draft, benchmark plan, and CI definitions.
- Historical privacy-minimized ZIP install tests were run on macOS arm64 and
  Linux arm64. No `TEST_RESULT.json` from those runs is committed, and
  `RELEASE_STATUS.json` records both install matrices as not passed, so no test
  count from them is quoted here. They predate Qualification v2 and do not
  satisfy its required second-machine scientific reproduction.
- The 1,095-word paper compiles with the official Open Journals Inara image and
  its four-page draft rendering has been visually checked.

## Required before `local_release_candidate`

- Adopt and execute every source-locked Qualification v2 case and pass all six
  Mark 1 release-blocking scopes. The current v2 state is
  `blocked_panel_incomplete`; no v2 scientific execution has occurred.
- Reproduce the successful scientific invariants on a second machine.
- Keep alpha-helical MembraneOrient visibly experimental and non-blocking until
  its own unchanged development and held-out gates pass.
- Pass the exact-mapped ABL-family StateAtlas held-out gate; other protein
  families remain prototype-only.
- Complete the third-party license and redistribution audit.
- Resolve or explicitly narrow the ActState generic catalytic-residue screen;
  membership in a broad residue set is not a position-specific chemistry test
  and cannot by itself establish that an annotated site is disrupted.
- Resolve the SF-CSA reciprocal-best-hit integration: the pipeline currently
  writes RBH status after structural classification, so
  `probable_same_function` is not reached by normal end-to-end execution.

## Required before `submission_eligible`

- Create the public repository only after explicit approval.
- Establish genuine public development history and public issue/review channels.
- Record successful independent installation and meaningful research use.
- Resolve contributor identities, affiliations, ORCIDs, and citation metadata.
- Recover exact versions for every AI tool used and finalize the required
  tool/scope/human-review disclosure.
- Approve conflict-of-interest and funding statements.
- Complete human scientific/editorial review of the compiled paper and replace
  journal-assigned draft metadata during submission.
- Make a release and archive only after explicit approval.

The current authoritative state is `pre_public_preparation`; passing local tests
must never be presented as JOSS acceptance or external scientific validation.

## External readiness review

An independent readiness audit against the JOSS review checklist is recorded in
[`docs/JOSS_READINESS_REVIEW_2026-08-26.md`](../docs/JOSS_READINESS_REVIEW_2026-08-26.md).
It confirms this file's blocking items and adds packaging, documentation, and
claim-sourcing findings not tracked here. Its verdict matches `PREPUBLIC_AUDIT.md`:
not submission-eligible.
