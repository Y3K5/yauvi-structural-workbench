# Pre-public JOSS and software audit

Audit date: 2026-08-25  
Release state: `pre_public_preparation`  
Publication authorization: **not granted**

## Audited scope

The JOSS candidate is the YAUVI Structural Workbench reviewer distribution, not
the private vaccine campaigns and not the complete YAUVI research platform. The
wheel contains StructQC, MembraneOrient, StateAtlas, SiteContext, ActState,
AssemblyContext, SF-CSA, the structural analysis store, the source registry, and
the reviewer CLI. Archived implementations, `subproteo`, `vaxpipe`, vaccine
campaign data, docking, molecular dynamics execution, and private sequence data
are outside this release.

## Fresh software evidence

- Current macOS Python 3.12 workspace and offline dependency set exercised.
- Current reviewer-scope offline result: **526 passed and 6 network/adapter
  tests deselected** on Python 3.12.7.
- Canonical module CLIs and SF-CSA deterministic fixture were included in the
  run; suites ran in separate processes to avoid legacy pytest basename
  collisions.
- The root wheel and all eight advertised standalone package/control-plane
  wheels built successfully offline with the project license included.
- The root wheel smoke-installed into a temporary virtual environment using the
  current local dependency set; all advertised package imports and a clean
  `yauvi analysis create` command passed. A dependency-fresh OS/Python matrix
  remains outstanding.
- JavaScript syntax checks passed for the loopback UI and public showcase.
- The pre-public JOSS-material validator passed for six workflow definitions.
- A historical privacy-minimized ZIP was installed and tested from extraction
  on macOS arm64 and in a clean Linux arm64 Python container. Both runs passed
  with zero failures. No `TEST_RESULT.json` was committed from either run, and
  `RELEASE_STATUS.json` records both install matrices as not passed, so the
  per-suite counts are not quoted here as evidence. The skips reflected optional
  MDAnalysis, the omitted private portal consumer, and omitted downloaded
  qualification-source files; they are not converted into passes. This predates
  Qualification v2 and is not the required second-machine reproduction of the
  current scientific invariants.
- The 1,095-word manuscript compiled successfully with the official Open
  Journals Inara container and all four rendered pages were visually reviewed.
  Draft DOI, volume, issue, page, editor, reviewer, and submission-date fields
  remain journal-assigned placeholders.

Software-test success is not external biological validation.

## Defects corrected during this audit

1. The reviewer CI tried to run the full private `yauvi_platform` test surface
   after installing a wheel that intentionally contains only
   `yauvi_platform.structural_workbench`. The CI is now scoped to the published
   structural test and the separate loopback-controller suite.
2. The recorded 655-test baseline mixed 175 unrelated private platform tests
   into the JOSS reviewer count. The clean reviewer/controller baseline is now
   previously recorded as 508 passed, 5 deselected, and 1 skipped. Qualification
   v2 guards bring the current local reviewer selection to 526 passed and 6
   deselected.
3. The build backend rejected the newer bare SPDX license form. The root and
   standalone packages now use backend-compatible Apache-2.0 metadata, the root
   wheel includes the complete license file, and offline wheel builds pass.
4. MembraneOrient documentation described OPM/PPM as experimental and implied
   antibody accessibility. The wording now identifies OPM/PPM as an external
   computational reference and outside-facing residue labels as geometric
   nominations only.
5. The deterministic ZIP writer flattened executable permissions, so the
   synthetic Foldseek and DIAMOND fixtures failed after extraction. POSIX file
   modes are now preserved and recorded in the bundle manifest.
6. One MembraneOrient palette test depended on a private portal stylesheet that
   is intentionally excluded from the reviewer package. The test now runs when
   that integration consumer exists and is explicitly skipped otherwise; the
   standalone scientific palette tests remain active.

## Scientific findings that remain open

1. **MembraneOrient:** Mark 1 is explicitly narrowed to beta barrels. The
   alpha-helical `tm_helix_axis_v2` path requires checksum-bound spans, reports
   unsupported sides as unresolved, and remains experimental and non-blocking.
   Its complete development and held-out panel is not yet adopted or executed.
2. **StateAtlas:** Reference Set v2 now requires exact UniProt ABL1 242-495
   mappings with at least 90 percent coverage and rejects non-default static
   selections. The ABL held-out Qualification v2 panel is not yet adopted or
   executed, so the scope remains prototype.
3. **ActState catalytic-residue screen: narrowed 2026-09-02.** The generic
   residue set is still not role-aware, and it no longer decides a label on its
   own. `active_site_disrupted` now requires a position-specific expected
   residue from an experimentally validated reference; without one a
   non-competent residue is reported as contradicting evidence and caps the
   label at `indeterminate`. The narrowing also closed the defect in the other
   direction: a substitution inside the competence set, which the generic screen
   could not see, is caught wherever an expectation is supplied.
4. **ActState occupancy:** a non-solvent heteroatom is detected, but identity is
   not yet proven against the declared cofactor in every case.
5. **SF-CSA RBH integration: fixed 2026-09-01.** Reciprocal-best-hit evidence
   is computed before structural classification and reaches the documented
   `probable_same_function` gate by measurement. `rbh` is an explicit
   keyword-only argument that a curator cannot set through a manifest field.
6. **SF-CSA parser consistency:** blank and non-numeric similarity fields do not
   yet share one explicit invalid-input policy.
7. Four passed v1 public qualification cases are still named cases, not broad
   stratified workflow-general accuracy evidence. Qualification v2 freezes the
   expanded panel and currently fails closed as `blocked_panel_incomplete`.

No MD, docking, affinity, binding-pose, immunogenicity, native exposure,
clinical, or experimental biological result was generated or inferred by this
audit.

## JOSS gates that remain blocked

- No Git repository or verifiable development history exists in this workspace.
- The required public repository and more than six months of genuine public,
  iterative development have not started and are not authorized here.
- Independent installation and meaningful research-use evidence are not yet
  recorded.
- The third-party license and redistribution audit is incomplete; downloaded
  benchmark source files must not be silently redistributed.
- Exact versions for every earlier AI assistant must be recovered or the
  disclosure narrowed with journal/editor guidance.
- Funding, conflict-of-interest, authorship, affiliation, and final citation
  metadata require human approval.
- The compiled paper still needs human scientific/editorial review and final
  journal metadata.
- Public release, archival deposit, repository creation, and submission require
  explicit authorization.

## Readiness judgement

The software is suitable for **private cross-platform testing and continued
hardening**. It is not a local release candidate and not JOSS submission-ready.
The most useful next evidence is an independent macOS/Linux install report,
followed by resolution of the ActState and SF-CSA interpretation contracts and
the two partial external scientific workflows.
