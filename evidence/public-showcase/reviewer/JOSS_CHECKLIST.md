# JOSS preparation checklist

This checklist separates retained historical evidence, the changed local
candidate, work that requires the author, and JOSS eligibility gates. The
current draft and records do not authorize a public release or submission.

## Manuscript and retained evidence

- The paper includes Summary, Statement of need, State of the field and
  build-versus-contribute justification, Software design, Validation and
  research use, Research impact statement, Limitations, AI usage disclosure,
  Conflicts of interest and funding, Acknowledgements, and References.
- The current body is 1,631 whitespace-delimited words after excluding YAML
  front matter, the `# References` section, and image-caption lines; section
  headings are included. This is within the 750–1,750-word range; the count does
  not replace review of the final compiled manuscript.
- Qualification v2 collection 2.11 adopts 94 of 110 records across six panels:
  ABL 14, StructQC 16, functional-site context 16, assembly 16, SF-CSA 16, and
  beta-barrel membrane orientation 16. The alpha-helical membrane stratum is
  unadopted and non-blocking.
- Public workflow run [35295451246](https://github.com/Y3K5/yauvi-structural-workbench/actions/runs/35295451246)
  completed on 18 September 2026 for commit
  `7981148e70c5eaa6424e3608f09cd65bcfb35834`. Its seven CI runners passed all
  five release-blocking panels and covered macOS and Ubuntu on Python 3.10,
  3.11, and 3.12, plus Ubuntu arm64 on Python 3.12. See the retained, sanitized
  [sanitized CI observation](../evidence/preparation-2026-09-20/CI_OBSERVATION.json).
  This is evidence for that public commit, not the changed local
  candidate.
- The retained aggregate `EXECUTION_SUMMARY.json` reports 67/110 across five
  executed panels and predates the later SF-CSA result. It is a historical
  execution total, not the adopted-record count or a replacement for the later
  CI observation.
- The current local offline suite reports 684 passed, 15 skipped, 5 deselected,
  and zero failures or errors on macOS arm64 / Python 3.12.0. Its machine-readable
  record is [`OFFLINE_TESTS.json`](../evidence/preparation-2026-09-20/OFFLINE_TESTS.json).
  A separate case-safety/input-validation suite reports 5 passed in
  [`PUBLIC_CASE_SAFETY_TESTS.json`](../evidence/preparation-2026-09-20/PUBLIC_CASE_SAFETY_TESTS.json);
  keep this count separate from the full offline suite.
- The current wheel was freshly installed into an isolated environment. The
  installed distribution reproduced the public case twice with 12 output files
  byte-identical. The installed-wheel browser check exercised the synthetic
  StructQC flow with zero console errors. See
  [`CANDIDATE_VERIFICATION.json`](../evidence/preparation-2026-09-20/CANDIDATE_VERIFICATION.json)
  and the [research-use handoff](docs/JOSS_RESEARCH_HANDOFF.md). These checks
  establish software behavior and reproducibility, not scientific qualification
  or whether the author considers the case qualifying research use.
- The manuscript describes the public workflow results as historical and makes
  no claim that they qualify the changed local candidate. The current 1,631-word
  draft compiled successfully with Inara (exit 0), and all six pages were
  visually inspected. The review rendering is
  [paper-review-2026-09-20.pdf](../build/joss-preparation-2026-09-20/paper-review-2026-09-20.pdf);
  its build record is
  [PAPER_BUILD.json](../evidence/preparation-2026-09-20/PAPER_BUILD.json).
  Draft placeholders remain. The author still needs to review scientific
  claims, citations, disclosures, metadata, and final wording; rebuild the PDF
  after any resulting changes.

## Local preparation and author handoff

- Review the completed public-data structural-portfolio case study in the
  [research-use handoff](docs/JOSS_RESEARCH_HANDOFF.md) and decide whether it
  documents qualifying developer research use. The report, locked protocol,
  sources, outputs, and limitations are recorded under `../examples/structural-portfolio/`.
  This is not independent-human use or biological validation.
- The current wheel's clean installation and repeated public-case reproduction
  are recorded above. If the artifact selected for any later release differs,
  repeat the installation and case checks against that exact artifact, recording
  its hash, inputs, environment, outputs, and discrepancies. This does not
  retroactively change the historical CI record.
- Have the author review scientific claims, citations, author metadata,
  acknowledgements, funding, conflicts, and the AI disclosure. The current
  manuscript explicitly records that human review is pending.
- Resolve remaining third-party attribution questions and ensure notices match
  the sources actually used. Keep any unresolved licensing decision visible.
- Review the successful draft rendering and build record above, complete the
  author reviews, then rebuild and inspect the exact manuscript and metadata
  selected for any later release.

## JOSS eligibility gates

- Demonstrate more than six months of active public development. The recorded
  public history began on 27 August 2026; the calendar threshold is after
  27 February 2027, conditional on continued substantive public activity.
- Document research use of the software. The developer's own documented use can
  satisfy this requirement; independent adoption is not a prerequisite. The
  public-data case study is complete, but the author's assessment and review of
  whether it meets this criterion remain pending. The exact software revision
  and research-use record must be traceable.
- Confirm journal-facing author, affiliation, ORCID, citation, funding, conflict,
  acknowledgements, and AI-use statements, and complete human review.
- Select and identify the exact release and archive record only after the
  author approves the exact outgoing artifact. Public release, archival deposit,
  and JOSS submission remain separate approval boundaries.

Software tests, benchmark execution, qualification, research use, journal
eligibility, and author approval are distinct claims. Passing one does not
establish another.
