# Research-use and human-review handoff

This handoff records one completed local public-data worked analysis and the
remaining author decisions. No participant pilot or independent-human review has
occurred. Keep all work local until the author reviews the exact files and
approves a specific release destination.

## Selected internal research-use target

**Public-data structural-portfolio case study** — the local worked analysis is
complete. It compares chain A in public RCSB PDB entries 2POW and 4RN4 against
the UniProt P00918 reference sequence using StructQC. The [report](../../examples/structural-portfolio/PUBLIC_RESEARCH_REPORT.md),
[protocol](../../examples/structural-portfolio/PUBLIC_RESEARCH_PROTOCOL.md), and
[reproduction guide](../../examples/structural-portfolio/PUBLIC_RESEARCH_GUIDE.md)
record the inputs, selection, environment, outputs, and claim limits. The author
must review whether this actual public-data workflow constitutes documented
developer research use for JOSS. It is not independent-human use or biological
validation. A runnable example, synthetic demonstration, qualification panel,
or software test alone is not research use.

The case study's author must record:

- the question and why the selected workflow is relevant;
- exact software revision or built artifact hash, input accessions and hashes,
  source terms, configuration, commands, environment, and output hashes;
- observations, missing evidence, exclusions, failures, and unresolved results;
- what the results support, what they do not support, and what independent
  experimental or expert evidence would still be needed;
- whether any result is suitable to describe as research use, with the author
  making that judgment.

The exact built wheel was also installed in an isolated environment on macOS
arm64 / Python 3.12.0. The installed distribution ran the public case twice;
12 compared output files were byte-identical. A separate installed-wheel browser
check completed the synthetic StructQC example flow with zero console errors.
These checks support reproducibility and software behavior only. They do not
establish scientific qualification or resolve whether the case constitutes
research use. The recorded state and wheel identity are in
[`CANDIDATE_VERIFICATION.json`](../../evidence/preparation-2026-09-20/CANDIDATE_VERIFICATION.json),
with case-specific environment and outputs in
[`PUBLIC_CASE_RUN_ENVIRONMENT.json`](../../evidence/preparation-2026-09-20/PUBLIC_CASE_RUN_ENVIRONMENT.json)
and [`PUBLIC_CASE_CASE_RESULTS.json`](../../evidence/preparation-2026-09-20/PUBLIC_CASE_CASE_RESULTS.json).

Do not use private or confidential research inputs for this public-bound case.
Do not invent a result or fill a gap with a favorable default. Preserve raw
outputs and failed runs; explain them in a separate interpretation record.

## Human tasks and owners

| Task | Owner | Evidence to retain |
|---|---|---|
| Review and interpret the completed public-data case study | Author | Confirm whether the question and outputs reflect research use; review exact revision/input identities, environment, results, exclusions, and claim limits recorded in the linked report |
| Review scientific accuracy, citations, figures, manuscript length, authorship and affiliation metadata | Author | Dated review record tied to the exact paper and release candidate |
| Confirm funding and conflicts statements and review the AI disclosure | Author | Explicitly approved wording; disclosure scope and known gaps stated plainly |
| Reproduce the exact candidate in a clean environment | Author or named independent technical reviewer | Artifact/input hashes, OS/Python/CPU, dependencies, commands, full outcomes, failures, and discrepancies |
| Decide whether to run a usability pilot and select participants | Author | Approved protocol and participant selection; do not infer consent or recruit without the author's decision |
| Review any outgoing release set for privacy, licensing, provenance, and exact contents | Author | Resolved file manifest, findings, excluded/held files, and approval for that exact destination |

The existing [pilot protocol](PILOT_PROTOCOL.md) proposes three computational
researchers and three experimental biologists, with a public or synthetic input
bundle. Its targets are at least five of six participants completing their task
without developer intervention, all six identifying the main interpretation
limit, and no uncorrected misunderstanding that converts a structural metric
into biological proof. These are proposed usability targets, not results or
JOSS eligibility criteria. Contact and participant choice remain the author's
responsibility. Do not require confidential inputs; keep participant identities,
consent records, and raw notes outside this repository.

## JOSS gates and supporting evidence

JOSS requires more than six months of active public development and documented
research use. The developer's own documented use can meet the research-use
requirement; independent adoption is not required. The completed public-data
case is available for the author's assessment, but that assessment is pending.
The recorded public history
began on 27 August 2026, so the six-month threshold is after 27 February 2027 and
depends on continued substantive activity. See the current [submission
requirements](https://joss.readthedocs.io/en/latest/submitting.html).

The public CI observation, qualification panels, installed-wheel reproduction,
and optional participant pilot answer different questions; none substitutes
for author review of the research-use record or establishes biological
validation. The recorded wheel is a local preparation artifact, not an approved
release candidate. Repeat verification against the exact artifact selected for
any later release if it differs from this wheel.

## Privacy and release boundary

Use public or synthetic data for materials intended to accompany the software.
Keep private research, participant identity, contact details, consent, and
unreleased outputs in their approved local locations. Before any later release,
the author must review the resolved outgoing file set, generated artifacts,
metadata, and evidence for private paths, identifying details, unsupported
claims, licensing conflicts, and secret-shaped strings. Surface findings for
the author's decision; do not silently remove uncertain material. Local
preparation is not authorization to publish, archive, submit, or contact anyone.
