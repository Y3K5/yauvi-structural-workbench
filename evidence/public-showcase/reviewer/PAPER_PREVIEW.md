---
title: "YAUVI Structural Workbench: evidence-bounded, reproducible structural protein analysis"
tags:
  - Python
  - structural bioinformatics
  - protein structure
  - reproducibility
  - provenance
authors:
  - name: Yuvraj Patel
    orcid: 0009-0002-2276-7336
    affiliation: 1
affiliations:
  - name: Independent Researcher, Ohio, United States
    index: 1
date: 6 September 2026
bibliography: paper.bib
---

# Summary

Researchers need several evidence types to interpret protein structures. YAUVI
Structural Workbench combines a local Python suite and browser interface for six
questions: coordinate identity and quality, membrane orientation, conformational
resemblance, functional-site context, assembly interfaces, and structural or
sequence relationships. Each analysis is an installable command-line package.
The workbench validates typed configurations and SHA-256 inputs, runs registered
commands, and presents JSON and tabular evidence in printable reports.

Evidence dimensions remain separate: unknown provenance stays unknown;
active-like resemblance is not biochemical activity; an annotated site is not
observed catalysis; and fold similarity does not transfer function. Missing
references, mappings, runtimes, and validation remain visible.

# Statement of need

Structural analyses span coordinate validation, residue numbering, sequence
annotations, conformations, assemblies, surface area, and similarity searches.
Author and label residue identities in PDBx/mmCIF are essential for tracing
measurements to sequence positions. wwPDB and MolProbity provide model-quality
evidence [@wwpdb; @molprobity], and Gemmi supports macromolecular formats
[@gemmi], but evidence often remains split across tools, with gaps easy to
misread.

YAUVI defines one reproducible case boundary while leaving calculations in
standalone packages. It serves structural biologists, bioinformaticians, and
trainees who need inspectable evidence and explicit claim limits. Deterministic
reports support method development, target characterization, and comparison of
experimental with predicted structures without generating wet-lab protocols.

# State of the field and build-versus-contribute justification

Existing tools provide stronger specialist calculations: Gemmi parses
coordinates [@gemmi], wwPDB and MolProbity validate models [@wwpdb; @molprobity],
OPM/PPM positions membranes [@opm], and Foldseek and DIAMOND search structures
and sequences [@foldseek; @diamond]. YAUVI invokes or imports them through named
adapters, preserving scores, parameters, units, versions, and limitations.

Biotite integrates sequence and structure representations, database access, and
external applications [@biotite]; ProDy supports structural dynamics and
ensembles [@prody]. Both are foundations and alternatives. YAUVI's comparison is
architectural; no speed or accuracy advantage has been benchmarked.

YAUVI's contribution is a case and reporting contract shared across six
independent workflows: exact input identity and residue mapping, provenance,
fail-closed preflight, separate evidence dimensions, deterministic reports, and
claim ceilings. The browser workbench makes these contracts accessible across
engines. This adds packaging and adapter maintenance; reusable numerical or
parsing improvements should still go upstream.

# Software design

![Shared execution and provenance architecture. Four status dimensions remain separate.](figures/architecture.png){#fig:architecture width=100%}


The workbench stores input bytes content-addressably and creates immutable run
records. Browser actions and command-line commands call the same case service,
which keeps the interface from becoming a second scientific implementation.
Parameter edits preserve earlier case revisions. A bounded local worker records
queued and running jobs; server restarts expose interrupted work for deliberate
resubmission. Failed attempts remain inspectable and cannot satisfy a completed-run
cache lookup. Reuse requires matching input, source-code, parameter, and runtime
identities, together with intact output checksums. Task definitions declare the scientific question, accepted artifact
types, format validator, source assistance, missing-evidence consequence,
outputs, and claim ceiling. A source finder links input roles to official RCSB
PDB, wwPDB validation, AlphaFold DB, UniProt, SIFTS, M-CSA, PDB CCD, ChEBI, and
OPM/PPM resources. Network access is disabled by default. When explicitly
enabled, only registered artifact types and public identifiers can be acquired;
cache acquisition and adoption into an analysis are separate operations.

A preparation step, `structprep`, precedes these analyses where crystallographic
coordinates carry crystallisation additives, alternate locations, or assembly
copies. It derives new coordinates rather than editing in place, records what it
removed and why, and refuses to remove a component that contacts the polymer
unless that component is named explicitly, so a silent change to interface
geometry is not possible. It declares that its output is not ready for force-field
work until a step that states a force field and pH has run. It answers no
scientific question itself and is not one of the six workflows.

![Six structural evidence questions. StructQC binds identity, provenance and a SHA-256 before any analysis interprets the coordinates; the five analyses consume that manifest independently and their evidence is never combined.](figures/six-evidence-questions.png){#fig:questions width=100%}

StructQC establishes coordinate provenance and residue identity before composed
workflows. Modified chemical components are preserved separately from explicit
parent-residue sequence normalization; normalization does not imply identical
chemistry. SIFTS provides a reference resource for sequence-to-structure mapping
[@sifts]. MembraneOrient separates experimental beta-barrel and
alpha-helical helix-axis paths, with OPM/PPM retained as an external comparison
standard [@opm]. StateAtlas uses exact declared residue equivalences for the
candidate ABL-family Mark 1 scope, Kabsch alignment, RMSD/RMSF, deterministic
clustering, and two-sided experimental references. SiteContext and ActState keep annotation,
observed chemistry, and geometric competence separate. Declared cofactors and
observed non-solvent groups cannot establish occupancy without exact component
identity and site proximity; that evidence remains unavailable where the current
ActState reader cannot establish it. M-CSA annotations support curated site
interpretation and retain their evidence boundaries [@mcsa]. AssemblyContext reports
heavy-atom contacts, stoichiometry evidence, and method-specific solvent
accessible surface area, including named FreeSASA calculations [@freesasa]. SF-CSA executes Foldseek [@foldseek] and DIAMOND
[@diamond] as separate structural and sequence legs against checksum-pinned
reference universes.

Completed or scientifically incomplete runs emit `REPORT_DATA.json`, printable
`REPORT.html`, a deterministic `RAW_EVIDENCE.zip`, checksums, and the canonical
run manifest. Display rounding never modifies the underlying scientific values.

# Validation and research use

The suite uses synthetic offline fixtures, schema tests, fail-closed boundary
tests, transformation and ordering invariance tests, deterministic output
comparisons, controller security tests, and standalone package tests. External
scientific qualification is separate from software correctness. Qualification
v2 freezes scope-specific strata, development and held-out splits, evidence
requirements, and numerical gates. Collection 2.11 adopts 94 of 110 records:
14 ABL state
records, 16 each for coordinate quality, functional-site context, assembly
interfaces, and SF-CSA, plus 16 beta-barrel membrane records. The alpha-helical
membrane stratum remains unadopted and non-blocking. Reference classifications
for ABL are informed by KinCore [@kincore]. A [public workflow run
35295451246](https://github.com/Y3K5/yauvi-structural-workbench/actions/runs/35295451246)
on 18 September 2026 (commit
`7981148e70c5eaa6424e3608f09cd65bcfb35834`) recorded passing results for all five
release-blocking panels: 14 ABL cases and one control; 16 coordinate-quality cases
and two controls; 16 functional-site cases and one control; 16 assembly cases;
and 16 SF-CSA cases. Its seven successful runners covered macOS and Ubuntu on
Python 3.10–3.12 and Ubuntu arm64 on Python 3.12. The retained aggregate execution
report records 67 of 110 cases across five executed panels and predates the
SF-CSA result; this is not the adopted-record count or a summary of that later
run. These are historical results for the cited public commit, not qualification
of the changed local candidate. Membrane orientation is research-only and
non-blocking; its beta-barrel accuracy gate failed, so no accuracy claim is made.

SF-CSA's no-false-promotion bound is labelled definitional, not scientific: its
analogy and unrelated strata are defined as different superfamily and the
classifier promotes only within a group, so the bound holds by construction.
Adoption certifies that panel's recall and evidence-separation gates only.

Hardening fixed three reproduction-checker defects: failing panels could be
aggregated as success, missing panels could discard runner summaries, and
blocking scopes could never receive a reproduced verdict. The checker now
validates panels independently and requires workflow coverage, source and
protocol identities, exact case IDs, controls, and per-case results. Missing or
incompatible evidence blocks qualification. Software tests, execution,
independent review, and author approval remain separate; no scope is claimed
qualified for the changed build.

![Historical qualification evidence. Counts describe curated panels and retain reference, experimental and unexecuted distinctions; they are not toolkit accuracy.](figures/qualification-history.png){#fig:qualification width=100%}

# Research impact statement

The public benchmark collection demonstrates workflows on curated structural
cases; it is not research impact by itself. A local public-data case study
compares two deposited human carbonic anhydrase II structures against the
UniProt reference sequence, with locked inputs, a predeclared protocol, and
reported limits. Its report and reproduction instructions accompany the
software. The author must review the case and decide whether it constitutes
documented developer research use for JOSS; it is not independent use or
biological validation. No private research inputs or results are included here.
JOSS also requires more than six months of active public development;
independent adoption is not a prerequisite. This draft and historical software
results do not establish submission eligibility.

# Limitations

The workbench does not replace experimental validation, density inspection,
biochemical assays, native-surface measurements, or expert curation of
conformational references. Optional external executables and databases retain
their own licenses and must be acquired independently. Predicted models cannot
establish activity, assemblies may be context-dependent, solvent-accessible area
is method-specific, and sequence or structure similarity alone cannot transfer
mechanism.

# AI usage disclosure

OpenAI Codex assisted with implementation, tests, interface text, documentation,
and manuscript drafting. Anthropic Claude assisted during earlier development.

Models were recovered from local session records rather than reconstructed from
memory, and every session below is one whose transcript references this
repository. OpenAI: `gpt-5.6-sol`, `gpt-6-astra`, `gpt-reserve`, `gpt-5.6-luna`
and `gpt-5.4-mini`, across nineteen sessions between 25 August and 7 September
2026. Anthropic: `claude-opus-5`, across twenty-one sessions between 29 August
and 7 September 2026.

Two limits are stated rather than smoothed over. First, this repository's public
history began on 27 August 2026 and development preceded it, so assistance
before that date is not captured by a repository-scoped search; records from the
surrounding workspace over the same period additionally show `claude-sonnet-5`
and `claude-opus-4-8` from 25 July 2026, which cannot be attributed to this work
specifically and are therefore not claimed as such. Second, session counts
measure which models were invoked, not how much any of them contributed. No
version has been inferred where a record was absent.

Earlier drafts record human review; review of the current changes remains
pending. Automated tests do not substitute for that review.
The human author remains responsible for originality, accuracy, licensing,
ethical and legal compliance, and all claims. AI output is not scientific evidence.
The author will handle editor/reviewer conversations without AI assistance,
except translation where journal policy permits it.

For this preparation on 20 September 2026, OpenAI Codex (`gpt-6-astra` and
`gpt-5.6-luna`) assisted with manuscript, checklist, roadmap, and research-handoff
updates. This entry covers this preparation only; it does not reconstruct model
use from 8–19 September.

# Conflicts of interest and funding

The author declares no competing interests. This work received no funding: no
grant, institutional, or commercial support was provided at any stage, and the
author conducted it as an independent researcher.

# Acknowledgements

The project depends on the maintainers and data curators of its open scientific
dependencies and public reference resources. Those projects must be cited
independently when their methods or data are used.

# References
