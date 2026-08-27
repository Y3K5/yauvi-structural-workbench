---
title: "YAUVI Structural Biology Platform — Mark 1: evidence-bounded, reproducible structural protein analysis"
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
date: 25 August 2026
bibliography: paper.bib
---

# Summary

YAUVI Structural Biology Platform — Mark 1 is a local Python software suite and browser
interface for six common structural-protein questions: coordinate identity and
quality, membrane orientation, conformational resemblance, functional-site
context, biological assembly interfaces, and comparative structural/sequence
relationships. Each analysis remains an independently installable command-line
package. The workbench creates typed configurations, verifies files by SHA-256,
invokes the registered command, and renders its JSON and tabular evidence into a
readable, printable report.

The central design constraint is that evidence dimensions are not collapsed into
a universal protein or function score. Unknown provenance remains unknown;
active-like describes coordinate resemblance rather than biochemical activity;
an annotated catalytic site is not observed catalysis; and a shared fold is not
exact functional transfer. Missing references, mappings, runtimes, or validation
reports remain visible as missing or scientifically incomplete evidence.

# Statement of need

Structural bioinformatics typically requires moving between coordinate
validation, residue numbering systems, sequence annotations, reference
conformations, biological assemblies, surface areas, and similarity searches.
Author and label residue identities in PDBx/mmCIF are particularly important
when a measurement must trace back to an exact sequence position. Community
validation systems such as wwPDB and MolProbity provide essential model-quality
evidence [@wwpdb; @molprobity], while Gemmi provides robust macromolecular format
support [@gemmi]. However, these results are often reviewed in separate tools,
and missing evidence can be mistaken for a favorable or negative result.

YAUVI provides one reproducible case boundary while leaving scientific
calculations in standalone packages. It is intended for structural biologists,
bioinformaticians, and trainees who need inspectable evidence and explicit claim
limits rather than an opaque ranking. The same deterministic reports can support
method development, target characterization, comparison of experimental and
predicted structures, and planning of follow-up validation without generating a
wet-lab protocol.

# State of the field and build-versus-contribute justification

Existing programs already provide stronger specialist calculations than a new
implementation should attempt to replace. Gemmi supplies standards-aware
coordinate parsing [@gemmi]; wwPDB and MolProbity supply community model
validation [@wwpdb; @molprobity]; OPM/PPM supplies a richer membrane-positioning
reference system [@opm]; and Foldseek and DIAMOND supply established structure
and sequence searches [@foldseek; @diamond]. YAUVI therefore invokes or imports
these methods through named adapters and preserves their scores, parameters,
units, versions, and limitations rather than reimplementing them.

The distinct contribution is the evidence boundary joining these methods: exact
input identity, residue mapping, explicit provenance, fail-closed preflight,
separation of scientific dimensions, deterministic reports, and claim ceilings
shared across six independently runnable workflows. Extending a single existing
specialist package would not provide this cross-method case and reporting
contract, while replacing the specialist algorithms would reduce scientific
quality and maintenance sustainability.

# Software design

The workbench stores input bytes content-addressably and creates immutable run
records. Task definitions declare the scientific question, accepted artifact
types, format validator, source assistance, missing-evidence consequence,
outputs, and claim ceiling. A source finder links input roles to official RCSB
PDB, wwPDB validation, AlphaFold DB, UniProt, SIFTS, M-CSA, PDB CCD, ChEBI, and
OPM/PPM resources. Network access is disabled by default. When explicitly
enabled, only registered artifact types and public identifiers can be acquired;
cache acquisition and adoption into an analysis are separate operations.

StructQC establishes coordinate provenance and residue identity before composed
workflows. MembraneOrient separates a beta-barrel path from an experimental
alpha-helical helix-axis path, with OPM/PPM retained as an external comparison
standard [@opm]. StateAtlas uses exact declared residue equivalences for the
candidate ABL-family Mark 1 scope, Kabsch alignment, RMSD/RMSF, deterministic
clustering, and two-sided experimental references. SiteContext and ActState keep annotation,
observed chemistry, and geometric competence separate. AssemblyContext reports
heavy-atom contacts, stoichiometry evidence, and method-specific solvent
accessible surface area. SF-CSA executes Foldseek [@foldseek] and DIAMOND
[@diamond] as separate structural and sequence legs against checksum-pinned
reference universes.

Completed or scientifically incomplete runs emit `REPORT_DATA.json`, printable
`REPORT.html`, a deterministic `RAW_EVIDENCE.zip`, checksums, and the canonical
run manifest. Display rounding never modifies the underlying scientific values.

# Validation and research use

The suite uses synthetic offline fixtures, schema tests, fail-closed boundary
tests, transformation and ordering invariance tests, deterministic output
comparisons, controller security tests, and standalone package tests. External
scientific qualification is reported separately from software correctness. The
first checksum-locked public collection includes wwPDB validation, OPM membrane
strata, KinCore-labeled two-sided conformational references, an M-CSA enzyme
case, a deposited assembly evaluated with FreeSASA, and a CATH-labeled SF-CSA
mini-database searched by Foldseek and DIAMOND. Four public cases pass their
predeclared gates and two remain partial. Qualification v2 separately freezes
scope-specific strata, development and held-out splits, evidence requirements,
and unchanged gates. Its public cases have not yet been adopted or executed. No
release qualification is claimed until all six Mark 1 release-blocking scopes
pass v2 and reproduce independently.

# Research impact statement

At this pre-public stage, the project has not recorded independent adoption,
published research use, or six completed external benchmark gates. The software
is therefore not presented as JOSS submission-eligible. Current reproducible
benchmark records expose both successful cases and scientific limitations; they
do not substitute for documented use in real structural-biology analyses,
independent installation feedback, and public issue-driven refinement.
Aspirational utility is not counted as realized research impact in the release
state.

# Limitations

The workbench does not replace experimental validation, density inspection,
biochemical assays, native-surface measurements, or expert curation of
conformational references. Optional external executables and databases retain
their own licenses and must be acquired independently. Predicted models cannot
establish activity, assemblies may be context-dependent, solvent-accessible area
is method-specific, and sequence or structure similarity alone cannot transfer
mechanism.

# AI usage disclosure

OpenAI Codex using GPT-5-family coding models assisted with code generation,
refactoring, interface copy, documentation, test scaffolding, and manuscript
drafting during private development. Anthropic Claude also assisted during
earlier private development, but its exact model/version record has not yet been
recovered; this is a release blocker that must be resolved before submission.
The human author reviewed and edited the assisted outputs, ran the recorded
software checks, selected the scientific boundaries and sources, and made the
core design decisions. The human author remains responsible for originality,
accuracy, licensing, ethical and legal compliance, and every manuscript claim.
AI-generated suggestions are never treated as scientific evidence or benchmark
results. AI tools will not be used for author-editor or author-reviewer
conversations except where a journal policy explicitly permits translation.

# Conflicts of interest and funding

Conflict-of-interest and funding statements have not yet been finalized for
submission. They must be supplied and approved by every listed author before the
paper can leave pre-public preparation.

# Acknowledgements

The project depends on the maintainers and data curators of its open scientific
dependencies and public reference resources. Those projects must be cited
independently when their methods or data are used.

# References
