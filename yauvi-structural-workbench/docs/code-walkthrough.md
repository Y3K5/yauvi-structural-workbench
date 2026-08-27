# Code walkthrough

This guide explains the reviewer distribution from the outside in. It describes
what the code actually does and keeps each scientific claim within the evidence
that produced it.

## One analysis, end to end

1. `yauvi_structural_workbench.cli` accepts a local analysis request.
2. `yauvi_platform.structural_workbench.StructuralAnalysisStore` validates the
   analysis type, file roles, extensions, sizes, and SHA-256 identities.
3. Preflight checks required evidence, exact mappings, optional runtimes, and
   analysis-specific parameters. Missing required evidence blocks; it is not
   converted into a favorable score.
4. The store invokes one registered command-line module. Browser text cannot
   provide an arbitrary command, URL, or filesystem path.
5. The module writes raw JSON/TSV evidence. The store preserves those files,
   records the command and versions, and produces deterministic
   `REPORT_DATA.json`, `REPORT.html`, `RAW_EVIDENCE.zip`, `CHECKSUMS.json`, and
   `RUN_MANIFEST.json` files.
6. The report renders evidence and non-claims. Display rounding never changes
   the raw values.

## Reviewer package map

| Layer | Main code | Responsibility |
|---|---|---|
| Reviewer CLI | `structural-workbench/src/yauvi_structural_workbench/cli.py` | Local create/add/validate/run/export commands and loopback server launcher |
| Analysis store | `platform/src/yauvi_platform/structural_workbench/store.py` | Content-addressed inputs, preflight, registered execution, immutable runs, deterministic reports |
| Source assistance | `platform/src/yauvi_platform/structural_workbench/sources.py` and `sources/src/yauvi_sources/` | Registered public artifacts, acquisition policy, checksum cache, explicit adoption |
| Coordinate QC | `structqc/src/structqc/` | Coordinate identity, residue mapping, provenance declaration, PAE and imported validation |
| Membrane frame | `Membrane Orientor/memorient/src/memorient/` | Context-declared orientation geometry and residue zones |
| State comparison | `state-atlas/src/state_atlas/` | Sequence-mapped Kabsch comparisons, RMSD/RMSF, clustering, and bounded state resemblance |
| Site context | `site-context/src/site_context/` | Exact annotation-to-coordinate mapping, ligand/cofactor context, and separate pocket evidence |
| Activity-state screen | `activity-state/src/actstate/` | Separate completeness, geometry, occupancy, conformation, and assembly signals; one legacy label remains under review |
| Assembly context | `assembly-context/src/assembly_context/` | Contacts, interface residues, assembly evidence, method-specific SASA and burial |
| Structure/function comparison | `sf-csa/src/sf_csa/` | Frozen Foldseek and DIAMOND search universes, separate evidence legs, bounded interpretation vocabulary |

## What each workflow measures

### StructQC

The parser inventories models, chains, author and label residue identifiers,
observed residues, missing atoms, B factors or predicted confidence, and exact
file identity. When a reference FASTA is supplied, sequence alignment creates a
traceable reference-position map. A provenance JSON declaration controls whether
coordinates are treated as experimental, predicted, or unknown; filename
sniffing is only a warning. Imported wwPDB/MolProbity-like metrics remain
external validation evidence. None of this establishes native conformation or
function.

### MembraneOrient

The selected biological context chooses the applicable geometry: beta-barrel
normal, transmembrane-helix belt, anchor-relative frame, or soluble SASA-only
mode. Coordinates are rotated into the inferred frame and residues receive
geometric zone/facing labels. Random-rotation checks measure numerical stability.
The output is a modeled coordinate frame; it is not intact-cell exposure,
antibody access, topology-assay evidence, or an experimental membrane position.

### StateAtlas

Static models, multi-model coordinates, or optional MD trajectories are mapped
to experimental reference structures by sequence and aligned with the Kabsch
algorithm. Per-frame RMSD, RMSF, clustering, unresolved fractions, and distances
to two-sided active/inactive reference sets are retained. `active_like` and
`inactive_like` mean structural resemblance only; biochemical activity requires
independent evidence.

### SiteContext and ActState

SiteContext maps reference-sequence annotations onto exact coordinate residues
and keeps curated catalytic roles, observed components, geometry, and predicted
pockets separate. ActState combines five named signals without adding them into
a universal score. Its generic catalytic-residue membership screen is not
position-specific chemistry and is recorded as a release blocker; the legacy
`active_site_disrupted` label must not be interpreted as observed loss of
catalysis.

### AssemblyContext

The module consumes an explicit biological assembly, not a title-derived guess.
It reports inter-chain heavy-atom contacts, interface residues, stoichiometry
evidence, and solvent-accessible surface changes. An unapplied mmCIF assembly
operator blocks the run. Burial is coordinate- and method-specific; it does not
prove binding, obligate oligomerization, native exposure, or function.

### SF-CSA

Foldseek and DIAMOND are invoked as separate external runtimes against
checksum-pinned structure and proteome universes. Structural and sequence rows
remain separate. A closed vocabulary distinguishes exact self-control, bounded
mechanism-class evidence, structural analogy, candidate divergence, and
unresolved evidence. Reciprocal-best-hit status is currently computed too late
to reach its documented promotion gate during ordinary execution, so that gate
is a pre-public integration decision rather than a validated capability.

## Why the reports are reproducible

- Inputs are identified by SHA-256, not only filenames.
- Acquisition and adoption are separate actions.
- Completed run directories are immutable and content-derived.
- JSON keys, TSV ordering, ZIP entry order, and ZIP timestamps are normalized.
- Absolute workspace paths are scrubbed from scientific outputs.
- Missing evidence, runtime versions, parameters, evidence class, and non-claim
  text travel with the result.

## Intended uses

The suite can support coordinate intake review, experimental-versus-predicted
model comparison, membrane-frame hypothesis generation, conformational-reference
triage, catalytic-site mapping, assembly-interface inspection, and bounded
structure/sequence comparison. It is designed for research planning, teaching,
method development, and traceable analysis review. It is not a vaccine validator,
clinical tool, docking package, MD engine, wet-lab protocol, or substitute for
experimental structural and biochemical validation.
