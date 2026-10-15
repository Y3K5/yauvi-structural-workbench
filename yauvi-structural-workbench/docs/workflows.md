# Structural workflows

Six workflows retain their recorded Mark 1 qualification scope. Bio-Orient,
DockPrep, and PoseEvidence are experimental additions. They do not extend the
frozen Mark 1 qualification panels.

## Structure QC

Inventories coordinate identity, models, chains, residue numbering, missing
atoms, provenance, reference mapping, PAE, and imported validation. It does not
establish native conformation or function.

Local chemistry checks compare atom-mapped chiral volumes to checksum-locked
CCD references. Missing references, atoms and ambiguous conformers remain
unevaluated; discrepancies produce review flags without changing coordinates.

## Membrane orientation

For beta barrels, places a coordinate model using the `barrel_normal` path and
labels residue positions. This path remains experimental under the current qualification manifest. The separate
alpha-helical `tm_helix_axis_v2` path requires checksum-bound transmembrane
spans and remains experimental; unsupported sidedness is unresolved. Neither
path proves intact-cell exposure. OPM/PPM remains an external qualification
reference.

## Conformational resemblance

Aligns query coordinates or trajectory frames to at least two independently
supported experimental references per state, reporting RMSD, RMSF, clustering,
and interpretable populations. For Reference Set v2, alignment and clustering
use only the exact declared domain while the resulting transform is applied to
the complete frame. The Mark 1 candidate scope is ABL-family resemblance over an exact
SIFTS or explicit UniProt ABL1 242-495 residue map with at least 90 percent
coverage. That scope remains prototype until its Qualification v2 held-out gate
passes. Labels describe resemblance, not activity.

The showcase separately compares the two locked ABL endpoints 2GQG and 2HYY
with a proper Cα rigid-body fit over 252 shared mapped residues. It displays
per-residue displacement, actual endpoint φ/ψ/χ1 and N–Cα–C geometry, and an
explicitly illustrative Cα interpolation. Missing/modified positions are not
filled or bridged. These derived measurements and their source hashes can be
downloaded; the interpolation is not a molecular-dynamics trajectory and the
original classification gate does not qualify its intermediate shapes.

## Functional-site evidence

Maps declared residues, roles, ligands, metals, and cofactors to exact
coordinates. UniProt annotations, M-CSA curation, observed chemistry, and pocket
predictions remain separate evidence panels.

## Assembly and interfaces

Measures contacts, interface residues, stoichiometry evidence, SASA, and burial
in a declared assembly. Incomplete assemblies produce lower bounds.

## SF-CSA comparison

Runs checksum-pinned structural and sequence search legs using Foldseek and
DIAMOND. Exact protein, homolog, shared fold, analogy, and unresolved evidence
remain distinct; no similarity result is automatically transferred as function.

The showcase illustrates two fold spaces: ferredoxin-like electron carriers
(P00198 / P00208) and OB-fold cold-shock proteins (P32081 / P0A9X9). It separates
frozen structural/sequence comparisons from independent experimental literature
supporting a shared functional role. Same-fold counterexamples show why a fold
assignment alone is insufficient. See the [guided onboarding examples](PILOT_PARTICIPANT_GUIDE.md#cross-species-comparison-two-guided-examples)
and [case references](SF_CSA_SHOWCASE_CASES.json).

## Bio-Orient (experimental)

Connects declared residue patches to measured neighboring objects and compatible
membrane, chemical, site and recorded-state evidence. Every finding is bound to
exact coordinates, model, assembly, chain copy and conformer. Distances, contacts,
sampled exposure and buried SASA are geometric observations. Structural,
orientational and functional availability stay separate; missing biological
evidence remains unknown. No dynamic probability is inferred from a static model.

Both viewers show patch/component selection, normals, labeled membrane planes,
compartments, chemistry flags and source records. Membrane assignments retain
supported, predicted, conflicting or unsigned states.

See [the input contract](../../software/bio-orient/README.md),
[software checks](../../software/bio-orient/VERIFICATION.md) and
[the public AQP1 example](../../software/bio-orient/KNOWN_COMPLEX_CHECK.md).
These do not add a qualified scope to Mark 1.

## Docking preparation (experimental)

`dockprep` builds checksum-bound local docking jobs for small molecules,
protein-protein pairs, and protein-peptide pairs. It validates exact residue
mapping, chemistry, source terms, conformational state, and human approval
before execution. AutoDock Vina, Meeko, and HADDOCK3 are external runtimes.
Runtime absence creates a named blocked evidence leg. It does not substitute
another engine.

Vina scores and HADDOCK scores retain their engine meaning and are never called
affinity or biological activity. A docking score is not affinity, engagement,
activity, efficacy, or mechanism.

The small-molecule Vina path runs from the installed `vina` binary and Meeko's
`mk_prepare_receptor` / `mk_prepare_ligand` commands. The macromolecular
HADDOCK3 leg requires a local runtime and the canonical runner
`refinement_workflow/run_haddock_matrix.py` (`YAUVI_HADDOCK_MATRIX_RUNNER`).
That runner and the pinned container image `protein-platform-haddock3:2026.7.0`
are not distributed with this repository. Anyone without them gets
`canonical_haddock3_matrix_runner` reported as missing evidence.

## Pose evidence (experimental)

`pose-evidence` imports local Vina, HADDOCK3, or read-only ClusPro poses and
reports method-appropriate clusters, contacts, clashes, geometry proxies,
restraint satisfaction, native-site overlap, state sensitivity, and repeated
seed consistency. Engine scores and units remain separate.

Docking scores remain engine-specific and are not affinities. Pose evidence is
not affinity, engagement, activity, or mechanism. Hydrogen bonds, salt bridges,
clashes, and metal coordination are geometry proxies. Imported ClusPro results
are read-only and are never uploaded by this module. Buried solvent-accessible
surface is left unevaluated unless an explicit SASA adapter is invoked; its
absence is reported and is not filled in.
