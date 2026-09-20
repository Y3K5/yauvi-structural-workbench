# Six structural workflows

## Structure QC

Inventories coordinate identity, models, chains, residue numbering, missing
atoms, provenance, reference mapping, PAE, and imported validation. It does not
establish native conformation or function.

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
