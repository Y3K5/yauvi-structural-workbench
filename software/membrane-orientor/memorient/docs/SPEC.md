# memorient — module contract specification

**Purpose.** `memorient` takes a predicted or experimental protein structure and decides,
down to the residue, how it sits in a membrane: which way is "out", where the bilayer is,
which residues face lipid / lumen / the extracellular world, and which are reachable by an
antibody. It is the geometry-and-accessibility stage of a reverse-vaccinology pipeline, but
it is standalone and context-aware — the same call works for a gram-negative outer-membrane
β-barrel, a eukaryotic single-pass receptor, or a soluble secreted protein, because the
biological assumptions live in a **context object**, not in the code paths.

This document is the module-by-module contract. It is written to be handed to Claude Code (or
any implementer) as the authoritative description of what each module must provide and
guarantee. Every signature below is the one the package actually ships (`memorient` v0.2.0).

---

## 0. Design invariants (read first)

These hold across every module and are the reason the design is shaped the way it is.

1. **Context, not conditionals.** Organism/membrane assumptions are carried by a
   `MembraneContext` (thickness prior, orientation method, which metrics are meaningful, LPS
   or not). Numeric modules read the context; they never branch on hard-coded organism names.
   A metric is computed **only if the context declares it** — a soluble protein never gets a
   membrane slab, a symmetric plasma membrane never gets an LPS-shielding score.

2. **The extracellular set is the invariant, not the coordinates.** Re-orienting a rotated
   copy of a structure will not reproduce the exact same rotation matrix (PCA in-plane axes
   are degenerate for a symmetric barrel). The guarantee memorient makes and tests is that the
   **set of residues labelled extracellular** is stable across arbitrary input rotations —
   measured by Jaccard overlap in `five_fold_validate` (≥0.95 on real OMPs).

3. **Geometry proposes, localization disposes.** The membrane fit says whether a surface is
   *geometrically* exposed. Whether it is *biologically* reachable (host-antibody-accessible)
   can be vetoed by an injected `LocalizationCall` (the P1 seam) — e.g. a protein that is
   periplasmic is not accessible however exposed its loops look. Default is pass-through.

4. **Pure NumPy/SciPy + BioPython.** No external membrane-orientation binaries. SASA is
   computed in-package (Shrake–Rupley), validated against BioPython (Pearson r=0.999).

---

## 1. `contexts.py` — the biological configuration registry

Stdlib-only (no numpy), so it can be imported anywhere without the compute stack.

### Enumerations
- `MembraneModel`: `ASYMMETRIC_LPS`, `SYMMETRIC_PHOSPHOLIPID`, `NONE`.
- `OrientationMethod`: `BARREL_NORMAL`, `TM_HELIX_BELT`, `ANCHOR_RELATIVE`, `SASA_ONLY`.
- `Metric`: `AROMATIC_GIRDLE`, `LIPID_PORE_GAP`, `HYDROPHOBIC_BELT`, `POSITIVE_INSIDE`,
  `LPS_SHIELDING`, `ROTATION_INVARIANCE`.

### Dataclasses
- `ThicknessPrior(mean: float, sd: float)` — frozen. `.penalty(d)` returns `((d-mean)/sd)²`,
  the quadratic penalty added (negatively) to the membrane-fit objective.
- `MembraneContext(name, description, membrane_model, orientation_method, thickness_prior,
  metrics: frozenset[Metric], has_membrane_sides: bool, lps_shielding: bool)`.
  Properties: `.has_bilayer` (membrane_model ≠ NONE), `.is_asymmetric` (LPS model).

### Functions
- `get_context(name: str) -> MembraneContext` — raises `KeyError` on unknown name.
- `list_contexts() -> tuple[MembraneContext, ...]`.
- `default_context() -> MembraneContext` → `gram_negative_om`.

### Registered contexts (the contract for what ships)
| name | model | method | thickness (Å) | notable metrics |
|---|---|---|---|---|
| `gram_negative_om` *(default)* | ASYMMETRIC_LPS | BARREL_NORMAL | 13.0 ± 2.0 | girdle, lipid/pore gap, belt, LPS-shielding |
| `eukaryotic_pm` | SYMMETRIC_PHOSPHOLIPID | TM_HELIX_BELT | 15.0 ± 2.0 | belt, girdle, positive-inside |
| `tm_receptor` | SYMMETRIC_PHOSPHOLIPID | TM_HELIX_BELT | 15.0 ± 2.0 | belt, girdle, positive-inside (**no** lipid/pore gap) |
| `gram_positive_surface` | NONE | ANCHOR_RELATIVE | — | SASA-based |
| `soluble_secreted` | NONE | SASA_ONLY | — | SASA-based |

---

## 2. `geometry.py` — structure IO + canonical framing

### `Structure` dataclass (parallel per-residue arrays)
`ca` (N×3 float), `resids` (N int), `resnames` (N str), `chains` (N str),
`sc_vec` (N×3, Cα→side-chain-centroid unit vectors, used by facing/girdle),
`plddt` (N float or None), `atoms` (list of per-residue atom dicts for SASA), `source` (str).
Methods: `.sequence` (property, one-letter str), `.copy()`, `.transformed(R)` (rotate all
coordinates + side-chain vectors by 3×3 `R`), `len(structure)` → residue count.

### Functions
- `rotation_matrix_to_z(v) -> 3×3` — Rodrigues rotation carrying unit vector `v` onto +Z;
  antiparallel case handled as `diag([1,-1,-1])`.
- `principal_axes(coords) -> (axes, eigvals, centroid)` — PCA.
- `canonical_rotation(coords) -> (R, centroid, info)` — deterministic PCA frame: axis **signs**
  fixed by third-moment skew, third axis = cross product (guarantees det +1), `info` reports
  `gap01`/`gap12` eigenvalue gaps and a `degenerate` flag when either gap < 0.02.
- `canonicalize(structure) -> (Structure, info)` — apply `canonical_rotation`.
- `load_structure(path, chain=None, model=0) -> Structure` — PDB/mmCIF via BioPython.
  **Hardened**: on a header-parser crash (malformed REMARK/BIOMT, as in OPM files) it falls
  back to parsing ATOM/HETATM lines only. `MSE`→`MET`. Non-standard residues skipped.
- `structure_from_string(text, fmt="pdb", chain=None)`, `to_pdb_string(structure)`,
  `write_pdb(structure, path)`.

**Contract note.** Raw-coordinate frame reproducibility is only guaranteed for structures with
non-degenerate PCA axes; the durable guarantee is at the extracellular-set level (§7).

---

## 3. `sasa.py` — solvent-accessible surface area

In-package Shrake–Rupley. Bondi radii (`_BONDI`), `PROBE_RADIUS = 1.40 Å`, Tien(2013) MaxASA
for RSA normalization.
- `atom_sasa(coords, radii, n_points=240, probe=1.4) -> per-atom SASA`.
- `compute_sasa(structure, n_points=240, probe=1.4, heavy_only=True) -> {"sasa", "rsa", "atom_sasa"}`
  (per-residue absolute SASA, RSA, and per-atom SASA arrays).
- `total_sasa(structure, **kw) -> float`.

**Contract note.** Fibonacci sphere points are generated in the lab frame, so per-residue SASA
of a *rotated* molecule carries sampling noise that shrinks as `n_points` rises; total SASA is
rotation-stable. Rotation-sensitive callers (validation) use `n_points ≥ 160`.
Validated vs BioPython on crambin (1CRN): total within 0.1%, per-residue Pearson r = 0.999.

---

## 4. `barrel.py` — membrane slab fit + classifier

### `fit_membrane(structure, ctx, n_scan=80, polish=True) -> MembraneFit`
Finds the bilayer slab that maximizes the objective
```
J = 1.0·ΔKD(lipid − pore) + 0.7·hydrophobic_belt + 0.5·aromatic_girdle − 0.25·thickness_prior(d)
```
where ΔKD is the Kyte–Doolittle hydrophobicity contrast between lipid-facing and lumen-facing
embedded residues. Search: 80 Fibonacci-sampled candidate normals + 3 PCA axes + a
hydrophobic-moment seed; Nelder–Mead polish on the top 6; embedded set refit per candidate.
`MembraneFit` fields: `normal` (unit), `center`, `half_thickness`, `centroid`, `score`,
`components` (per-term breakdown), `embedded_mask`, `n_embedded`, `inner_frac`, `delta_kd`.

### `classify_membrane_protein(structure, ctx, fit=None, dkd_barrel=1.0, dkd_surface=0.6, inner_frac_max=0.25, n_embedded_min=40) -> MembraneClass`
Returns `label ∈ {barrel, surface, soluble, tm_helix}`, `confidence`, `reasons`, `fit`.
Thresholds are arguments, defaults shown. Barrel requires strong ΔKD, enough embedded
residues, and a low interior fraction (a barrel is hollow; a globular blob is not).

---

## 5. `membrane.py` — zones, facing, accessibility, context metrics

Constants: zones `extracellular` / `extracellular_interface` / `hydrophobic_core` /
`periplasmic_interface` / `periplasmic`; accessibility `antibody_accessible` / `lps_shielded` /
`lipid_embedded` / `pore_lumen_facing` / `periplasmic` / `buried_interior`.
`INTERFACE_WIDTH = 4.0`, `LPS_BUFFER = 6.0`, `RSA_EXPOSED = 0.20`.

- `project_membrane(structure, fit, ctx, ec_sign, rsa=None) -> MembraneProjection` — per-residue
  signed membrane depth (`ec_depth`, positive = extracellular), `zone`, `facing`
  (lipid/pore via side-chain vector), `accessibility`. **Raises** if `ctx.has_bilayer` is False
  (soluble proteins have no membrane frame to project into).
  `MembraneProjection` fields: `ec_depth`, `zone`, `facing`, `accessibility`, `ec_sign`.
- `context_metrics(structure, fit, ctx, proj) -> dict[str, float]` — evaluates **only** the
  metrics named in `ctx.metrics`. Never returns a metric the context does not declare.

---

## 6. `labeler.py` — extracellular-side call + per-residue labels

### `call_extracellular_side(structure, fit, ctx) -> SideCall`
Decides which membrane side is extracellular by weighted vote. `SideCall` fields: `ec_sign`
(±1), `votes` (per-signal ±1/0), `scores`, `agreement`, `confidence`.
Signals and weights:
- **loop_architecture** (1.0): extracellular loops are longer; face score =
  `median(loop lengths) + 8.0·long_loop_fraction` — median, *not* mean, so one giant soluble
  arm (a POTRA domain, a TonB plug) cannot invert the call.
- **terminus** (0.7): classic OMP topology puts both N- and C-termini periplasmic.
- **positive_inside** (0.25, β-barrel tiebreaker): Arg/Lys enrichment on the inner side.

For α-helical contexts (`tm_helix_belt`) the weights invert to make **positive-inside the
primary signal** (weight 1.0), terminus a tiebreaker — the physically correct rule for
single-pass helices (glycophorin A's C-terminal RRLIKK → cytoplasmic).

### `label_residues(structure, proj, rsa, ctx, fit=None) -> LabelSet`
`LabelSet.labels` is a list of `ResidueLabel(resid, resname, chain, zone, facing,
accessibility, extracellular, rsa, ec_depth)`. `LabelSet.surface_set` is the set of resids that
are **antibody-accessible** (exposed + extracellular-facing + not LPS-shielded).
A residue's `extracellular` flag is **strict**: `zone == extracellular` AND
`rsa ≥ RSA_EXPOSED + 0.02` (the interface band does not count as extracellular proper).

---

## 7. `orientor.py` — unified entry point + validation

### `orient_structure(structure, context, localization=None, validate=True, n_points=240, n_validate_seeds=5) -> OrientationResult`
The one call most users need. Pipeline:
1. `canonicalize` the structure.
2. `compute_sasa`.
3. Route on `context.orientation_method`:
   - `BARREL_NORMAL` / `TM_HELIX_BELT` → `fit_membrane` + `classify` + `call_extracellular_side`
     → reframe so **+Z is extracellular, membrane core at origin** → `project_membrane` →
     `label_residues` → `context_metrics`.
   - `SASA_ONLY` / `ANCHOR_RELATIVE` → no membrane; `projection` is `None`, the surface set is
     the SASA-exposed residues.
4. Apply the `LocalizationCall` veto to `host_antibody_accessible`.
5. If `validate`, run `five_fold_validate` and attach the result.

`OrientationResult` fields: `context` (name str), `method`, `label`, `confidence`, `structure`
(oriented), `rsa`, `labels` (LabelSet), `host_antibody_accessible`, `localization`, `fit`,
`classification`, `side`, `projection`, `metrics`, `validation`.
Methods: `.summary() -> dict` (flat, JSON-friendly, includes `metric.*` keys), `.to_dict()`,
`.residue_table() -> list[dict]`, `.extracellular_resids() -> list[int]`,
`.write_pdb(path)`.

### `LocalizationCall(localization: str, surface_exposed: bool, source="default_passthrough", confidence=1.0)`
The P1 seam. Default pass-through trusts geometry; an injected call (e.g. from a subcellular
localization predictor) can force `host_antibody_accessible = False` regardless of geometry.

### `five_fold_validate(structure, context, localization=None, seeds=5, threshold=0.95, n_points=160) -> dict`
Orients `seeds` random rotations of the structure and measures the mean pairwise Jaccard
overlap of the extracellular residue sets against the reference orientation. Returns
`passed`, `mean_jaccard`, `jaccards`, `n_reference_extracellular`. This is the operational
definition of rotation-invariance (design invariant #2).

---

## 8. `viz.py` — 3D export

- `display_oriented(result) -> dict` — 3Dmol.js descriptor: oriented PDB text, per-residue
  accessibility colours, and a `membrane_slab` (leaflet z-bounds + LPS band for asymmetric
  contexts) that is **`None` when the context has no bilayer**. Fully JSON-serializable.
- `write_pymol_script(result, path)` — `.pml` colouring by accessibility, epitope surface as
  sticks, membrane leaflets as pseudo-atoms.
- `write_3dmol_html(result, path)` — self-contained HTML page (loads 3Dmol.js from CDN)
  rendering the oriented cartoon coloured by accessibility with translucent membrane slab.
  Colour palette shared across all three exporters (`ACC_COLORS`).

---

## 9. `cli.py` — command line

`python -m memorient <subcommand>`:
- `contexts [--json]` — list registered contexts.
- `describe <context>` — dump one context's configuration (exit 2 on unknown name).
- `orient <structure> [--context NAME] [--chain X] [--n-points N] [--no-validate]
  [--out-json F] [--out-pdb F] [--out-viz F] [--out-pymol F] [--out-html F] [--max-rows N]`
  — full pipeline; prints the per-residue table to stdout and the summary to stderr, writes
  any requested output files.

---

## 10. Correctness evidence (P4)

`examples/p4_benchmark.py` takes five gram-negative OMPs from OPM (coordinates oriented with
membrane normal = +Z), **un-orients** each with random rotations, re-fits with memorient, and
measures the angle between the recovered normal and OPM's known normal:

| PDB | protein | normal error | fitted d (Å) | OPM d (Å) | Jaccard |
|---|---|---|---|---|---|
| 1BXW | OmpA | 5.7° | 12.4 | 12.6 | 1.00 |
| 2POR | Porin | 0.1° | 13.5 | 11.7 | 1.00 |
| 1QD6 | OmpF-family | 8.4° | 13.0 | 12.0 | 0.97 |
| 2F1T | NspA | 9.3° | 11.6 | 12.0 | 1.00 |
| 1P4T | OmpLA/PldA | 13.7° | 12.1 | 12.5 | 1.00 |

**Mean normal angular error 7.4°** against an external experimental reference — correctness,
not just internal self-consistency. Fitted half-thickness tracks OPM within ~1.5 Å.

---

## 11. Test suite map

`pytest` (65 tests; mark `network` needs internet, RCSB/OPM):
`test_contexts` (10) · `test_geometry` (11) · `test_sasa` (7) · `test_barrel` (10) ·
`test_membrane` (6) · `test_labeler` (4) · `test_orientor` (8, incl. real-OMP strict
invariance) · `test_viz_cli` (9, incl. real-OMP CLI). Synthetic structure builders in
`tests/synthetic.py` (barrel / TM helix / soluble blob / ellipsoid + random rotation).
