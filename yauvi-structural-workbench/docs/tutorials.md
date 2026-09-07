# Six worked examples

These examples use synthetic coordinates and annotations. Their purpose is to
teach execution and interpretation limits; they are not research-use evidence
or external scientific qualification. All six were regenerated locally on
6 September 2026; their case-level evidence is recorded in
[the worked-example record](../implementation-evidence/2026-09-06/six-synthetic-examples.json).

## Installed offline start

Install the built distribution, then run these commands from an empty directory:

```sh
yauvi --workspace ./analysis-workspace example --analysis qc-example
yauvi --workspace ./analysis-workspace analysis validate --analysis qc-example
yauvi --workspace ./analysis-workspace analysis run --analysis qc-example
yauvi --workspace ./analysis-workspace analysis export --analysis qc-example --out ./qc-report
yauvi --workspace ./analysis-workspace workbench serve
```

In the browser, choose **Try an offline example → Load example**, inspect the required evidence,
run, load the structure and select a residue. Open the measurements and their
supporting records. The exported REPORT.html opens locally in a browser; the
archive preserves the underlying inputs and evidence. Do not interpret its
synthetic imported validation values as validation of a real protein.

For a missing-evidence example, select **Leave out validation evidence** before loading
it, or use `example --analysis qc-no-validation --without-validation`. Inspect
which dimensions remain unevaluated. A run may complete while scientific
evidence remains unavailable; those are separate states.

## Regenerate all six from the source package

From the package root in its installed environment, select a new output directory:

```sh
python tools/build_five_use_case_showcase.py --out /tmp/yauvi-worked-examples
```

The historically named builder also creates the sixth, SF-CSA ceiling example.
It refuses an existing output directory. The resulting SHOWCASE.html explains
the cases; SHOWCASE.json, inputs and runs retain measurements and command logs.
Each corresponding standalone command remains documented in
[the generated CLI reference](cli-reference.md).

## HUC-01: Can I safely interpret these coordinates?

Before mapping a variant or functional residue, are sequence identity, residue numbering, provenance, and validation bound to the exact model?

- Sequence coverage: 100.0%. Reference residues with mapped coordinates
- Sequence identity: 100.0%. Identity within the accepted mapping
- Validation: imported. Imported, not recomputed by StructQC
- Chain breaks: 0. Detected coordinate discontinuities

Interpretation limit: Native conformation, biological function, or experimental correctness beyond the imported validation record.

Exercise: identify the input and reference records, locate one residue or search hit, and explain what additional evidence would be needed for a biological claim.

## HUC-02: Which side of a membrane might a receptor expose?

Can a single-pass membrane protein be placed in a consistent coordinate frame and divided into membrane and sided residue sets?

- Structure label: tm_helix_experimental. Geometry route selected by the tool
- Orientation method: tm_helix_axis_v2. Named computational method
- Residues reviewed: 57. Coordinate-bound residue annotations
- Modeled surface set: 28. Geometry-derived candidates; not intact-cell evidence

Interpretation limit: Native topology, intact-cell exposure, antibody accessibility, expression, receptor function, or Mark 1 alpha-helical qualification.

Exercise: identify the input and reference records, locate one residue or search hit, and explain what additional evidence would be needed for a biological claim.

## HUC-03: Which experimental conformation does my model resemble?

When two bounded reference states are declared, which reference is geometrically closer and is the margin interpretable?

- Resemblance label: active_like. Bounded structural label
- Best RMSD: 0.000 Å. After declared C-alpha alignment
- Reference margin: 1.672 Å. Distance separation from the alternative reference
- Interpretable frames: 1/1. Unresolved frames would remain in the denominator

Interpretation limit: Biochemical activity, activation, inhibition, mechanism, efficacy, or a time-resolved transition pathway.

Exercise: identify the input and reference records, locate one residue or search hit, and explain what additional evidence would be needed for a biological claim.

## HUC-04: Are the declared functional residues structurally present?

Do curated residues map exactly, retain role-compatible identities, and overlap a separately identified pocket?

- Mapped sites: 3. Declared residues found in the structure
- Role compatible: 3. Identity fits the declared role vocabulary
- Maximum separation: 5.295 Å. Descriptive C-alpha geometry
- Pocket methods: fpocket. Scores stay specific to each named tool

Interpretation limit: Observed catalysis, ligand affinity, inhibition, druggability, physiological function, or clinical relevance.

Exercise: identify the input and reference records, locate one residue or search hit, and explain what additional evidence would be needed for a biological claim.

## HUC-05: Which residues become part of an oligomer interface?

In a declared two-chain assembly, which subject residues contact a partner and how much surface becomes buried?

- Assembly complete: true. Expected and observed chain inventory agree
- Contact residues: 3. Subject residues with a partner inside the cutoff
- Buried surface: 60.600 Å². Method-specific buried SASA
- SASA method: biopython_shrake_rupley_240_canonical_frame. FreeSASA is preferred when installed

Interpretation limit: Native oligomer abundance, binding affinity, intact-cell accessibility, physiological interaction, or mechanism.

Exercise: identify the input and reference records, locate one residue or search hit, and explain what additional evidence would be needed for a biological claim.

## Expert reproduction

Use the same accepted files and parameters in a new case revision. Compare raw
measurements and case-level decisions, not timestamps or host paths. The outer
CHECKSUMS.json binds the raw-evidence archive; the manifest inside that archive
binds its constituent files without a circular archive self-hash. An interrupted
or failed run is retained and cannot be reused as a completed result.

The browser worker runs one analysis at a time with a bounded queue. Each engine
invocation has a one-hour limit; cancellation terminates the owned process group.
Server restart marks unfinished jobs interrupted and requires deliberate retry.
