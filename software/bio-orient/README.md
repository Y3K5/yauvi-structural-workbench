# Bio-Orient

Bio-Orient connects declared structural sides to their represented surroundings
and recorded evidence. It is an experimental interpretation layer in the local
YAUVI Structural Workbench. Computed geometry is not biological validation.

## Where it fits

| Engine | Contribution | Boundary |
|---|---|---|
| StructQC | Identity, observed coverage, imported validation, local component chirality | Chemistry does not assign compartments |
| StructPrep | Explicit parent/child derivation | Declare patches on the child bytes and run QC again |
| MemOrient | Input-frame plane and mapped compartment markers | Geometric or sequence heuristics remain predictions |
| AssemblyContext | Existing 5 Å contact convention | Solvent burial is not physical contact area |
| SiteContext | Scoped functional records and chemical roles | Nearby atoms do not establish a role |
| StateAtlas / ActState | Recorded state outputs and their original limits | Static structures do not establish dynamic probabilities |
| Workbench | Runs QC, composes bound inputs, presents graph and exports | Execution, qualification and biology remain separate |

## Run locally

Install the root workbench distribution for the complete experience:

```bash
python -m pip install ".[dev]"
```

For the standalone CLI, resolve its engine dependencies from this checkout in
the same installation command; these local engines need not be published on PyPI:

```bash
python -m pip install ./software/structqc "./software/Membrane Orientor/memorient" \
  ./software/assembly-context ./software/bio-orient
```

Its CLI exposes `describe`, `validate`, and `run`.

```bash
bio-orient describe
bio-orient describe --structure model.cif --model 0 --assembly-id asu > scope.json
bio-orient validate --structure model.cif --patch-declaration patches.json
bio-orient run --structure model.cif --patch-declaration patches.json --out results
```

`describe --structure` reports the exact binding, chemical component identities,
chain copies and a declaration starter. Choose the residues deliberately and
record the declaration's source. `--model` is a zero-based index; `model_id` in the
binding is the deposited identifier. `asu` explicitly selects the deposited
asymmetric unit. Other assembly IDs must have deposited recipes. No crystal
neighbors are added implicitly.

The workbench registers `bio_orient`:

```bash
yauvi --workspace local-workspace analysis create --analysis sides --type bio_orient \
  --question "What surrounds these declared surfaces?" --subject-id "Selected structure"
yauvi --workspace local-workspace analysis add --analysis sides --role structure --file model.cif
yauvi --workspace local-workspace analysis add --analysis sides --role patch_declaration --file patches.json
yauvi --workspace local-workspace analysis validate --analysis sides
yauvi --workspace local-workspace analysis run --analysis sides
```

Attach optional CCD files individually as `chemical_component`, their manifest as
`chemical_reference`, bound outputs as `supporting_evidence`, and membrane topology
as `topology_evidence`. Private coordinates and cases belong in separate local
workspaces, outside public-bound software.

## Declaration contract, version 1.0

Copy `binding` from `describe --structure`. Every field must match the selected
coordinates. Chain IDs are exact copies, such as `A:identity` for an mmCIF label
chain, or the labels Gemmi generated for a PDB file. Never guess a copy name from
the author-chain letter alone.

```json
{
  "schema_version": "1.0",
  "binding": {
    "coordinate_sha256": "COPY_FROM_DESCRIBE",
    "model_id": "1",
    "assembly_id": "asu",
    "frame": "deposited_coordinates",
    "conformer_policy": "highest_mean_occupancy_per_component; chemistry with alternates remains unevaluated"
  },
  "sides": [{
    "side_id": "declared-domain",
    "residue_set": [{"chain_id": "COPY_FROM_DESCRIBE", "auth_seq_id": 42, "insertion_code": ""}],
    "source": {"id": "local-declaration", "citation": "Record explaining this residue selection"},
    "basis": "curated",
    "equivalence_group": "declared-group",
    "biological_roles": [],
    "functional_annotation": []
  }]
}
```

Selections can represent a domain, site, membrane region or measured interface;
their residue sets remain explicit. Missing residues stay in the declared
denominator. Each side reports its centroid, elemental/component composition,
solvent-accessible area and area-weighted outward direction. Directions below
0.25 coherence remain unresolved. This threshold is an experimental display
criterion, not biological qualification. Sampling reuses MemOrient's 240-point
Shrake–Rupley method with a 1.4 Å solvent probe in an intrinsic proper frame.

Protein and nucleic-acid chains and represented nonpolymers become neighbor
objects. Known component classes are preserved using Gemmi's residue table.
Unknown categories stay generic. `object_annotations` can label an exact object
as a protein, nucleic acid, ligand, glycan, lipid, solvent or chemical component;
each annotation requires `object_id`, `kind`, and source `id`/`citation`.
`context_objects` explicitly declare compartments, cavities, membranes or solvent
contexts with an object ID and cited source; optional centroids must be finite.

For each represented neighbor, edges report minimum heavy-atom distance, contacts
at 5 Å, and the increase in sampled patch SASA when that neighbor is removed.
That increase measures occlusion in the represented assembly. It is not physical
contact area, a probe-specific reachability measurement, or an affinity score.
Angles to neighboring centroids and directed patch-to-patch normal angles report
their geometry explicitly. Chemical complementarity, curvature and reachability
remain unevaluated unless compatible scoped method output is supplied.
Neighbor removal preserves the full-assembly sampling frame and measures the
selected atoms only. Objects beyond the sum of probe-expanded radii contribute
zero burial. A supported membrane frame does not resolve a patch spanning its
core or both faces; that patch's orientational availability remains unresolved.

Side count, declared equivalence groups and distinct biological roles are separate
descriptors. Multisided, homosided and heterosided relationships can coexist.
Symmetry does not establish equivalent function. Missing evidence produces unknown
descriptors rather than negative biological findings.

## Supporting-output envelope

Supply JSON with `binding`, `engine`, `method`, `output`, and `output_sha256`.
Allowed engines are `memorient`, `assembly_context`, `site_context`, `state_atlas`,
`actstate`, and `patch_measurements`. Method records require an `id` and
`source_sha256`. Hash the output using UTF-8 JSON with sorted keys, separators
`,` and `:`, and finite numbers. The envelope file itself is also checksum recorded.

Side-specific `output.side_records` require a `side_id` and the exact declared
`residue_set`. Records may supply curvature, mobility, dynamic-state evidence or
functional annotations. State outputs remain stored with their original scope
and limitations; their recorded values are not inferred by Bio-Orient.
`output.edge_records` name `source_side`, `adjacent_object`, exact `residue_set`, a
source and `interpretation_limit`. Functional availability requires these scoped
records; proximity never supplies it. Conflicting edge records remain visible.
Existing condition-specific vitality records are not redefined.

## Membrane contexts

Optional `membrane_context` and an explicit author `membrane_chain` request existing
MemOrient placement on the selected source model. Alternatively, a compatible
MemOrient envelope supplies `input_coordinate_membrane` in the deposited frame.
Topology inputs require the coordinate SHA-256 and exact Bio-Orient `binding`.
New alpha-helical contexts require mapped transmembrane spans. Each span retains
exact chain/author residue/insertion-code mappings and its source.

| Context | First / second compartment |
|---|---|
| `gram_negative_im` | Periplasm / cytoplasm |
| `er_membrane` | ER lumen / cytosol |
| `mitochondrial_outer` | Cytosol / intermembrane space |
| `mitochondrial_inner` | Intermembrane space / matrix |
| `organelle_membrane` | Two explicitly supplied distinct names |

Compartment ordering is a serialization convention. It is not a biological
assignment of +Z. Topology `side_markers` contain `chain_id`, `auth_seq_id`, optional
`insertion_code`, `compartment`, and evidence `basis` (`curated`, `experimental`,
`predicted`, or `derived`). A cited source is mandatory. Markers must map uniquely
outside the modeled slab. Conflicts remain conflicting; missing/uninformative
markers remain unknown. Prediction evidence remains predicted. Signed exported
markers are re-mapped against selected coordinates before use. Expanded-copy
ambiguity requires explicitly mapped evidence rather than automatic transfer.
Organelle contexts do not inherit positive-inside heuristics.

## Local chemistry references

References are local CCD mmCIF files identified by exact component ID and SHA-256.
The [wwPDB CCD](https://www.wwpdb.org/data/ccd) defines component connectivity,
stereochemical assignments and ideal coordinates. Gemmi parses the local files;
Bio-Orient does not retrieve or submit components automatically.

```json
{"schema_version":"1.0","components":[
  {"component_id":"ALA","path":"ALA.cif","sha256":"SHA256_OF_THIS_LOCAL_FILE"}
]}
```

Standalone reference paths are relative to the manifest directory. Workbench
attachments are resolved by checksum and materialized into a derived manifest;
the original uploaded manifest and files are preserved. The reader uses
`_chem_comp_atom` atom IDs, elements, `pdbx_stereo_config`, ideal Cartesian columns,
and `_chem_comp_bond` connectivity. Atom-mapped signed volumes are compared to
reference geometry without assigning observed absolute R/S labels. Near-planar
centers, missing atoms/references, alternative conformers, ambiguous definitions
and unsupported chemistry remain unevaluated. D-residues and modified components
retain their chemical identifiers. Imported validator findings remain separate.
Biology-preserving transformations require finite orthonormal matrices with
determinant +1; reflections are rejected. Chemistry discrepancies are review
flags, with no coordinate repair or automatic downstream hold.

## Results and viewers

`BIO_ORIENT.json` contains the complete graph; `BIO_ORIENT_SIDES.tsv` and
`BIO_ORIENT_EDGES.tsv` provide tabular summaries. `BIO_ORIENT_LAYER.json` carries
the viewer data, and `RUN_MANIFEST.json` records inputs, methods, local reference
identities, software hashes, runtime versions and output hashes. The workbench
adds its QC output and includes these files in reports and evidence exports.

In Findings, select a patch, then a neighbor, or select a chemical component and
atom. Directions, compartment evidence, coverage, chemistry flags, availability
and source records remain inspectable. The connected biological viewer uses the
same graph renderer. Existing biological cases can add a `bio_orient` evidence
record using their existing locked source/run binding; the adapter verifies the
graph's atom and copy identities before exposing it. Existing cases require no
migration, and legacy unsigned placements remain unknown.

Automatic patch discovery, per-frame dynamics, fold handedness, turnover,
mechanistic strain and immunological scoring are later extensions.

See [delivery checks](VERIFICATION.md) and the [AQP1 tetramer check](KNOWN_COMPLEX_CHECK.md)
for measured verification scope.
