# Regions, review and replay

Coordinates and a complete reference record go in. The Workbench verifies their identities, maps exact sequence correspondence and source features, and measures represented component contacts. Source-bound JSON records, a residue table, a local 3D viewer and the usual report/evidence export come out.

This is an experimental local workflow. Static geometry does not establish affinity, a reaction history, a transition probability or biological activity. It adds no qualified scientific scope.

## Start without downloading anything

In the Workbench, choose **Try an offline example → Regions and component contacts**, load it, check readiness and run. In **Findings**, choose **Open region explorer**. The example is synthetic software data with a missing residue, two alternate conformers and an incomplete component; it carries no biological claim.

The equivalent commands are:

```bash
yauvi --workspace ./my-workspace example --analysis region-demo --regions
yauvi --workspace ./my-workspace analysis validate --analysis region-demo
yauvi --workspace ./my-workspace analysis run --analysis region-demo
yauvi --workspace ./my-workspace workbench open
```

The demonstration remains scientifically incomplete because it has no experimental validation or declared biological provenance. Its region measurements can still be inspected.

## Supply another case

Choose **Region explorer** when creating an analysis. Attach deposited mmCIF coordinates and a complete local UniProt entry JSON. Downloads require the existing explicit retrieval action; this workflow makes no automatic external requests.

```bash
yauvi --workspace ./my-workspace analysis create --analysis regions --type region_explorer --question 'Which recorded regions are represented?'
yauvi --workspace ./my-workspace analysis add --analysis regions --role structure --file ./coordinates.cif
yauvi --workspace ./my-workspace analysis add --analysis regions --role reference_record --file ./reference.json
yauvi --workspace ./my-workspace analysis validate --analysis regions
yauvi --workspace ./my-workspace analysis run --analysis regions
```

Model index defaults to zero, the first deposited model. Assembly defaults to `asu`, the deposited asymmetric unit. Set `model`, `assembly_id` and `conformer` through Parameters when needed. Conformer choices are `auto` (highest mean occupancy per residue), `blank`, or a deposited alternate ID combined with blank atoms. Missing atoms do not fall back to another alternate. Alternate correlations remain unknown.

Mapping uses unique exact contiguous correspondence. A deposited sequence difference may be excluded only when its accession, chain, position and represented component match its recorded declaration and the remaining sequence correspondence is unique. Exact segments are mapped around that position; the changed component is never normalized into an annotated residue. Other mismatches remain unmapped.

For verified explicit mappings or source-backed assertions, attach `region_declaration`. Its schema is:

```json
{
  "schema_version": "1.0",
  "binding": {"coordinate_sha256": "EXACT_COORDINATE_SHA256", "sequence_sha256": "EXACT_REFERENCE_SEQUENCE_SHA256"},
  "sequence_mappings": [{"entity_id": "1", "label_start": 1, "reference_start": 1, "length": 100}],
  "assertions": []
}
```

Mapping segments must match exactly; overlapping segments and unknown entities are rejected. Assertion records use the existing biological-case contract: `id`, `property`, `value`, evidence `state`, `basis`, a locked `source_id`, `scope` and `interpretation_limit`. Available source IDs are `reference`, `coordinates` and, when supplied, `declaration`. Source citations can be recorded within the locked declaration. A state label additionally requires exact structure, coordinate checksum, sequence checksum, model, assembly, `frame: deposited_coordinates` and `conformer_selection` scope. Preserve construct descriptions and source limits in the value. Occupancy never supplies a state label.

## Inspect the result

The viewer separates **mapped annotations**, **measured contacts** and **literature context**. Region cards show residue/atom selections, source references, missing coverage and exact binding. The binding includes the case manifest and complete locked source inventory as well as coordinate, reference, mapping, atom, method, model, assembly, frame and conformer identities. Select a region, zoom to it, inspect the component, or reset to the whole structure. Any represented atom can be selected for a labeled geometric distance. Keyboard focus and Enter activate region cards. Changing case, structure, model, assembly or conformer clears incompatible selections.

Contacts use AssemblyContext's existing **5 Å heavy-atom convention**, excluding hydrogens, deuterium and zero-occupancy atoms from measurement. A closest protein/component atom pair is retained per contacted residue. Distance is not a hydrogen-bond assignment, physical contact area or affinity estimate. Deposited component bonds and covalent links are retained where unambiguous; chemical names alone do not establish substrate/drug/adduct roles.

Chemical atom coverage uses deposited component atom lists, not density or independent chemistry validation. Unknown coverage remains unknown. Curated compartment annotations do not create membrane planes. Only a compatible recorded placement may appear; unknown sides remain unsigned. Legacy engine/placement evidence is withheld for explicitly selected alternate scenarios because those records lack the required conformer identity.

`?no_webgl=1` requests evidence-only viewing. Cards, sources, atom selection, distances and selection exports remain available. The bundled 3Dmol runtime and license are served locally; its existing bytes are unchanged.

Outputs include `REGION_VIEW.json`, `REGION_MANIFEST.json`, `REGIONS.tsv` and a self-contained `CASE/` directory, inside the normal run's raw evidence bundle. Original StructQC observations are retained separately for the compatible asymmetric-unit/automatic-conformer view. The exact reference JSON sequence is extracted to a separately checksum-recorded FASTA for StructQC; original coordinate chemical identities are unchanged. Missing provenance or validation does not become a clean QC result.

## Record a separate review

Use `yauvi evidence review --case-dir CASE --selection selection.json --record review.json --out reviews`.

`selection.json` contains `structure_id`, `model_id`, `assembly_id` and `conformer`. The review record contains the exact region/view `binding`, `target_id`, SHA-256 of the canonical original target, `source_ids`, `reviewer`, `rationale` and `adjudication`. Allowed adjudications are `supported_discrepancy`, `reference_applicability_limitation` and `unresolved`. An optional `target_pointer` selects an individual observation within a bound engine record using a JSON pointer; its checksum must match that specific observation. Canonical JSON uses sorted keys, compact separators, UTF-8, no non-finite numbers and one trailing newline.

Review receipts are content-addressed and append-only, stored outside the original case. They retain the original target and locked sources. This records a human adjudication; it does not clear, repair or rewrite the original finding. Do not attach another view's review to a new structure or conformer.

## Build and replay a portable region kit

```bash
yauvi evidence bundle --case-dir ./CASE --selection ./selection.json --out ./region-kit
yauvi evidence replay --bundle ./region-kit --receipt ./new-replay-receipt.json
```

The bundle explicitly copies the case manifest, its declared source files, expected view records, frozen method copies, dependency/method locks and instructions. Paths inside the bundle remain relative. Original case provenance is preserved byte-for-byte. Inspect the allowlist before moving a kit outside a private workspace; creating it does not authorize sharing it.

Replay requires the exact installed methods, Python major/minor, OS/architecture and scientific dependency versions recorded in the kit. A missing or mismatched runtime fails clearly; there is no automatic installation or network request. Bundled code is never executed. Region mapping, selected atom identities and 5 Å contact measurements are recalculated. Other engine outputs are **checksum verification only** and listed separately in the receipt. Default numerical comparison uses an absolute tolerance of `1e-6`, zero relative tolerance, and exact comparison for nonnumeric identities. Changed coordinates, models, copies, frames, conformers, mappings or references cannot reuse another view's selection.

Repeat in a fresh directory on the same machine to check portability. That is not independent researcher or second-machine reproduction. Recalculation of other engine methods, arbitrary replay scripts, docking, dynamics, energies and automatic pocket discovery remain outside this recipe.
