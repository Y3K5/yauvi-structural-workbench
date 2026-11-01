# Protein function records and unknown investigations

Local, hash-pinned protein annotations and supplied evidence go in. Workbench
shows a record for every catalogue member, separates evidence from proposals,
and saves versioned questions about the unknowns. No scientific engine runs.

Open **Protein function records & unknowns** in the analysis workspace sidebar,
or `/function-records/`. The release contains code and invented software checks.
It contains no proteome, study library, research configuration or biological
result. Without a supplied catalogue it reports that configuration is missing.

## Try the invented example

From this source checkout, with the existing project installed:

```sh
python examples/function_records_demo.py function-records-demo
export YAUVI_FUNCTION_RESEARCH_CONFIG="$(pwd)/function-records-demo/config.json"
yauvi --workspace function-records-demo/workspace workbench open --port 8968 --label "Invented function records"
```

The generator refuses to overwrite an existing directory. Its identifiers,
hash placeholders and descriptions are invented; they are not reference
sequences, predictions, observations or qualified labels. No optional evidence
configuration is supplied by the example. Clear any previously set
`YAUVI_PROTEIN_FUNCTION_EVIDENCE` before using an independent example.

The `examples/` helper is source-checkout material, not an installed CLI command.
The browser and Python capability are included in the wheel.

## Read a record

Filter by reference organism, original category, annotation/accession or evidence
gap. Each profile shows exact identity, function annotation or retained
hypothesis, separate location prediction/context, typed source relationships,
saved-form coverage, original holds and explicit unknowns. Expand assertions to
inspect their source hash, JSON pointer and interpretation limits. Export the
exact record as JSON.

An annotation, family label or related structure does not establish native
activity. A predicted location does not establish exposure or orientation.
Complex membership is distinct from pairwise binding. Missing partners,
structures and condition-specific evidence are unevaluated coverage, not
biological negatives. Original categories and scores are preserved.

## Investigate an unknown

Choose a gap and record a question, proposed explanation, alternative explanation
and evidence that could distinguish or refute them. Link exact source assertions,
or leave support explicitly absent. Seven evidence tracks cover identity,
domains, homology, location, partners, functional observations and conditions.

Save explicitly. Plans use the native `studies/…/revisions/` and content-addressed
object layout, with immutable revision hashes and exact source snapshots. Reopen
from the protein's investigation list, inspect the native record, and export a
readable dossier plus JSON and checksums. Earlier revisions can be retrieved by
their exact hashes using the native API. Stale revisions or changed source
snapshots refuse edits. Start a reviewed successor rather than replacing evidence.

Every hypothesis remains `user_proposed_not_measured`. This profile has no
measurement-label intake, general study import, scientific dispatch, model
training, sequence generation or experiment submission. Saving a plan does not
update a function annotation. Unsaved drafts survive navigation only in the tab;
save before closing or reloading.

## Supply a local catalogue

Set `YAUVI_FUNCTION_RESEARCH_CONFIG` to a resident, bounded JSON configuration:

```json
{
  "schema": "function-research-config-v1",
  "catalogue_path": "ABSOLUTE_LOCAL_CATALOGUE_PATH",
  "catalogue_sha256": "EXACT_64_CHARACTER_SHA256"
}
```

Paths are resolved locally and symlinked inputs are refused. The catalogue
contract is `functional-protein-catalog-v1`: exact dataset and source-ledger
hashes, organism/category definitions, unique protein identities with sequence
hashes and lengths, and explicit `location` and `structure` evidence states.
The invented generator and tests show the minimum shape. Optional profiles,
typed relationship records and references supply context; identities are never
guessed from a similar name.

Optional `YAUVI_PROTEIN_FUNCTION_EVIDENCE` points to a local hash-pinned review
configuration bound to the same dataset and catalogue. Supported compatibility
schemas cover annotation refinements and retained domain/reference reviews.
Every row must match exact key/sequence and declared accession, organism and
length. A missing or changed optional source stays visibly unavailable and
holds new plan writes. These readers ship without any source payload.

## Architecture and boundaries

`function_catalogue.py` validates local inputs. `function_records.py` derives
records and gap filters. `function_investigations.py` owns proposal validation
and revision rules. `function_record_store.py` provides bounded custody in the
existing native layout. The shared browser renderer uses same-origin guarded
routes from `function_record_routes.py`.

This public profile exposes function records and investigations. A fuller local
Research Home may inject its existing general StudyStore into the same core;
the broader designer, scientific execution adapters, case-specific workflows
and study libraries are not included by this update. No second platform,
remote storage service or external dependency is introduced.

Software tests cover exact identities, preserved missingness, stale refusals,
save/reopen/revision/export and guard behavior. They do not validate protein
function, biological interpretation, experiment feasibility or clinical use.
