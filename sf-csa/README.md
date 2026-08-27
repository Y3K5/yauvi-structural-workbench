# sf-csa

Structure–Function Comparative Species Analysis.

Compares a checksum-pinned protein model against a frozen experimental structure
database and a local proteome universe, and keeps the two answers **apart**. A
Foldseek title can nominate an architecture; it can never, by itself, become a
function claim.

```bash
sf-csa build-manifests --spec campaign.json --out config/
sf-csa run      --queries config/target_manifest.json \
                --databases config/database_manifest.json --output results/
sf-csa verify   --output results/ --databases config/database_manifest.json
sf-csa describe                 # machine-readable IO contract
sf-csa validate --queries ... --databases ...
sf-csa fetch --plan             # what raw files are needed, and where from
```

## The two legs stay separate

Structural similarity (Foldseek, against a frozen PDB snapshot) and sequence
homology (DIAMOND, against a local proteome set) are reported as **separate
outputs** and are never merged into one similarity score. They fail differently:
a fold can be conserved where sequence is not, and a sequence match can span a
domain that carries no function.

## The interpretation vocabulary is closed

Six labels, and nothing else is ever emitted:

| label | what it means |
|---|---|
| `exact_function_supported` | exact accession identity — a control, not an independent experiment |
| `probable_same_function` | reciprocal best hit plus whole-architecture match plus compatible mechanism |
| `same_mechanism_class` | compatible architecture and independently named mechanism |
| `structural_analogy_only` | fold similarity with no concordant function evidence |
| `candidate_functional_divergence` | shared framework, differing specificity |
| `unresolved_or_conflicted` | threshold not met, or the evidence is contested |

A hit that fits none of them is `unresolved_or_conflicted` — not a new category.

**Pre-public implementation note:** reciprocal-best-hit status is computed in
the sequence leg after structural rows have already been classified. Normal
end-to-end execution therefore does not currently feed RBH evidence into the
`probable_same_function` promotion gate. The result fails conservatively toward
`same_mechanism_class`, but the label/documentation contract must be reconciled
before release.

## The campaign is data, not code

This is the change that made the module independently runnable. Target
accessions, their organisms, strains, structures and mechanism groups used to be
module-level dictionaries in `build_manifests.py`, which is why the descriptor
recorded as a limitation that *"organism specificity still lives in
build_manifests.py and must be rewritten per organism."*

Now a **campaign spec** names them and a generic builder consumes it:

```json
{
  "root": ".",
  "sequence_manifest": "results/SEQUENCE_MANIFEST.tsv",
  "decision_ledger":   "results/SELECTION_LEDGER.tsv",
  "targets": [
    {
      "accession": "EX0001",
      "common_name": "ExampleProtein",
      "organism": "Examplus fictus", "strain": "T1",
      "uniprot_accession": "EX0001",
      "mechanism_group": "example_family",
      "protein_specific_boundary": "architecture does not establish substrate.",
      "structure_path": "results/EX0001.pdb",
      "source_proteome_path": "proteomes/example.faa"
    }
  ],
  "database": { "pdb_database": "db/structdb", "proteome_globs": ["proteomes/*.faa"] }
}
```

`build-manifests` resolves every path, verifies each sequence against the ledger,
checksums the structures and the database, and refuses on drift. **The package
ships no target dictionaries** — a test asserts it.

## Four tables that are configuration, not algorithm

| manifest key | what it decides |
|---|---|
| `mechanism_families` | regex patterns mapping a PDB title to a mechanism group |
| `contested_groups` | groups that may never be promoted to a shared-function claim |
| `divergence_sets` | groups sharing a framework but not a substrate |
| `title_traps` | titles that must never carry a hit up to a function claim |

Their built-in defaults are **periodontal-pathogen biology**. A campaign against
other organisms that does not override them is using the wrong table, so
`sf-csa validate` reports when the defaults will be used rather than letting them
be inherited silently. Each release records which tables it ran with, so two
releases built against different tables are visibly not comparable.

## Fail-closed

- No release present → the platform module **blocks**. It never launches a
  PDB-wide sweep implicitly.
- A query with no local structure → `candidate_missing_structure`, never a
  structural negative.
- Sequence-ledger drift, a missing structure, an unversioned database → the build
  stops. An unversioned database cannot be cited, so it may not be used.
- `verify` audits a release against the manifest it was **configured** with, not
  against the shape the release recorded for itself.

## Dependencies

None. `core.py` is standard-library only. Foldseek and DIAMOND are external
binaries resolved from PATH and version-pinned in the database manifest — they
are runtimes, not Python dependencies, and preflight is fail-closed on their
absence.

    pip install sf-csa
    pip install 'sf-csa[sources]'   # adds yauvi-fetch integration

## Tests

    pytest                              # offline, no external binary needed
    pytest -m "adapter"                 # adds the foldseek/diamond checks

## Inside the platform

Registered as `catalogs/modules/sf_csa.yaml`, adapter at
`platform/src/yauvi_platform/modules/native/sf_csa.py`. That descriptor's output
paths are project-scoped on purpose: it records how the module is wired into one
campaign. The module's independence lives here, in the package and its CLI.
