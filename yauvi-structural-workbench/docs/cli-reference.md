# CLI reference

Generated from each command's own `--help`, so it cannot drift from the code.
Regenerate with `python tools/build_cli_reference.py`.

Common conventions across the scientific modules:

- `describe` prints the module contract as JSON: inputs, outputs, and the claim ceiling.
- `validate` checks inputs without producing an evidence record.
- `run` performs the analysis and writes a result bundle to `--out`.
- `fetch` resolves registered public sources. Acquisition never adopts a file
  into an analysis automatically, and stays disabled without an explicit
  reference-fetch flag.
- Exit codes: `0` completed, `1` scientifically incomplete, `2` invalid input or
  configuration. `1` is a real result, not a failure to be retried around.


| Command | Purpose |
|---|---|
| [`yauvi`](#yauvi) | Structural analysis case store: create, add inputs, validate, run, export. |
| [`structqc`](#structqc) | Coordinate trust: completeness, provenance class, imported validation. |
| [`memorient`](#memorient) | Membrane orientation and sidedness labelling. |
| [`state-atlas`](#state-atlas) | Conformational-state resemblance against declared references. |
| [`site-context`](#site-context) | Functional-site roles, cofactors, ligands, and pockets. |
| [`actstate`](#actstate) | Activity-state evidence assembly. |
| [`assembly-context`](#assembly-context) | Biological assembly, stoichiometry, contacts, and burial. |
| [`sf-csa`](#sf-csa) | Structure- and sequence-based functional comparison. |
| [`yauvi-fetch`](#yauvi-fetch) | Registered public-source acquisition and staging. |


---

## yauvi

Structural analysis case store: create, add inputs, validate, run, export.

```
usage: yauvi [-h] [--workspace WORKSPACE] {analysis,workbench} ...

Local, evidence-bounded structural protein analysis.

positional arguments:
  {analysis,workbench}
    analysis            Create, validate, run, and export structural analyses.
    workbench           Serve the loopback-only browser workbench.

options:
  -h, --help            show this help message and exit
  --workspace WORKSPACE
                        Structural Workbench repository root.
```

### `yauvi analysis create`

```
usage: yauvi analysis create [-h] --analysis ANALYSIS --type TYPE --question
                             QUESTION [--subject-id SUBJECT_ID]

options:
  -h, --help            show this help message and exit
  --analysis ANALYSIS
  --type TYPE
  --question QUESTION
  --subject-id SUBJECT_ID
```

### `yauvi analysis add`

```
usage: yauvi analysis add [-h] --analysis ANALYSIS --role ROLE --file FILE

options:
  -h, --help           show this help message and exit
  --analysis ANALYSIS
  --role ROLE
  --file FILE
```

### `yauvi analysis validate`

```
usage: yauvi analysis validate [-h] --analysis ANALYSIS

options:
  -h, --help           show this help message and exit
  --analysis ANALYSIS
```

### `yauvi analysis run`

```
usage: yauvi analysis run [-h] --analysis ANALYSIS

options:
  -h, --help           show this help message and exit
  --analysis ANALYSIS
```

### `yauvi analysis export`

```
usage: yauvi analysis export [-h] --analysis ANALYSIS --out OUT

options:
  -h, --help           show this help message and exit
  --analysis ANALYSIS
  --out OUT
```


---

## structqc

Coordinate trust: completeness, provenance class, imported validation.

```
usage: structqc [-h] {describe,validate,fetch,run} ...

positional arguments:
  {describe,validate,fetch,run}

options:
  -h, --help            show this help message and exit
```

### `structqc describe`

```
usage: structqc describe [-h]

options:
  -h, --help  show this help message and exit
```

### `structqc validate`

```
usage: structqc validate [-h] --structure STRUCTURE [--subject-id SUBJECT_ID]
                         [--reference-fasta REFERENCE_FASTA]
                         [--provenance PROVENANCE] [--pae PAE]
                         [--validation-report VALIDATION_REPORT]
                         [--require-external-validation] [--model MODEL]
                         [--chain CHAIN]

options:
  -h, --help            show this help message and exit
  --structure STRUCTURE
  --subject-id SUBJECT_ID
  --reference-fasta REFERENCE_FASTA
  --provenance PROVENANCE
  --pae PAE
  --validation-report VALIDATION_REPORT
  --require-external-validation
  --model MODEL
  --chain CHAIN
```

### `structqc fetch`

```
usage: structqc fetch [-h] --plan

options:
  -h, --help  show this help message and exit
  --plan
```

### `structqc run`

```
usage: structqc run [-h] --structure STRUCTURE [--subject-id SUBJECT_ID]
                    [--reference-fasta REFERENCE_FASTA]
                    [--provenance PROVENANCE] [--pae PAE]
                    [--validation-report VALIDATION_REPORT]
                    [--require-external-validation] [--model MODEL]
                    [--chain CHAIN] --out OUT

options:
  -h, --help            show this help message and exit
  --structure STRUCTURE
  --subject-id SUBJECT_ID
  --reference-fasta REFERENCE_FASTA
  --provenance PROVENANCE
  --pae PAE
  --validation-report VALIDATION_REPORT
  --require-external-validation
  --model MODEL
  --chain CHAIN
  --out OUT
```


---

## memorient

Membrane orientation and sidedness labelling.

```
usage: memorient [-h] {contexts,describe,orient,run,validate,fetch} ...

``memorient`` command line.

positional arguments:
  {contexts,describe,orient,run,validate,fetch}
    contexts            list membrane contexts
    describe            the module interface, or one membrane context if named
    orient              orient a structure and label residues
    run                 common-contract run with deterministic output names
    validate            check an input structure without orienting it
    fetch               what raw files this module needs, and where from

options:
  -h, --help            show this help message and exit
```

### `memorient contexts`

```
usage: memorient contexts [-h] [--json]

options:
  -h, --help  show this help message and exit
  --json      emit JSON
```

### `memorient describe`

```
usage: memorient describe [-h] [--json] [context]

positional arguments:
  context     context name (see `memorient contexts`)

options:
  -h, --help  show this help message and exit
  --json      emit JSON
```

### `memorient orient`

```
usage: memorient orient [-h] [--context CONTEXT] [--chain CHAIN]
                        [--topology-evidence TOPOLOGY_EVIDENCE]
                        [--out-json OUT_JSON] [--out-pdb OUT_PDB]
                        [--out-viz OUT_VIZ] [--out-pymol OUT_PYMOL]
                        [--out-html OUT_HTML] [--max-rows MAX_ROWS]
                        pdb

positional arguments:
  pdb                   path to a PDB or mmCIF file

options:
  -h, --help            show this help message and exit
  --context CONTEXT, -c CONTEXT
                        membrane context (default: gram_negative_om)
  --chain CHAIN         restrict to one chain id
  --topology-evidence TOPOLOGY_EVIDENCE
                        checksum-bound transmembrane-span JSON
  --out-json OUT_JSON   write full result JSON here
  --out-pdb OUT_PDB     write oriented PDB here
  --out-viz OUT_VIZ     write 3Dmol display JSON here
  --out-pymol OUT_PYMOL
                        write PyMOL .pml script here
  --out-html OUT_HTML   write self-contained 3Dmol.js HTML viewer here
  --max-rows MAX_ROWS   max residue rows to print (0 = all)
```

### `memorient run`

```
usage: memorient run [-h] --structure STRUCTURE [--context CONTEXT]
                     [--chain CHAIN] [--topology-evidence TOPOLOGY_EVIDENCE]
                     --out OUT

options:
  -h, --help            show this help message and exit
  --structure STRUCTURE
                        path to a PDB or mmCIF file
  --context CONTEXT, -c CONTEXT
  --chain CHAIN
  --topology-evidence TOPOLOGY_EVIDENCE
  --out OUT
```

### `memorient validate`

```
usage: memorient validate [-h] [--chain CHAIN] pdb

positional arguments:
  pdb            path to a PDB or mmCIF file

options:
  -h, --help     show this help message and exit
  --chain CHAIN  restrict to one chain id
```

### `memorient fetch`

```
usage: memorient fetch [-h] [--plan]

options:
  -h, --help  show this help message and exit
  --plan      print the declared sources only
```


---

## state-atlas

Conformational-state resemblance against declared references.

```
usage: state-atlas [-h] {describe,fetch,validate,run} ...

positional arguments:
  {describe,fetch,validate,run}

options:
  -h, --help            show this help message and exit
```

### `state-atlas describe`

```
usage: state-atlas describe [-h]

options:
  -h, --help  show this help message and exit
```

### `state-atlas fetch`

```
usage: state-atlas fetch [-h] --plan

options:
  -h, --help  show this help message and exit
  --plan
```

### `state-atlas validate`

```
usage: state-atlas validate [-h] --reference-set REFERENCE_SET
                            [--alignment-map ALIGNMENT_MAP]

options:
  -h, --help            show this help message and exit
  --reference-set REFERENCE_SET
  --alignment-map ALIGNMENT_MAP
```

### `state-atlas run`

```
usage: state-atlas run [-h] --manifest MANIFEST --reference-set REFERENCE_SET
                       [--alignment-map ALIGNMENT_MAP] [--structure STRUCTURE]
                       [--topology TOPOLOGY] [--trajectory TRAJECTORY]
                       [--chain CHAIN] [--selection SELECTION]
                       [--stride STRIDE] [--pbc {none,unwrap}]
                       [--cluster-cutoff-A CLUSTER_CUTOFF_A]
                       [--collective-variables COLLECTIVE_VARIABLES] --out OUT

options:
  -h, --help            show this help message and exit
  --manifest MANIFEST
  --reference-set REFERENCE_SET
  --alignment-map ALIGNMENT_MAP
  --structure STRUCTURE
  --topology TOPOLOGY
  --trajectory TRAJECTORY
  --chain CHAIN
  --selection SELECTION
  --stride STRIDE
  --pbc {none,unwrap}
  --cluster-cutoff-A CLUSTER_CUTOFF_A
  --collective-variables COLLECTIVE_VARIABLES
  --out OUT
```


---

## site-context

Functional-site roles, cofactors, ligands, and pockets.

```
usage: site-context [-h] {describe,fetch,validate,run} ...

positional arguments:
  {describe,fetch,validate,run}

options:
  -h, --help            show this help message and exit
```

### `site-context describe`

```
usage: site-context describe [-h]

options:
  -h, --help  show this help message and exit
```

### `site-context fetch`

```
usage: site-context fetch [-h] --plan

options:
  -h, --help  show this help message and exit
  --plan
```

### `site-context validate`

```
usage: site-context validate [-h] --manifest MANIFEST --structure STRUCTURE
                             --annotations ANNOTATIONS
                             [--component-map COMPONENT_MAP]
                             [--pocket-result POCKET_RESULT]

options:
  -h, --help            show this help message and exit
  --manifest MANIFEST
  --structure STRUCTURE
  --annotations ANNOTATIONS
  --component-map COMPONENT_MAP
  --pocket-result POCKET_RESULT
```

### `site-context run`

```
usage: site-context run [-h] --manifest MANIFEST --structure STRUCTURE
                        --annotations ANNOTATIONS
                        [--component-map COMPONENT_MAP]
                        [--pocket-result POCKET_RESULT] --out OUT

options:
  -h, --help            show this help message and exit
  --manifest MANIFEST
  --structure STRUCTURE
  --annotations ANNOTATIONS
  --component-map COMPONENT_MAP
  --pocket-result POCKET_RESULT
  --out OUT
```


---

## actstate

Activity-state evidence assembly.

```
usage: actstate [-h] [--version] {run,validate,describe,fetch} ...

Classify whether evidence supports a protein being in a working state.

positional arguments:
  {run,validate,describe,fetch}
    run                 assess every protein and write the results
    validate            check inputs without running
    describe            print the machine-readable IO contract
    fetch               plan the raw files this module needs

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```

### `actstate describe`

```
usage: actstate describe [-h]

options:
  -h, --help  show this help message and exit
```

### `actstate fetch`

```
usage: actstate fetch [-h] [--plan]

options:
  -h, --help  show this help message and exit
  --plan      print the declared sources only
```

### `actstate validate`

```
usage: actstate validate [-h] [--in INPUT] [--annotation ANNOTATION]
                         [--fasta FASTA] [--structures STRUCTURES]

options:
  -h, --help            show this help message and exit
  --in INPUT            input directory, or an annotation table
  --annotation ANNOTATION
                        UniProt-style annotation TSV/CSV
  --fasta FASTA         sequences, if not in the annotation table
  --structures STRUCTURES
                        directory of PDB/mmCIF files
```

### `actstate run`

```
usage: actstate run [-h] [--in INPUT] [--annotation ANNOTATION]
                    [--fasta FASTA] [--structures STRUCTURES] --out OUTPUT
                    [--chain CHAIN] [--fold-state FOLD_STATE]
                    [--expected-residues EXPECTED_RESIDUES]
                    [--reference-comparison REFERENCE_COMPARISON]
                    [--max-separation MAX_SEPARATION]

options:
  -h, --help            show this help message and exit
  --in INPUT            input directory, or an annotation table
  --annotation ANNOTATION
                        UniProt-style annotation TSV/CSV
  --fasta FASTA         sequences, if not in the annotation table
  --structures STRUCTURES
                        directory of PDB/mmCIF files
  --out OUTPUT          output directory
  --chain CHAIN         restrict geometry to one chain
  --fold-state FOLD_STATE
                        JSON of fold_state records, keyed by accession
  --expected-residues EXPECTED_RESIDUES
                        JSON of expected catalytic residues, keyed by
                        accession then position. Required to reach
                        active_site_disrupted.
  --reference-comparison REFERENCE_COMPARISON
                        JSON of reference-state comparisons, keyed by
                        accession
  --max-separation MAX_SEPARATION
                        active-site cluster bound in angstrom (default 16.0)
```


---

## assembly-context

Biological assembly, stoichiometry, contacts, and burial.

```
usage: assembly-context [-h] {describe,fetch,validate,run} ...

positional arguments:
  {describe,fetch,validate,run}

options:
  -h, --help            show this help message and exit
```

### `assembly-context describe`

```
usage: assembly-context describe [-h]

options:
  -h, --help  show this help message and exit
```

### `assembly-context fetch`

```
usage: assembly-context fetch [-h] --plan

options:
  -h, --help  show this help message and exit
  --plan
```

### `assembly-context validate`

```
usage: assembly-context validate [-h] --manifest MANIFEST --isolated ISOLATED
                                 --assembly ASSEMBLY --subject-chain
                                 SUBJECT_CHAIN --relationship
                                 {exact_protein,homolog_assembly,architecture_analogy,unresolved}
                                 [--reference-id REFERENCE_ID]
                                 [--assembly-id ASSEMBLY_ID]
                                 [--expected-chains EXPECTED_CHAINS]

options:
  -h, --help            show this help message and exit
  --manifest MANIFEST
  --isolated ISOLATED
  --assembly ASSEMBLY
  --subject-chain SUBJECT_CHAIN
  --relationship {exact_protein,homolog_assembly,architecture_analogy,unresolved}
  --reference-id REFERENCE_ID
  --assembly-id ASSEMBLY_ID
  --expected-chains EXPECTED_CHAINS
                        comma-separated chain ids
```

### `assembly-context run`

```
usage: assembly-context run [-h] --manifest MANIFEST --isolated ISOLATED
                            --assembly ASSEMBLY --subject-chain SUBJECT_CHAIN
                            --relationship
                            {exact_protein,homolog_assembly,architecture_analogy,unresolved}
                            [--reference-id REFERENCE_ID]
                            [--assembly-id ASSEMBLY_ID]
                            [--expected-chains EXPECTED_CHAINS] --out OUT

options:
  -h, --help            show this help message and exit
  --manifest MANIFEST
  --isolated ISOLATED
  --assembly ASSEMBLY
  --subject-chain SUBJECT_CHAIN
  --relationship {exact_protein,homolog_assembly,architecture_analogy,unresolved}
  --reference-id REFERENCE_ID
  --assembly-id ASSEMBLY_ID
  --expected-chains EXPECTED_CHAINS
                        comma-separated chain ids
  --out OUT
```


---

## sf-csa

Structure- and sequence-based functional comparison.

```
usage: sf-csa [-h] [--version]
              {run,verify,build-manifests,describe,validate,fetch} ...

Structure-Function Comparative Species Analysis

positional arguments:
  {run,verify,build-manifests,describe,validate,fetch}
    run                 run the comparison and write a release
    verify              audit an existing release
    build-manifests     build checksum-pinned manifests from a campaign spec
    describe            print the machine-readable IO contract
    validate            check the manifests without running
    fetch               what raw files this module needs, and where from

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```

### `sf-csa describe`

```
usage: sf-csa describe [-h]

options:
  -h, --help  show this help message and exit
```

### `sf-csa validate`

```
usage: sf-csa validate [-h] --queries QUERIES --databases DATABASES

options:
  -h, --help            show this help message and exit
  --queries QUERIES
  --databases DATABASES
```

### `sf-csa fetch`

```
usage: sf-csa fetch [-h] [--plan]

options:
  -h, --help  show this help message and exit
  --plan      print the declared sources only
```

### `sf-csa run`

```
usage: sf-csa run [-h] --queries QUERIES --databases DATABASES --output OUTPUT

options:
  -h, --help            show this help message and exit
  --queries QUERIES
  --databases DATABASES
  --output OUTPUT
```

### `sf-csa verify`

```
usage: sf-csa verify [-h] --output OUTPUT --databases DATABASES

options:
  -h, --help            show this help message and exit
  --output OUTPUT
  --databases DATABASES
```

### `sf-csa build-manifests`

```
usage: sf-csa build-manifests [-h] --spec SPEC --out OUT

options:
  -h, --help   show this help message and exit
  --spec SPEC  campaign spec JSON
  --out OUT    directory to write the manifests into
```


---

## yauvi-fetch

Registered public-source acquisition and staging.

```
usage: yauvi-fetch [-h] [--registry REGISTRY] [--cache CACHE]
                   {plan,get,stage,verify,where,sources} ...

Plan, acquire, and verify the raw input files a module declares.

positional arguments:
  {plan,get,stage,verify,where,sources}
    plan                report what a module needs and what is present
    get                 retrieve the sources policy permits
    stage               adopt a file you acquired by hand
    verify              re-hash cached files against their manifests
    where               print the cached path for a source
    sources             list the declared registry

options:
  -h, --help            show this help message and exit
  --registry REGISTRY   path to catalogs/sources.yaml
  --cache CACHE         source cache directory (default:
                        ~/.cache/yauvi/sources)
```

### `yauvi-fetch sources`

```
usage: yauvi-fetch sources [-h] [--channel CHANNEL] [--json]

options:
  -h, --help         show this help message and exit
  --channel CHANNEL  filter by channel, e.g. localization
  --json
```

### `yauvi-fetch plan`

```
usage: yauvi-fetch plan [-h] --for MODULE [--manifest MANIFEST] [--json] [-v]

options:
  -h, --help           show this help message and exit
  --for MODULE         module id, e.g. subproteo
  --manifest MANIFEST  explicit path to the module's sources.yaml
  --json               machine-readable output
  -v, --verbose
```

### `yauvi-fetch get`

```
usage: yauvi-fetch get [-h] --for MODULE [--manifest MANIFEST]
                       [--source-id SOURCE_ID] [--arg SOURCE_ID=VALUE]
                       [--run-dir RUN_DIR] [--no-probe]

options:
  -h, --help            show this help message and exit
  --for MODULE          module id, e.g. subproteo
  --manifest MANIFEST   explicit path to the module's sources.yaml
  --source-id SOURCE_ID
                        fetch only this source
  --arg SOURCE_ID=VALUE
                        identifier for a source that needs one, e.g.
                        uniprot_proteomes=UP000005640
  --run-dir RUN_DIR     append acquisitions to this platform run ledger
  --no-probe            skip the up-front reachability check and attempt every
                        retrieval
```

### `yauvi-fetch stage`

```
usage: yauvi-fetch stage [-h] [--note NOTE] [--run-dir RUN_DIR] source_id path

positional arguments:
  source_id
  path

options:
  -h, --help         show this help message and exit
  --note NOTE        how it was obtained
  --run-dir RUN_DIR
```

### `yauvi-fetch verify`

```
usage: yauvi-fetch verify [-h] [--source-id SOURCE_ID] [-v]

options:
  -h, --help            show this help message and exit
  --source-id SOURCE_ID
  -v, --verbose
```

### `yauvi-fetch where`

```
usage: yauvi-fetch where [-h] source_id

positional arguments:
  source_id

options:
  -h, --help  show this help message and exit
```

