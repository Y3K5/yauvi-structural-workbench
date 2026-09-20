# Structural portfolio examples

The supported Mark 1 example in this directory is a reproducible public-data
case study. It asks how two deposited human carbonic anhydrase II structures
map to the UniProt P00918 reference sequence, using the shipped StructQC
module. It locks source accession, chain, model, downloaded bytes, and analysis
scope, and writes all downloaded data and run outputs outside this source
directory.

Start with [`PUBLIC_RESEARCH_GUIDE.md`](PUBLIC_RESEARCH_GUIDE.md). The protocol
was locked before the primary comparison; the results are descriptive
coordinate evidence, not a ligand, activity, affinity, or biological
validation claim. The earlier 1CA2 smoke check is disclosed and excluded.

From this directory, run the public case with a fresh output path:

```sh
./run_example.sh /private/tmp/yauvi-ca2-public-run-001
```

This downloads locked public inputs to a sibling temporary directory and
writes the StructQC evidence and comparison report under the requested output
directory. Network access is needed for acquisition; once the three hash-locked
files are present, the analysis is offline. The script requires StructQC to be
installed, or its source package to be on `PYTHONPATH`.

## Legacy synthetic fixture

The other files in this directory (`query.pdb`, `assembly.pdb`, the annotation,
reference, variant, and configuration fixtures) contain invented data. They
remain as historical synthetic material only. Their former ten-module runner
targeted packages that are not part of the current Mark 1 distribution; those
files do not demonstrate biological function, and they are not the supported
Mark 1 example. The `run_example.sh` entry point now runs the public CA II case.
No legacy fixture data was deleted or used in the public case study.
