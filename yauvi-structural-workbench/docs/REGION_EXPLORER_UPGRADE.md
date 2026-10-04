# Region explorer upgrade candidate

This upgrade adopts the CCD quotation correction, the normal Region explorer
workflow and reusable separate review/region-replay tools. It retains exact
coordinate and reference identities, missing coverage, unknown states and
original engine findings. See [the usage guide](REGION_EXPLORER.md).

The isolated candidate is the existing Git checkout plus the explicitly selected
upgrade files. Its reviewer suite passed **872 tests**, with **15 skipped** and
**5 deselected**, recorded in [the candidate receipt](REGION_EXPLORER_CANDIDATE_TESTS.json).
The earlier [local-tree receipt](REGION_EXPLORER_TESTS.json) records 902 passes.
The 30-test difference is an older AssemblyContext test expansion in the local
tree that is outside this upgrade. These are different checked file sets; the
counts should not be combined or treated as biological validation.

The candidate wheel was built and installed in a separate temporary directory.
Its bundled synthetic region example produced the expected scientifically
incomplete QC result while generating region evidence and a report. A relocated
region replay passed. QC records were checksum-verified; the region recipe did
not recalculate QC or all other scientific engines. The existing 3Dmol runtime
and license remained byte-identical. Runtime package contents were screened.

Private study coordinates, sequences, reports, replay payloads, research Findings
and local pages are excluded. No scientific qualification is added. No public
release, tag, commit or push is established by this document.
