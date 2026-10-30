# Record listing and local QC integrity

The native Workbench lists saved analysis cases; Biological Orient lists its
separate local case records. Empty folders do not constitute saved cases. Existing
public readers already ignore folders without their committed manifest, so an
unwritten folder cannot prevent other saved cases from being listed.

The Biological Orient list keeps a malformed committed manifest visible as an
invalid case. Direct reads and existing run/input membership checks retain their
errors. This change does not introduce automatic repair, migration or deletion.

## What goes in and comes out

Invented saved cases, empty creation folders and synthetic QC coordinates go in →
existing native listing, execution, revision and membership checks → visible saved
cases, explicit invalid biological-case status and source-bound QC reports. A
parameter edit must not silently give an older run the new input identity, and an
input attached to one case must not be exposed through an unrelated case.

QC source files must have resident local bytes. The packaged example reference
FASTA and StructQC typing marker were already present in this distribution and
served as verified local donors for a separate installation repair. Their contents
are unchanged here. The wheel guard independently requires the bundled Workbench
QC reference and StructQC marker, supplementing package-data declarations.

## Public and private boundaries

The public package does not ship the broader Personal Research UI or StudyStore.
Its readers use their own existing analysis and biological-case contracts. The
private StudyStore repair is not copied as an unused module or represented as a
new public Research feature. No private study, input, report, custody manifest,
source placeholder, screenshot or recovery path is included.

## Repeat the invented checks

Run `python -m pytest -p no:cacheprovider -q software/platform/tests/test_record_qc_regressions.py`
in the declared project environment. The public reviewer workflow also runs
`software/platform/tests/test_structural_workbench.py` and
`software/platform/tests/test_biological_case.py`; run them when checking the
underlying contracts. `tools/test_workbench_wheel.py` checks omission of each
required asset is rejected. Tests use temporary invented libraries and existing
local engines; no download, external account or real research library is required.

Keep old scientific records and unavailable originals intact when repairing local
source availability. Prefer a verified existing local copy; a checksum or passing
software test does not establish biological correctness or scientific qualification.
This update changes no release/qualification gate and performs no data acquisition.
