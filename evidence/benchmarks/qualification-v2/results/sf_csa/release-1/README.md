# release-1 — superseded sf-csa evidence, retained deliberately

**This release does not describe the current sf-csa panel.** It was produced on
2026-09-01 under the query selection that collection 2.7 replaced the following
day. Read it as a record of what was withdrawn, not as a result.

It is kept rather than deleted because the adoption protocol says to preserve
what you withdrew. Deleting it would remove the evidence behind a documented
selection change; leaving it unlabelled would let it be mistaken for current.
This file is the label.

## What changed under it

Collection 2.7 replaced **five of the twelve queries**:

| retired | replaced by |
|---|---|
| P00193, P23370, P45850, P00138, P00147 | P00208, P00818, O00625, Q9A980, P43934 |

The reason is recorded in `PANEL_MANIFEST.json` under `threshold_revisions`,
collection 2.7. In short: the previous selection was made against a gate that
checked whether an entry's *organism* had a reference proteome, rather than
whether the entry is *in* one. Five of the twelve were in no reference proteome
at all and a sixth pointed at a proteome containing zero occurrences of it, so
the sequence leg was absent for thirteen of sixteen records — including two
exact self-matches. The queries here are not wrong structures; they are
structures whose sequence evidence could not be computed.

`P26394` is also present, and its proteome pointer was corrected in 2.7
(UP000002695 → UP000001014) rather than the query being replaced.

Seven of the twelve queries carried forward unchanged: P00198, P32081, P0A9X9,
P00817, P26394, P37610, Q07688.

## The eight superposition.html files do not render

Each `targets/<accession>/superposition.html` embeds two structures as PDB text
and then loads a viewer:

    <script src="../../assets/vendor/3Dmol-min.js"></script>

**No release has ever contained that file.** `sf_csa.core.write_superposition_html`
hardcodes the reference and nothing writes the asset, so every one of these
pages opens blank. The embedded coordinates are intact and readable in a text
editor; only the rendering is missing.

They are not patched here. `CHECKSUMS.json` in this directory covers all 112
files including these eight, and editing frozen evidence to make it look better
is the opposite of what the checksums are for. The generator is the thing worth
fixing, and a fix there changes future releases without rewriting this one.

## What is trustworthy here

The comparison matrix, the per-target TSV evidence, the release manifest and the
checksums are all internally consistent and verify against `CHECKSUMS.json`. As
a record of *what the pipeline did on 2026-09-01 with that selection*, this
directory is sound. It is superseded as a scientific result, not corrupted as an
artefact.

## This file is not in CHECKSUMS.json, on purpose

`CHECKSUMS.json` records the 112 files the pipeline produced and publishes.
This README is documentation written afterwards, not evidence the run generated,
so adding it to the digest set would blur the line the file exists to draw. The
published directory therefore holds 113 files against 112 recorded digests, and
that is intended rather than an omission.

A local checkout will show more: `work/` holds pipeline scratch and is excluded
by `.gitignore`, so it is neither published nor digested.

## Where the current evidence lives

The sf-csa panel holds no adopted records. Its curation is in
`ADOPTION_DRAFT_SFCSA.json`, its findings in `SF_CSA_PREADOPTION_FINDINGS.md`,
and it does not execute in CI: Finding 8 records that ten of its reference
proteomes are locked against UniProt stream queries that no sha256 can lock, and
those entries were withdrawn from `SOURCE_LOCK.json` on 2026-09-04.
