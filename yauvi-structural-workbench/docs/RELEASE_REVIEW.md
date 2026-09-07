# Release and publication review

Current state: local development preview; publication authorization false.

Before release, inspect the actual wheel/source archive contents and their hashes,
not only the edited source. Review bundled viewer notices and all dependency,
reference-data and executable licenses. The selected structural registry is a
public default; it does not grant rights to redistribute third-party databases.

The explicit execution-evidence projection selects status/summary records and
accounts for replaced output-directory fields. The prepared ledger standardizes
82 fields, including 51 disclosed home-directory paths from three original
execution-status files. Original provenance and historical failures are preserved.
Review the exact projection, its omissions and its destination before publishing.
Existing public history is not rewritten by this work.

Human decisions still required: exact outgoing artifacts and destination;
authorship/contributions; review of assisted output. Record tools,
versions, dates, assistance and human review without inferring missing versions.

## PUB-04 decision record — 7 September 2026

**Funding and conflicts — decided.** Yuvraj Patel declares no competing interests
and no funding of any kind. Written into `paper/paper.md`. This was the whole of
the outstanding decision; the gate is closed.

**Authorship — deliberately deferred, not unresolved.** The author list stays as
`CITATION.cff` and `paper/paper.md` record it (Yuvraj Patel, ORCID
0009-0002-2276-7336) until the USE-03 pilot participants and REP-03 independent
reproductions exist, at which point credit is revisited. JOSS expects authorship
to track substantial contribution, so testers are expected to land in
Acknowledgements rather than the author list unless their involvement goes
further. Recorded here so a reviewer sees a decision with a trigger, not a gap.

**AI model/version records — recovered from records.** Recovered by reading local
session logs, not reconstructed from memory, and no version was inferred where a
record was absent. OpenAI `gpt-5.6-sol`, `gpt-6-astra`, `gpt-reserve`,
`gpt-5.6-luna`, `gpt-5.4-mini` across nineteen sessions referencing this
repository, 25 August to 7 September 2026. Anthropic `claude-opus-5`,
`claude-sonnet-5`, `claude-opus-4-8`, 25 July to 7 September 2026. The Anthropic
range is scoped to the containing workspace rather than this repository alone and
therefore overstates what touched this work; that limit is stated in the paper
rather than smoothed over. Counts are invocation volume, not contribution.

**Licensing and third-party audit — evidenced against built artifacts.** Audited
the wheel and sdist actually produced, not the edited source, per the requirement
above. `tools/verify_structural_workbench_wheel.py` passes: 87 files, canonical
structural namespaces only, no archived or out-of-scope code. Ten canonical
top-level namespaces and none of the forbidden private ones. Three license or
notice files in each artifact. Neither artifact contains a home-directory path.

    wheel  de0331953694664f9a201da64e6a5a66a5eeedd0eefe4c38af06264531307a21
    sdist  181abe25762a244489ae5883a837cc091688978e71c41d6bbc3b185021a34590

Runtime dependency licenses, read from installed distribution metadata: gemmi
MPL-2.0; numpy and scipy BSD; PyYAML MIT; requests Apache-2.0; biopython under
the Biopython License Agreement. All permissive and compatible with distributing
this work under Apache-2.0. FreeSASA, Foldseek, DIAMOND and MDAnalysis are
invoked, never redistributed; the only `foldseek` and `diamond` files in the
repository are this project's own Python stubs under `tools/fixtures/sfcsa/`,
which are test doubles and not the tools. Third-party reference data stays out
of the repository by `.gitignore`, with `SOURCE_LOCK.json` shipping in its place.

**Reference-data provider terms — reviewed 7 September 2026.** The provider list
was taken from the `SOURCE_LOCK.json` files, which record what is actually
acquired, rather than from the prose in `NOTICE.md`. Six providers appear there:
wwPDB/RCSB PDB (133 artifacts), UniProt (50), OPM (24), AlphaFold DB (23),
M-CSA (16) and CATH (2). SIFTS, ChEBI and the PDB CCD do not appear in a source
lock but are referenced in code, so they were reviewed too. Each was checked
against its own published terms; the results are tabulated in `NOTICE.md`.

Every one is public-domain or attribution-only: CC0 1.0 for the wwPDB archive
and the CCD, CC BY 4.0 for UniProt, AlphaFold DB, M-CSA, CATH, SIFTS and ChEBI.
None restricts redistribution and none is copyleft, so none conflicts with
distributing this work under Apache-2.0. The point is largely moot for the
repository itself, which redistributes none of it.

Two findings came out of the review rather than confirming the existing record:

1. **CATH was being used and was not declared.** Two CATH artifacts sit in the
   qualification source lock while `NOTICE.md` did not name CATH at all. Added.
   CATH's CC BY 4.0 requires named attribution to its authors, which is now
   recorded.
2. **OPM does not state a licence on its own site.** Secondary sources report
   CC BY 3.0. Recorded as unconfirmed rather than asserted: the attribution
   obligation is treated as binding, the version as unverified. Confirming it
   with the University of Michigan maintainers is the one open item, and it
   affects wording, not permission — CC BY 3.0 and 4.0 both permit this use.

The CC BY licences oblige *this project* to attribute, not only its users. Where
benchmark collections and showcases present derived data, the provider is named
in the presenting record.

The draft paper describes the local changed build, which must be independently
installed/tested before that browser claim appears in a released manuscript.
The existing paper.pdf is retained as historical; paper-working-draft.pdf is the
new reviewed rendering. Tie the final manuscript to a specific released version
and its qualifying evidence, not to an unpinned working directory.

More than six months of active public development and real research use are
required under the checked [JOSS submission requirements](https://joss.readthedocs.io/en/latest/submitting.html).
Given the reviewed start of 27 August 2026, the earliest calendar checkpoint is
after 27 February 2027, conditional on actual history and all other gates.
The [AI policy](https://joss.readthedocs.io/en/latest/submitting.html#ai-usage-policy)
requires disclosure; Yuvraj handles editor/reviewer conversations, with only the
policy's translation exception. This is not a promised submission or acceptance.

Prepare final release notes, citation metadata, DOI archive materials, figures,
checksums and the visually reviewed paper only after their respective gates pass.
Yuvraj performs every push. This document and successful tests grant no outgoing
authorization.
