# Start YAUVI Structural Workbench

Mark 1 is the primary integrated experience for the six structural-analysis
workflows in the YAUVI Structural Workbench. The platform name does not rename
the standalone scientific packages or their command-line interfaces.

## Fastest local start

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

That single command installs the whole distribution and puts ten console
scripts on `PATH`. Confirm the install:

```bash
structqc describe
```

Then follow [`docs/quickstart.md`](docs/quickstart.md) for a complete offline
StructQC analysis you can run and inspect, and
[`docs/cli-reference.md`](docs/cli-reference.md) for every command.

Reference acquisition is disabled unless a command is given
`--allow-reference-fetch`. That flag enables only registered public-accession
providers; it does not authorize arbitrary URLs, uploads, or publication.

## Scope of this distribution

The current local development build packages both the command-line tools and a
standalone loopback browser. Start it with `yauvi --workspace ./my-analysis workbench open`, which starts the application and opens a browser on it; `workbench serve` starts the same application without opening one.
Both interfaces use the same structural analysis store and scientific engines.
The published commit recorded in the implementation baseline predates this extraction;
publication of these changes remains pending review.

## What to share

**Name:** YAUVI Structural Workbench

**One-sentence description:** A local-first structural bioinformatics platform
that turns protein coordinate files into inspectable, checksum-bound evidence
across six analysis workflows.

**Status statement:** This public-development project has a changed local
candidate. Qualification v2 collection 2.11 adopts 94 of 110 required records.
Public workflow run 35295451246 recorded passing outcomes for all five
release-blocking panels on commit
`7981148e70c5eaa6424e3608f09cd65bcfb35834`; those results apply to that public
commit, not this changed candidate. The older retained aggregate execution
summary reports 67/110 across five panels and predates the later SF-CSA result.
Membrane orientation remains a non-blocking, research-only scope. See the
[qualification evidence](../evidence/benchmarks/qualification-v2/README.md),
[JOSS checklist](JOSS_CHECKLIST.md), and [research-use handoff](docs/JOSS_RESEARCH_HANDOFF.md).
The manuscript and current candidate still require author review; no submission
eligibility or biological validation is claimed.

**Non-claim:** It is not a clinical tool, a biochemical activity assay, or a
universal protein-scoring system. Passing local tests is not JOSS acceptance and
not external scientific validation.

Use `PLATFORM_IDENTITY.json` as the source of truth for display and share labels.
Use `CITATION.cff` for software citation metadata.
