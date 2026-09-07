# Start YAUVI Structural Biology Platform — Mark 1

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

That single command installs the whole distribution and puts nine console
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
standalone loopback browser. Start it with `yauvi --workspace ./my-analysis workbench serve`.
Both interfaces use the same structural analysis store and scientific engines.
The published commit recorded in the implementation baseline predates this extraction;
publication of these changes remains pending review.

## What to share

**Name:** YAUVI Structural Biology Platform — Mark 1

**One-sentence description:** A local-first structural bioinformatics platform
that turns protein coordinate files into inspectable, checksum-bound evidence
across six analysis workflows.

**Status statement:** This is a public-development project with local, unpublished
hardening changes. The historical collection 2.9 execution includes five executed
panels: four required panels passed 62 cases and four controls across six runners.
Membrane orientation passed 5/16 in five runners and 4/16 in one. SF-CSA was not
executed. These are baseline measurements, not qualification of the changed code.
See [implementation progress](IMPLEMENTATION.md) and the retained release evidence.
No independent scientific approval or submission eligibility is claimed.

**Non-claim:** It is not a clinical tool, a biochemical activity assay, or a
universal protein-scoring system. Passing local tests is not JOSS acceptance and
not external scientific validation.

Use `PLATFORM_IDENTITY.json` as the source of truth for display and share labels.
Use `CITATION.cff` for software citation metadata.
