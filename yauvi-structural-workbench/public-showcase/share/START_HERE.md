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

This is the command-line distribution. The loopback browser workbench is **not**
included: its controller imports private control-plane modules that are outside
the published boundary. Every scientific capability is reachable from the CLIs,
which is where the calculations live in any case — the browser layer only ever
orchestrated and reported.

## What to share

**Name:** YAUVI Structural Biology Platform — Mark 1

**One-sentence description:** A local-first structural bioinformatics platform
that turns protein coordinate files into inspectable, checksum-bound evidence
across six analysis workflows.

**Status statement:** Mark 1 is a pre-public scientific build. In the historical
v1 public qualification collection, four cases passed and two remain partial.
The current release gate is Qualification v2, whose audit state is
`blocked_panel_incomplete`: **four of six panels are adopted and executed, two
are not, and no scope has reproduced on an independent second machine.** 53 of
110 required cases have been executed and passed. No scope is qualified for
release.

**Non-claim:** It is not a clinical tool, a biochemical activity assay, or a
universal protein-scoring system. Passing local tests is not JOSS acceptance and
not external scientific validation.

Use `PLATFORM_IDENTITY.json` as the source of truth for display and share labels.
Use `CITATION.cff` for software citation metadata.
