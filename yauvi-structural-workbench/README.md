# YAUVI Structural Biology Platform — Mark 1

**YAUVI Structural Biology Platform — Mark 1** is the primary integrated,
local-first experience for the **YAUVI Structural Workbench** Python suite and
its six evidence-bounded structural-protein analyses:

1. Structure provenance and coordinate quality (`structqc`)
2. Membrane orientation (`memorient`)
3. Reference-bounded conformational resemblance (`state-atlas`)
4. Functional-site and catalytic-competence evidence (`site-context` with `actstate`)
5. Biological assembly and interface context (`assembly-context`)
6. Separate structural and sequence comparison (`sf-csa`)

Each package remains independently installable. The root distribution provides
one reviewer install without copying scientific implementations.

The platform display identity is frozen in [PLATFORM_IDENTITY.json](PLATFORM_IDENTITY.json).
Start and share guidance is in [START_HERE.md](START_HERE.md). Package names,
CLI commands, module IDs, and scientific contracts remain unchanged.

> **Scientific boundary:** the suite reports measurements, provenance,
> uncertainty, missing evidence, and explicit non-claims. It does not calculate a
> universal protein, activity, interface, or function score. Structural
> resemblance is not biochemical activity; modeled membrane orientation is not
> intact-cell exposure; docking is not part of this release.

## Install from the repository

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

The default scientific workflows remain offline. Add the `source-fetch` extra
(`".[dev,source-fetch]"`) only if you intend to use registered public-accession
retrieval; it installs an HTTP client and enables nothing by itself.

## Run an analysis

```bash
structqc run \
  --structure software/structqc/examples/model.pdb \
  --reference-fasta software/structqc/examples/reference.fasta \
  --provenance software/structqc/examples/provenance.json \
  --validation-report software/structqc/examples/validation.json \
  --out qc-demo
```

Every scientific capability is reachable this way. See
[docs/cli-reference.md](docs/cli-reference.md) for all ten commands and
[docs/quickstart.md](docs/quickstart.md) for a guided run.

Public-accession retrieval stays disabled unless a command is given
`--allow-reference-fetch`, and then only for registered providers. No arbitrary
URL, private sequence upload, or caller-supplied filesystem path is accepted.
Acquired files enter a checksum cache and must be explicitly adopted into a new
analysis revision.

### Scope of this distribution

The current local development build packages both the command-line tools and a
standalone loopback browser. Start it with `yauvi --workspace ./my-analysis workbench open`, which starts the application and opens a browser on it; `workbench serve` starts the same application without opening one.
Both interfaces use the same structural analysis store and scientific engines.
The published commit recorded in the implementation baseline predates this extraction;
publication of these changes remains pending review.

## Five-minute offline example

See [Quickstart](docs/quickstart.md). The shipped example uses synthetic
coordinates and requires no public or private data.

## Five tested human use cases

The [plain-language showcase](../evidence/showcase/five-human-use-cases/SHOWCASE.html)
contains five actual CLI executions covering coordinate trust, membrane
sidedness, conformational resemblance, functional-site mapping, and assembly
interfaces. Every benefit is paired with an explicit non-claim, and the raw
JSON/TSV evidence remains inspectable. Rebuild and verify it with:

```bash
python tools/build_five_use_case_showcase.py --replace
python tools/verify_five_use_case_showcase.py
```

## Public evidence showcase

The [public narrative microsite](../evidence/public-showcase/index.html) explains all six
workflows for non-specialists, exposes five executed synthetic analyses and a
checksum-bound SF-CSA process-boundary case, and then presents the separate
six-workflow public qualification collection. It includes accepted file types,
biological context, failed checks, raw JSON/TSV evidence, source locks, a
reviewer quickstart, and the gated JOSS publication roadmap. The page can be
opened directly from disk or served at `/public-showcase/` by the local
controller. It has no analytics, remote assets, automatic network retrieval, or
upload path.

```bash
python tools/build_five_use_case_showcase.py --replace
python tools/verify_public_showcase.py
```

## Documentation

- [Installation](docs/install.md)
- [Six workflows](docs/workflows.md)
- [Files and official sources](docs/files-and-sources.md)
- [Scientific methods and limitations](docs/methods-and-limitations.md)
- [Benchmarks and qualification](docs/benchmarks.md)
- [Reproducibility](docs/reproducibility.md)
- [Reviewer quickstart](docs/reviewer-quickstart.md)
- [Code walkthrough](docs/code-walkthrough.md)
- [Cross-platform testing and file sharing](docs/cross-platform-testing.md)
- [Pre-public audit](PREPUBLIC_AUDIT.md)
- [JOSS readiness](RELEASE_STATUS.json)
- [JOSS publication roadmap](JOSS_PUBLICATION_ROADMAP.json)
- [JOSS preparation checklist](JOSS_CHECKLIST.md)
- [Research-use and human-review handoff](docs/JOSS_RESEARCH_HANDOFF.md)
- [Recorded offline baseline](BASELINE.json)
- [Independent public qualification v1](../evidence/benchmarks/qualification-v1/README.md)
- [Mark 1 Qualification v2 panel and current gaps](../evidence/benchmarks/qualification-v2/README.md)

## Current release status

The repository is public and remains in development. Local implementation work is
tracked in [IMPLEMENTATION.md](IMPLEMENTATION.md). Collection 2.11 adopts 94 of
110 required records across six panels. Public workflow run 35295451246 recorded
passing results for all five release-blocking panels on commit
`7981148e70c5eaa6424e3608f09cd65bcfb35834`; the retained aggregate execution
summary reports 67/110 across five panels and predates the later SF-CSA result.
These records describe historical runs, not the changed local candidate, which
requires fresh review and qualification. Membrane orientation remains
non-blocking and research-only, with no accuracy claim. The public-data
structural-portfolio case study is complete; the author must review it and decide
whether it documents qualifying developer research use. The manuscript and
current candidate still require review, and the project must meet the
more-than-six-month active public-development criterion. Developer use can meet
the JOSS research-use criterion; independent adoption is not required. See the
[research-use handoff](docs/JOSS_RESEARCH_HANDOFF.md).

The pre-public audit recorded two code-interpretation decisions that had to be
resolved or explicitly narrowed before submission. Both are now narrowed:
ActState's generic catalytic-residue screen can no longer reach
`active_site_disrupted` without a position-specific expected residue
(2026-09-02), and SF-CSA computes reciprocal-best-hit evidence before structural
classification so it reaches the `probable_same_function` gate by measurement
rather than by a manifest field (2026-09-01). The audit's remaining open items
stand. The local hardening now marks declared cofactor occupancy unavailable when
exact identity and site proximity are missing; it does not treat heteroatom presence
as a supported occupancy signal.

## License and citation

Original YAUVI code is staged under Apache-2.0; public release remains blocked
until the third-party-asset audit passes. See [CITATION.cff](CITATION.cff) and
[NOTICE.md](../NOTICE.md). Public source data and optional runtimes retain their own
licenses and are never silently redistributed.
