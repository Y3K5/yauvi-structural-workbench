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
  --structure structqc/examples/model.pdb \
  --reference-fasta structqc/examples/reference.fasta \
  --provenance structqc/examples/provenance.json \
  --validation-report structqc/examples/validation.json \
  --out qc-demo
```

Every scientific capability is reachable this way. See
[docs/cli-reference.md](docs/cli-reference.md) for all nine commands and
[docs/quickstart.md](docs/quickstart.md) for a guided run.

Public-accession retrieval stays disabled unless a command is given
`--allow-reference-fetch`, and then only for registered providers. No arbitrary
URL, private sequence upload, or caller-supplied filesystem path is accepted.
Acquired files enter a checksum cache and must be explicitly adopted into a new
analysis revision.

### Scope of this distribution

This is the command-line distribution. The loopback browser workbench is not
included: its controller imports private control-plane modules that sit outside
the published boundary. The browser layer only ever orchestrated and reported —
the calculations live in the CLIs.

## Five-minute offline example

See [Quickstart](docs/quickstart.md). The shipped example uses synthetic
coordinates and requires no public or private data.

## Five tested human use cases

The [plain-language showcase](showcase/five-human-use-cases/SHOWCASE.html)
contains five actual CLI executions covering coordinate trust, membrane
sidedness, conformational resemblance, functional-site mapping, and assembly
interfaces. Every benefit is paired with an explicit non-claim, and the raw
JSON/TSV evidence remains inspectable. Rebuild and verify it with:

```bash
python tools/build_five_use_case_showcase.py --replace
python tools/verify_five_use_case_showcase.py
```

## Public evidence showcase

The [public narrative microsite](public-showcase/index.html) explains all six
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
- [Recorded offline baseline](BASELINE.json)
- [Independent public qualification v1](benchmarks/qualification-v1/README.md)
- [Mark 1 Qualification v2 panel and current gaps](benchmarks/qualification-v2/README.md)

## Current release status

This workspace is in **pre-public preparation**. It is not a JOSS release
candidate and not submission-eligible: version-control history, independent
installation evidence, public research use, and the expanded Qualification v2
panels remain incomplete. Historical v1 evidence contains four passed named
public cases and two partial cases. The current v2 audit is
`blocked_panel_incomplete`; it has not executed the scientific panels and is not
a claim of workflow-general accuracy. Publication, repository creation, and
JOSS submission require separate approval.

The pre-public audit recorded two code-interpretation decisions that had to be
resolved or explicitly narrowed before submission. Both are now narrowed:
ActState's generic catalytic-residue screen can no longer reach
`active_site_disrupted` without a position-specific expected residue
(2026-09-02), and SF-CSA computes reciprocal-best-hit evidence before structural
classification so it reaches the `probable_same_function` gate by measurement
rather than by a manifest field (2026-09-01). The audit's remaining open items
stand, including ActState's occupancy caveat: a non-solvent heteroatom is
detected, but its identity is not proven against the declared cofactor.

## License and citation

Original YAUVI code is staged under Apache-2.0; public release remains blocked
until the third-party-asset audit passes. See [CITATION.cff](CITATION.cff) and
[NOTICE.md](NOTICE.md). Public source data and optional runtimes retain their own
licenses and are never silently redistributed.
