# yauvi-sources

The executable half of the evidence-source registry.

`catalogs/sources.yaml` declares every external database, reference panel, and
predictor the platform can compare a protein against — what each one is, how it
is reached, what its licence permits, and what it cannot tell you. It was written
to be read. This package makes it act.

    yauvi-fetch plan   --for subproteo     what this module needs, and what is here
    yauvi-fetch get    --for subproteo     retrieve what the licence permits
    yauvi-fetch stage  deg <path>          adopt a file you obtained by hand
    yauvi-fetch verify                     re-hash the cache against its manifests
    yauvi-fetch where  uniprot_proteomes   print a cached path, for use in --in
    yauvi-fetch sources                    list the registry

## What it will not do

The governing rule is the one `shared/runtime-registry.yaml` sets for runtimes
and this layer inherits for data: **fail closed**. A source that cannot be
obtained is reported as not obtained. It is never approximated, never
substituted, and never quietly skipped.

Five classes, decided from each entry's declared `status` and `access`:

| class | what happens |
|---|---|
| `open_fetchable` | downloaded and hashed |
| `license_gated` | **never downloaded** — instructions printed, exit non-zero |
| `table_only` | you run the tool; we validate the shape of your export |
| `runtime` | an executable, not a file — resolved by runtime preflight |
| `internal` | computed in-process; nothing to acquire |

`status` is consulted before `access`, because they answer different questions.
IEDB is the case that makes this matter: it is reachable over an API, so
transport alone would mark it downloadable — but its status is `table_only` and
the registry states plainly that no code path fetches it. Honouring transport
over status would start automating a source the platform deliberately does not.

An access mode with no policy entry **raises**. It does not default to
downloadable.

## The cache

    <cache>/<source_id>/<sha256>/<filename>
    <cache>/<source_id>/MANIFEST.json

Content addressing is what makes a run reproducible. Two files differing by one
byte occupy different directories, so a silently-updated remote release can never
overwrite the copy an earlier run was built on. Both stay retrievable.

Default location `~/.cache/yauvi/sources`; override with `YAUVI_SOURCE_CACHE` or
`--cache`.

`MANIFEST.json` records, per acquisition: digest, size, origin URL or staged
path, retrieval timestamp, and the upstream release string when the endpoint
reports one. That last field closes a gap the registry admits about itself —
of UniProt it says *"the release is whatever UniProt serves that day."* UniProt
does report its release in a response header, so the fetcher captures it.

Pass `--run-dir` to also append each acquisition to a platform run ledger
(`AppendOnlyRunStore`). That is a soft dependency: this package installs and
works without `yauvi-platform`, and says so rather than failing if asked to
record into a ledger it cannot reach.

## Declaring what a module needs

Each module ships a `sources.yaml` next to its code:

```yaml
schema_version: "1.0"
module_id: subproteo
requires:
  - source_id: uniprot_proteomes
    role: "target proteome, and every panel proteome the config names"
    required: true
  - source_id: deg
    role: "stage 2 essentiality reference; licence-gated, so staged by hand"
    required: false
```

Keeping the requirement list in the module rather than in the registry is what
lets a module be planned on its own. `yauvi-fetch plan --for memorient` works
from an install with no checkout of the workspace present — only the registry
itself is needed, via `--registry` or `YAUVI_SOURCES_REGISTRY`.

`required: false` means the pipeline runs without it and records the affected
channel as unevaluated — never as a pass.

## Exit codes

These commands are meant to be run from scripts.

| code | meaning |
|---|---|
| 0 | satisfied |
| 1 | a required source is absent, or verification failed |
| 2 | usage or configuration error |

## Offline behaviour

`get` probes for a route before attempting anything, with a bound enforced from
outside the resolver call. Name resolution is not covered by a socket timeout —
`getaddrinfo` blocks — so without this an offline machine stalls on every source
in turn instead of saying immediately that it is offline. Use `--no-probe` to
attempt retrieval regardless.

`plan` never touches the network at all, and works without the `fetch` extra.

## Install

    pip install yauvi-sources            # plan, stage, verify
    pip install 'yauvi-sources[fetch]'   # adds network retrieval

## Tests

    pytest                     # offline
    pytest -m network          # adds the live-endpoint checks
