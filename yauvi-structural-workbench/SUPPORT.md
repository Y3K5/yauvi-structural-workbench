# Support

## Start with the documentation

- [`START_HERE.md`](START_HERE.md) — install and first run
- [`docs/quickstart.md`](docs/quickstart.md) — a complete offline analysis you can run
- [`docs/cli-reference.md`](docs/cli-reference.md) — every command in all nine CLIs
- [`docs/methods-and-limitations.md`](docs/methods-and-limitations.md) — what results may and may not be taken to mean
- [`docs/reproducibility.md`](docs/reproducibility.md) — determinism and provenance contract

## Asking for help

Open an issue: https://github.com/Y3K5/yauvi-structural-workbench/issues

Useful reports include the command you ran, the exit code, the relevant
`RUN_MANIFEST.json`, and your platform and Python version. Exit code `1` means
*scientifically incomplete*, not failure — the run manifest names what evidence
was missing, and that list is usually the answer.

**Treat silent scientific corruption as high severity** — wrong residue
identity, a mismatched checksum, permissive missing-evidence behavior, or a
claim widened beyond its recorded ceiling. If the issue is a vulnerability
rather than a bug, use the private route in [`SECURITY.md`](SECURITY.md)
instead.

Do not attach unpublished sequences, private coordinates, credentials, or local
omics data to a public report. Reduce the problem to a synthetic or public
fixture first.

## What this project does not offer

This is a pre-public research build maintained by one person. There is no
service-level commitment, no guaranteed response time, and no scientific
consulting. Questions about interpreting a result for your own research are
welcome as issues but may go unanswered.
