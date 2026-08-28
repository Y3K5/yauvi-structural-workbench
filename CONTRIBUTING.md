# Contributing

Contributions must preserve evidence and claim boundaries.

- Missing evidence must never become favorable evidence.
- New input formats need a validator, format guide, fixture, and failure tests.
- New public retrieval needs an allowlisted provider, bounded transport,
  identifier validation, license record, checksum provenance, and a zero-network
  default test.
- Scientific methods need invariance tests and an external benchmark record;
  unit tests do not qualify a method.
- External binaries and datasets are adapters or checksum manifests, not
  silently bundled assets.
- Third-party benchmark data is never committed. Record it in the relevant
  `SOURCE_LOCK.json` with provider, URL, and SHA-256 so others acquire the
  identical bytes themselves.
- Bug fixes require regression tests and a changelog entry.

## Setup and tests

```bash
python -m pip install -e ".[dev]"
python tools/run_structural_workbench_tests.py
```

Use that runner, not bare `pytest`. Suites are invoked in separate processes
because several contain duplicate test-module and `conftest` basenames, and a
plain `pytest` invocation fails on the collision.

Expect 495 passed, 6 network/adapter deselected, 1 skipped.

If you change a CLI's arguments, regenerate the command reference and commit it:

```bash
python tools/build_cli_reference.py
```

CI runs `python tools/build_cli_reference.py --check` and fails if it is stale.

## Opening work

Open issues and pull requests at
https://github.com/Y3K5/yauvi-structural-workbench/issues. Discussing a design
change in an issue before implementing it is usually faster than a large
unsolicited pull request, particularly for anything touching an evidence
contract or a claim ceiling.

Conduct expectations are in [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md);
vulnerabilities go through [`SECURITY.md`](SECURITY.md), not public issues.
