# Contributing

Contributions must preserve evidence and claim boundaries.

- Missing evidence must never become favorable evidence.
- New input formats need a validator, format guide, fixture, and failure tests.
- New public retrieval needs an allowlisted provider, bounded transport,
  identifier validation, license record, checksum provenance, and zero-network
  default test.
- Scientific methods need invariance tests and an external benchmark record;
  unit tests do not qualify a method.
- External binaries and datasets are adapters or checksum manifests, not
  silently bundled assets.
- Bug fixes require regression tests and a changelog entry.

Set up with `python -m pip install -e ".[dev]"` and run `pytest -m "not network"`.
Open design issues publicly after repository publication; publication itself is
not authorized by this staging tree.
