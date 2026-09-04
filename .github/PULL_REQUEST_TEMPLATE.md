## What this changes

<!-- One or two sentences. Link the issue if there is one. -->

## Evidence and claim boundaries

Contributions must preserve these. Confirm each that applies:

- [ ] Missing evidence still cannot become favorable evidence.
- [ ] No claim is widened beyond its recorded ceiling.
- [ ] New input formats ship a validator, format guide, fixture, and failure tests.
- [ ] New public retrieval is allowlisted, identifier-validated, checksum-recorded,
      license-recorded, and has a zero-network default test.
- [ ] No third-party benchmark data is committed; sources are recorded in a
      `SOURCE_LOCK.json` with provider, URL, and SHA-256.
- [ ] Scientific method changes carry invariance tests and an external benchmark
      record. Unit tests alone do not qualify a method.
- [ ] Bug fixes carry a regression test and a changelog entry.

## Tests

```
python tools/run_structural_workbench_tests.py
```

Result: <!-- e.g. 526 passed, 6 deselected, 1 skipped -->

- [ ] If a CLI's arguments changed, I regenerated `docs/cli-reference.md`
      with `python tools/build_cli_reference.py`.
