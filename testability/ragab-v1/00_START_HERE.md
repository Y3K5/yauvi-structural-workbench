# RagAB testability: start here

[README.md](README.md) explains inputs, setup, checks, results and interpretation.
Public source files go in; fixed structural calculations and explicit failure
checks run; new private execution records and numerical comparisons come out.

Use `replay.py verify-kit` for shipped identities and `test_replay.py` for
synthetic failure controls. Scientific replay also requires checksummed source
files, the declared dependencies and a working native FreeSASA 2.1.3.

Do not alter frozen adapters, expected results or old study records to make a
run pass. Record a changed runtime or protocol as a new version.
