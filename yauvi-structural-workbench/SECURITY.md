# Security policy

## Scope

The workbench is local-first. Its optional loopback interface binds to the
loopback address, checks Host and Origin, requires a session token for
mutations, invokes registered commands only, and exposes no filesystem browser,
arbitrary URL fetcher, or external upload path.

Security-relevant classes for this project include path traversal, command
injection, arbitrary-file access, unsafe archive expansion, source-host bypass,
checksum substitution, and **silent scientific corruption** — wrong residue
identity, mismatched checksums, permissive missing-evidence behavior, or claims
widened beyond their recorded ceiling. The last class is treated with the same
severity as a memory-safety bug, because a quietly wrong evidence record is more
damaging here than a crash.

## Reporting a vulnerability

Report privately through GitHub's private vulnerability reporting:

**https://github.com/Y3K5/yauvi-structural-workbench/security/advisories/new**

That channel is private between you and the maintainer until an advisory is
published. Please do not open a public issue for an unpatched vulnerability.

Include the affected command or module, the input that triggers it, what you
observed, and what you expected. A reproducing fixture is more useful than a
description. Do not attach unpublished sequences, private coordinates,
credentials, or unpublished research data to any report.

## Expectations

This is a pre-public research build maintained by one person. There is no
guaranteed response window and no security-support commitment for any release.
Reports are acknowledged as promptly as is practical.

## Supported versions

No version carries a security-support guarantee. The project is at
`0.1.0.dev0` and has made no release.

## Third parties

External tools (FreeSASA, Foldseek, DIAMOND, MolProbity, mkdssp) and provider
services (wwPDB, UniProt, AlphaFold DB, OPM, CATH, M-CSA) retain their own
security policies. Report issues in those components to their maintainers.
