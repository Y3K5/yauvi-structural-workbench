# SF-CSA reviewed artifact adoption

The Claude Science artifact was treated as review material, not as a source of
biological truth. Its archive checksum is recorded in
`SFCSA_ARTIFACT_ADOPTION.json`, and every accepted component was compared with
the canonical SF-CSA source before use.

## What was reusable

- The archive's SF-CSA `core.py`, `cli.py`, and `manifests.py` are byte-identical
  to the canonical workspace package. They were retained in place rather than
  copied over themselves.
- The 22-test offline fixture is now a recorded regression suite and supplies
  the deterministic process boundary used by public case `HUC-06`.
- The SF-CSA interpretation-ceiling probe is byte-identical to the canonical
  probe and remains a fail-closed scientific-boundary check.
- The public case runs the canonical CLI, parses its actual JSON/TSV release,
  verifies its checksum manifest, and keeps structural and sequence evidence in
  separate tables.

## What was deliberately not copied

The flat artifact export includes a different `core.py` and several Phase B
provenance files without a complete package-aware change set. Those fragments
were not reconstructed or installed. Generated reports and images were also not
treated as source code or scientific qualification evidence. Campaign-specific
interpretation defaults remain outside the organism-neutral JOSS workflow.

## Current scientific boundary

The fixture executables are deterministic Foldseek- and DIAMOND-shaped test
doubles. They compute no alignments. `HUC-06` therefore demonstrates command
construction, parsing, classification ceilings, evidence separation, checksums,
and release auditing—not alignment accuracy or biological function.

The remaining qualification work is explicit:

1. Run real checksum-pinned Foldseek and DIAMOND mini-databases.
2. Resolve how reciprocal-best-hit results enter structural classification.
3. Make title-trap protection consistent at both direct and audit boundaries.
4. Normalize missing-value behavior.
5. Review provenance-envelope changes only from a complete source tree.

## Reproduction

```bash
PYTHONPATH=sf-csa/src python -m pytest tools/fixtures/sfcsa -q
python tools/build_sfcsa_showcase_case.py --replace
python tools/verify_sfcsa_showcase_case.py
python tools/build_five_use_case_showcase.py --replace
python tools/verify_public_showcase.py
```

No command above publishes, deploys, uploads, or submits the workbench.
