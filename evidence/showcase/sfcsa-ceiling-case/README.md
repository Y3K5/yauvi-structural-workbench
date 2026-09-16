# SF-CSA public showcase case

This generated case runs the canonical SF-CSA pipeline through the reviewed
offline fixture's deterministic Foldseek and DIAMOND test doubles. The stubs
compute no alignments. The case demonstrates pipeline wiring, parsing, evidence
separation, classification, fail-closed release verification, and reproducible
artifacts—not biological performance or external scientific qualification.

Rebuild and verify from the repository root:

```bash
python tools/build_sfcsa_showcase_case.py --replace
python tools/verify_sfcsa_showcase_case.py
```
