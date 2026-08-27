# Methods and limitations

The browser is an orchestration and reporting layer. Scientific calculations
remain in the standalone CLIs, and exit codes retain their common meaning: `0`
completed, `1` scientifically incomplete, and `2` invalid input/configuration.

No workflow emits a universal score. Numeric values retain their method, units,
input checksum, parameter record, and evidence class. Presentation precision
does not alter raw JSON or TSV values.

Important boundaries:

- Unknown provenance remains unknown.
- Missing reference sequence means completeness is unevaluated.
- Predicted confidence is not experimental validation.
- Active-like/inactive-like is not biochemical activity.
- Alpha-helical membrane placement is experimental and outside the Mark 1
  qualified scope.
- Annotated catalytic residues are not observed catalysis.
- Assembly burial is not intact-cell accessibility.
- Structural similarity and sequence homology are separate evidence legs.
- Missing runtimes or references disable their named evidence leg; they never
  create a favorable value.

Readiness is scope-specific. A tool can have a conditionally qualified
beta-barrel path and a prototype alpha-helical path at the same time. Software
tests, public named-case agreement, expanded-panel qualification, runtime
availability, and independent reproduction remain separate fields.

See each run's `RUN_MANIFEST.json`, `REPORT_DATA.json`, claim ceiling, and
missing-evidence list for the exact boundary that applied.
