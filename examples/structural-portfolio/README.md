# Synthetic structural portfolio example

This fixture contains invented coordinates, annotations, references, variants, and alignments. It supports software reproduction only and makes no biological claim.

Run from anywhere:

```sh
./examples/structural-portfolio/run_example.sh /tmp/yauvi-structural-example
```

The script accepts scientific exit code `1` as a completed analysis with unresolved evidence. Exit code `2` still stops the workflow. It exercises all ten standalone packages: the five structural core modules plus structural relationships, docking preparation, pose evidence, constrained design, and oral context. The final `STRUCTURAL_REPORT.html` is self-contained and uses the canonical `yauvi-fold/4` bridge.

The docking fixture prepares and analyzes local synthetic records only; it does not claim that AutoDock Vina or HADDOCK3 ran. Public scientific benchmark execution remains a separate checksum-pinned gate.
