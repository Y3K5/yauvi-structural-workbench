# Development tester candidate — 2026-09-19

This record describes local checks, not independent tester outcomes or scientific validation.

## Intended use

- Guided product testing of six worked examples, including uncertainty, provenance and limitations.
- Full-workbench installation, local navigation and the bundled synthetic coordinate-quality example.
- A second synthetic run that deliberately omits validation evidence.
- New scientific inputs only with a reviewed workflow-specific bundle and its required runtimes.

## Verification status

- Environment: macOS on Apple silicon, Python 3.12. Application wheel installed into a separate temporary directory without downloading dependencies; existing local scientific dependencies were used. This verifies packaged code, not a clean-machine installation.
- Focused workbench, launcher and endpoint-geometry checks: 30 passed, one local-socket check initially skipped by the sandbox. That exact socket check was subsequently rerun with local binding permitted and passed (31 distinct checks passed overall).
- Installed-wheel synthetic example: creation, readiness, execution and report export passed. The deliberately missing-validation example finished as `scientifically_incomplete` with the expected nonzero run exit code, and exported its report successfully. This is expected evidence handling, not a calculation crash.
- Browser checks: refreshed showcase and installed-workbench welcome page rendered; plain-language help opened; the incomplete-evidence finding named the missing geometry validation; the workbench showcase link loaded. Showcase animation started and arrow-key navigation changed panels. No browser console errors were observed in the checked views.
- Wheel inspection: 95 members, restricted to the intended structural application and scientific namespaces. The builder screens the actual payload and every wheel member for local paths, selected private-project names and credential patterns; see CONTENT_SCREEN.json. It verifies archived file hashes and ZIP integrity.
- The tester welcome HTML was reviewed as source. Direct local-file browser inspection was blocked by the browser tool's URL policy; welcome-page rendering, feedback download, small-screen layout and offline file-opening behavior remain tester checks. No workaround was used.

## Limits

- Windows and Linux fresh installation are not tested in this session.
- Dependency download on a fresh machine is not yet demonstrated. The package does not contain dependency wheels or Python.
- Native runtimes and reference databases for all six analysis workflows are not bundled.
- Independent user testing and independent scientific reproduction have not been performed by preparing this kit.
- The conformational-state illustration is limited to the two supplied ABL endpoints. It is not molecular dynamics.
- Membrane orientation is research-only. Fold similarity does not establish equivalent function.
- Use the checksum manifest to identify this candidate; a version number alone does not identify its exact bytes.

## Suggested acceptance criteria for a pilot

A tester can find the relevant example, identify its input and evidence, explain
one limitation, and save feedback without coaching. For the executable path, a
tester can launch, run the bundled synthetic case, explain the missing-evidence
variant and reopen an exported report. Record setup failures separately from
analysis failures. These are proposed criteria, not measured outcomes.
