# Bio-Orient delivery checks

Checked locally on 2026-09-30. These are software and geometry checks, not
biological qualification or an external release.

## Offline regression checks

GitHub preparation additionally passed **801 tests**, **15 skips**, **5
deselections**, and no failures/errors with dependency preflight enabled.
[GITHUB_PREPARATION_TESTS.json](GITHUB_PREPARATION_TESTS.json) records that run.
Six new controls check bundled engine ownership; installed-wheel preflight also
passed, including when source build metadata coexists with wheel metadata.
The JOSS staging, canonical wheel namespace, generated CLI reference, existing
public showcase and release digest checks passed. New CI jobs are prepared but
have not run on GitHub. These checks do not establish a clean-machine install.
The legacy showcase was rebuilt into a temporary directory and its validators
passed. Its builder projects only the six workflows with recorded qualification
narratives, so experimental additions cannot inherit those results.

The earlier development run below is retained with its original scope.

The distribution runner passed **795 tests**, with **15 skips**, **5 deselections**,
and no failures or errors. Bio-Orient contributes 42 synthetic controls. Coverage
includes reference-relative handedness, missing atoms/references, conformers,
proper transforms, exact view identities, scoped evidence, fixed-frame neighbor
occlusion, unresolved patches in a supported membrane frame, old cases and reports.
The machine-readable runner result is [VERIFICATION.json](VERIFICATION.json).

The command was `python tools/run_structural_workbench_tests.py --skip-preflight`
from the distribution root. Execution used ARM64 Python 3.11.15, NumPy 1.26.4,
SciPy 1.13.1, Biopython 1.84 and Gemmi 0.6.5. Test-only Hypothesis and pytest were
made available through temporary/source paths; canonical environments were not
modified. `--skip-preflight` was needed because standalone distribution metadata
was absent for the bundled `memorient`, `yauvi-assembly-context` and
`yauvi-structqc` names. This exception does not certify a fresh dependency install.

Skips cover eleven optional external StructPrep fixtures, one optional MemOrient
portal stylesheet, one optional StateAtlas fixture, one unavailable qualification
fixture and one socket check blocked in the test sandbox. Network/adapter-marked
checks were excluded. Each limitation remains separate from the passing checks.

## Installed packages

Both the workbench distribution and standalone Bio-Orient wheel were built from
copies of the final application sources and installed with source paths removed.
The temporary virtual environment shared the existing numerical dependencies;
it was not a clean machine or independent dependency installation.

Installed `describe`, `validate` and `run` commands worked. Workbench runs included
StructQC and the graph, tables, viewer layer and manifest in reports/evidence ZIPs.
A deliberately reflected synthetic center remained a visible review flag while
the workflow completed. Repeated outputs were byte-identical.

The real AQP1 graph also reproduced byte for byte through the installed CLI:
SHA-256 `3256211e505506b9a5f70c6bb95070c45d486fa72a95e24237a5a57c2ab11942`.
Its coordinates, patch declaration, assembly choice, CCD references and membrane
evidence were passed explicitly from the recorded run.

## Browser checks

Checked the running loopback preview, its actual assets and both viewers.
Unknown, supported, predicted and conflicting membrane records kept their states
and compartment labels consistent. Unsigned faces retained unsigned labels.
Patch, neighbor, component and atom selection, connected-view navigation, source
tracing and evidence export were exercised. The downloaded ZIP contained the graph,
both tables and the viewer layer. Browser console checks found no errors/warnings.

The large-complex check exposed two presentation faults: coordinate-only atoms
were styled as sticks without bond connections, and the legacy case-list response
could overwrite the connected graph's status. The viewer now renders visible
atom spheres plus mapped protein C-alpha traces, separates plane labels, outlines
the modeled planes, and provides reset, zoom, focus and background controls.
Readout rows present plain-language findings; exact fields remain expandable.
Coordinate and evidence identities are unchanged by these presentation choices.

## Real-complex scope

See [the AQP1 check](KNOWN_COMPLEX_CHECK.md). This is one public structure and
annotation agreement check. It does not qualify all membrane contexts, native
assembly completeness, ligand stereochemical methods, transport activity, probe
reachability, immunological scoring or dynamic probabilities. Private cases were
not added to the software distribution. No publication, deployment or push was
performed.
