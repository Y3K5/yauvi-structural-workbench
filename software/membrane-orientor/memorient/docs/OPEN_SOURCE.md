# Organism Orient and MemOrient in Structural Workbench

MemOrient calculates membrane geometry, rigid receptor placement and reviewed
local molecular preparation. Organism Orient is its explorer. Workbench retains
the input files, revisions, queued jobs, results and exports. MARSE remains a
separate ecological engine; its context import is read-only.

## What goes in and comes out

- A local single-model PDB and declared plane, sphere, cylinder or ellipsoid
  produce a rigidly placed receptor and a reopenable scene. Units are Angstroms;
  placement angles are radians in files and degrees in the browser. Original
  chains and all ATOM/HETATM records are retained; assembly identity is not inferred.
- A reviewed CHARMM36-family planar profile, complete standalone TOP and its
  checksum-matched GRO produce water, neutralizing ions and a restrained minimum
  only after numerical acceptance. A placed PDB alone cannot supply chemistry.
- A portable result ZIP imports an existing result without repeating computation.
  Every artifact is bound to a relative path and SHA-256. Failed/held/interrupted
  states remain distinct; loading a result never starts a job.

Open **Organism Orient · membranes** from the Workbench sidebar, or visit
`/organism-orient/` on the same local Workbench server. The invented helix and all
four geometric shapes work without a network connection or molecular runtime.
Save a scene to reopen its original PDB, placement settings, pinned placements
and camera. Protein banks have a separate save/import action.

Use **Save current placement to Workbench** for a durable case/job. The recorded
run list reopens earlier results and exports a portable ZIP. Use the import input
to record such a ZIP. Choose **Reviewed molecular preparation** only after
reviewing the profile and matching molecular files. CPU threads default to four
and are bounded to 1–16. Progress and cancellation use Workbench's existing queue.
Restarted jobs become interrupted and require deliberate resubmission.

Minimized molecular coordinates have locked shape controls. The display shows
heavy solute; water and added monovalent ions remain in the full GRO/TOP. Generic
displays preserve saved Cartesian coordinates, so periodic molecules can cross
the box boundary. Display generation is bounded to 100,000 heavy solute atoms;
larger accepted systems retain their full downloadable artifacts with a display hold.

## Optional local runtime and source files

The core installation's dependencies are unchanged. Install the standalone
package with its existing `compute` extra to use numeric geometry. No dependency,
structure or force field is installed/downloaded by launching the explorer.

Preparation uses externally supplied GROMACS and a Python environment that can
import MemOrient, NumPy and SciPy. Configure `YAUVI_MEMORIENT_PYTHON` and
`YAUVI_MEMORIENT_GMX` on the Workbench server before starting it. These are local
server configuration, never executable fields in a browser profile. The optional
Python 3.12 conversion requirements are recorded in `requirements-md-prep.txt`;
ParmEd is only needed for separate topology conversion. Do not downgrade the main
application environment to accommodate conversion tools.

The reviewed profile must name organism/strain/layer, sources, assumptions,
force field, protein molecule types, hydration policy and checksummed inputs.
Admission is `parameterized_inputs_reviewed`, `stop_after` is `minimization`, and
`equilibration_authorized` is false. The current protocol supports three-site
CHARMM TIP3P, SOD/CLA, orthorhombic planar systems, explicit special 1–4 terms,
reviewed contact cleanup and at most thirty final ten-step minimization blocks.
Acceptance requires Fmax below 1000 kJ/mol/nm, finite energy, no recorded
constraint/fatal messages, neutral charge, matching counts/order, water geometry
and periodic intermolecular heavy-atom contact checks. These are numerical gates.

Catalogue metadata includes E. coli K-12, B. subtilis and yeast source references.
Each coordinate entry has its own availability state. Explicitly acquired,
checksum-matched view PDBs may be placed in a user-owned cache named
`<view_pdb_sha256>.pdb`; use `memorient viewer --cache CACHE_DIRECTORY` for the
standalone explorer, or set `MEMORIENT_CACHE` before starting Workbench. A cache mismatch holds that entry. You can always load your
own PDB, which remains a user assertion rather than an authenticated catalogue source.
The development checkout also retains its original local catalogue separately.

For an existing accepted legacy system, `memorient result-bundle --run RUN_DIRECTORY
--display DISPLAY_DIRECTORY --out NEW_DIRECTORY` verifies and copies its original
coordinates, topology, scene, atom map, log and receipts. It creates a local
`MEMBRANE_RESULT.zip`. Originals are unchanged. The archive contains biological
data and may contain machine paths in original receipts/logs: it is not a public
release asset. `memorient viewer --manifest MEMBRANE_RUN.json` inspects a portable
result independently of the source checkout.

## Public distribution and licensing

Project code retains Apache-2.0. The vendored 3Dmol runtime retains its BSD and
incorporated-component notices. Exact bytes are identified in `3DMOL-NOTICE.md`;
the original upstream release tag is unrecorded. This is an explicit provenance
limitation, not an inferred version. The existing project NOTICE continues to apply.

Public packages contain code, catalogue metadata, viewer assets/notices and
invented demonstrations. They exclude biological coordinate snapshots, OPM
frames, lipid donors, CHARMM parameter archives, real generated models, molecular
run bundles and logs. Third-party data is not relicensed as project code.
References, hashes and attribution do not authorize redistribution of source bytes.

`tools/build_membrane_release.py --out NEW_DIRECTORY` in the Workbench source
creates explicit code/asset snapshots, builds wheels and source archives without
network access or dependency installation, and screens every archive member.
`RELEASE_MANIFEST.json` records exact paths, hashes and findings. A zero-match
screen is a review aid, not publication approval or a complete legal assessment.

Run the synthetic bundle and Workbench membrane tests with the built package on
the Python import path. The HTTP suite needs permission to bind loopback ports.
The existing wheel namespace guard must pass. Preparation tests must not execute
new molecular jobs by default; reference-run replay needs its own compute budget.

## Scientific limits

Geometric orientation does not establish native sidedness, abundance, receptor
accessibility or curved lipid compatibility. A fixed-box restrained energy minimum
does not establish equilibrium, stability, biological accuracy or native neighbors.
Other microbes need reviewed chemistry before molecular preparation. Whole-organism
atomic reconstruction, general curved atomistic packing and equilibration are not
provided by this release. No external upload is required for any local workflow.
