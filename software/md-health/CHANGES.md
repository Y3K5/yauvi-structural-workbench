## 2026-10-13 — Supplied GROMACS inputs are audited and run bounded, and a short clean run is not evidence of anything

First release. This module takes GROMACS inputs that already exist, derives an
immutable copy with recorded checksums, runs a deliberately short pilot, and
reports concrete setup failures and declared-band departures in the logs and
energy series. It does not build a system: it will not choose or repair a force
field, assign protonation, add solvent, or override a preprocessing warning.

The distinction the module exists to hold is that **a short run that finishes
without a detected problem is not equilibrium, convergence, or biological
validation.** Every report states that, because the failure mode here is a green
pilot being read as a result about a protein.

What it does:

- `audit`/`derive` inspect supplied `.gro`, `.top`, `.mdp` and include files,
  copy them to a child directory with parent checksums preserved, and refuse an
  unsafe or missing reference rather than resolving it.
- The runner calls `grompp` and `mdrun` explicitly, with the chosen inputs and
  no implicit flags, and keeps every raw attempt including the failures.
- Log and `.xvg` parsing reports exit codes, fatal preprocessing errors, finite
  and non-finite sample counts, and departures from a declared band. Missing
  evidence is reported as unknown rather than inferred.
- A paired trajectory check reads frame counts and diagnostics from a completed
  run.

Verified against a real engine, on a non-biological case. GROMACS 2026.3 was
installed and a self-contained Lennard-Jones particle case was run three ways:

- the baseline completed 1000 steps through 1 ps with 101 finite potential and
  temperature samples, and its paired trajectory check read 101 frames;
- a topology whose molecule count disagreed with the coordinate count (64
  coordinates against 63 molecules) was **rejected at preprocessing**, with no
  `mdrun` reached;
- deliberately overlapping atoms **stopped dynamics at step 0** with a potential
  of 3.06282e+33, and the report names the observed failure rather than
  reporting a completed run.

One mistake made and fixed inside that work: the first energy export selected
terms in a way that omitted temperature, and the report correctly came back as
`insufficient_evidence` rather than reporting a temperature it did not have. The
corrected export uses an explicit named selection and is recorded separately, with
the original attempt retained.

What this does **not** establish: no qualified argon or protein property, no
constraint or LINCS behaviour, no convergence, no result on any independent host,
and nothing biological. A toy particle case tests the software's handling of a
real engine, not the physics of a protein.

53 tests use invented non-biological fixtures and mocked engine calls, so they
run without GROMACS present; three exercise the real command-line diagnostics
against clean, fatal and partial outcomes.
