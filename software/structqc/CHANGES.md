## 2026-09-11 — a reference spread across six chains is not full coverage

Diagnostic change in `analyze`, recorded here because it changes what a manifest means
to a reader, not what any number is.

A single-copy reference aligned against several concatenated chains spreads across
them: reference position 1 lands in one chain, position 2 in the next.
`coverage_fraction` counts distinct mapped reference positions, so the concatenation
still scored **1.0** with no warning. On PDB 8DBK, a six-chain cryo-EM hexamer, a run
with no `--chain` reports `coverage_fraction: 1.0`, `mapped_residues: 318` and
`coordinate_residues: 1868` — a manifest that looks authoritative and whose per-chain
residue map is unusable. The same run with `--chain B` reports 0.993711 over 316
residues at `identity_fraction` 1.0.

That silence produced a wrong result downstream. A 30-entry survey of catalytic-loop
occupancy came back with an implausibly tidy "exactly one chain per entry resolves the
loop" pattern, which was this line distributing one reference across six chains. It was
caught only because the regularity was too neat to believe.

`analyze` now appends a warning whenever `chain` is None and more than one polymer
chain is present, naming the chains and saying what `coverage_fraction` actually
counted.

**Why a warning and not a refusal.** Refusing is the stronger fix and remains the right
end state. One adopted qualification record — the declared clean multichain control —
runs with no chain selected and carries a frozen `coverage_fraction`, and the panel's
immutability policy forbids editing a threshold in place: changing that number requires
a new collection version with a recorded revision. Making the ambiguity visible changes
no number and can ship now; the refusal should follow with the collection bump.

No output field changed, no threshold moved, and `warnings` was already part of the
contract. Tests unchanged: 8 passing.
