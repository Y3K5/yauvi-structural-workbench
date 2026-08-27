# YAUVI Mark 1 Qualification v2

This directory freezes the scientific scopes, panel strata, split counts, evidence requirements, and numerical gates for Mark 1.

It is intentionally fail-closed. `PANEL_MANIFEST.json` currently records the complete required panel but contains no adopted v2 cases. Existing v1 public files are candidate material only. They do not become v2 evidence until a curator records the exact source release, checksum, license, citation, split, expected result, mapping, and exclusion rationale in a new immutable manifest version.

Run the offline collection audit with:

```text
python benchmarks/qualification-v2/run_qualification.py
```

The runner performs no acquisition and no network requests. It validates the source lock, panel composition, development/held-out separation, required metadata, and release-scope policy. It writes deterministic JSON, TSV, HTML, and checksum artifacts under `results/`.

Current expected result: `blocked_panel_incomplete`. This is not a software failure and it is not a scientific qualification pass.

Scope boundaries:

- Beta-barrel MembraneOrient is the Mark 1 release-blocking membrane scope.
- Alpha-helical orientation remains experimental and non-blocking.
- StateAtlas Mark 1 scope is limited to ABL-family resemblance over exact UniProt ABL1 residues 242–495 mappings.
- Resemblance is not activity, orientation is not native exposure, a mapped site is not catalysis, an interface is not affinity, and similarity is not exact functional transfer.
