"""The module's declared interface, and the checks behind it.

Every module in the platform answers the same four commands — `run`, `describe`,
`validate`, `fetch` — so that they can be discovered, checked, and composed
without special-casing each one. This file holds memorient's half of that
contract.

`describe` output is what `catalogs/modules/membrane_orientation.yaml` is checked
against: a module that quietly stops emitting what it declares should fail a
test, not surprise someone reading the results.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

MODULE_ID = "membrane_orientation"

SOURCES_MANIFEST = Path(__file__).with_name("sources.yaml")


def describe() -> dict[str, Any]:
    """The machine-readable interface. Built from the live context registry."""
    from .contexts import list_contexts

    version = _version()
    return {
        # Two names, deliberately both reported. `module_id` is how the platform
        # addresses this module (catalogs/modules/membrane_orientation.yaml);
        # `package` is how pip and `yauvi-fetch --for` address it. Reporting only
        # one of them makes the other look like a typo.
        "module_id": MODULE_ID,
        "package": "memorient",
        "display_name": "Membrane Orientation and Accessibility",
        "version": version,
        "one_line": (
            "How does this structure sit in its membrane, and which residues face out?"
        ),
        "commands": ["run", "orient", "contexts", "describe", "validate", "fetch"],
        "inputs": [
            {
                "name": "structure",
                "format": "PDB or mmCIF",
                "required": True,
                "note": "one structure per invocation; chain selectable with --chain",
            },
            {
                "name": "topology_evidence",
                "format": "JSON",
                "required": False,
                "note": (
                    "coordinate-checksum-bound transmembrane spans; required for the "
                    "experimental tm_helix_axis_v2 method"
                ),
            },
        ],
        # These are the CLI's own per-structure outputs. They deliberately claim
        # no platform contract id: the contracts the descriptor names
        # (`oriented_antigens`, `orientation_manifest`) are run-level artifacts
        # produced by the pipeline stage over a whole candidate set, not by a
        # single `memorient orient` invocation. Claiming them here would assert a
        # shape this command does not produce.
        "outputs": [
            {"name": "MEMBRANE_ORIENTATION.json", "format": "json", "note": "full orientation result from `run`"},
            {"name": "RESIDUE_ORIENTATION.tsv", "format": "tsv", "note": "per-residue labels from `run`"},
            {"name": "RUN_MANIFEST.json", "format": "json", "note": "runtime, inputs, limitations, and output inventory"},
            {"name": "--out-json", "format": "json", "note": "full orientation result for one structure"},
            {"name": "--out-pdb", "format": "pdb", "note": "the structure, oriented in the bilayer frame"},
            {"name": "--out-viz", "format": "json", "note": "3Dmol display payload"},
            {"name": "--out-pymol", "format": "text", "note": "PyMOL script"},
            {"name": "--out-html", "format": "text", "note": "self-contained viewer"},
        ],
        "platform_contracts": {
            "note": (
                "produced by the pipeline stage over a candidate set, not by this CLI"
            ),
            "outputs": ["oriented_antigens", "orientation_manifest"],
        },
        "contexts": [c.name for c in list_contexts()],
        "runtimes": ["yauvi-python"],
        "optional_runtimes": [],
        "extras": {
            "compute": "numpy, scipy, biopython — required by `orient`",
            "examples": "adds matplotlib for the benchmark scripts",
        },
        "reproducible": True,
        "limitations": [
            "Orientation is geometric. It fits a bilayer to coordinates; it does not "
            "predict topology from sequence, and it cannot correct a structure that is "
            "wrong.",
            "A predicted monomer is oriented as supplied. Whether that monomer is the "
            "form the protein takes in the membrane is a separate question, answered by "
            "fold_state, not here.",
            "Solvent accessibility is computed on the coordinates given. It is not "
            "antibody accessibility on an intact cell.",
            "The thickness prior is per-context and is a prior, not a measurement; two "
            "runs under different contexts are not comparable on absolute z-coordinates.",
            "Mark 1 qualification is limited to beta-barrel membrane proteins. "
            "Alpha-helical orientation remains experimental.",
        ],
        "scientific_scopes": [
            {
                "scope_id": "beta_barrel",
                "scientific_state": "conditionally_qualified",
                "benchmark_collection": "yauvi-structural-public-qualification-v2",
                "release_blocking": True,
                "supported_subject_class": "beta-barrel membrane proteins",
            },
            {
                "scope_id": "alpha_helical",
                "scientific_state": "prototype",
                "benchmark_collection": "yauvi-structural-public-qualification-v2",
                "release_blocking": False,
                "supported_subject_class": "alpha-helical membrane proteins",
            },
        ],
    }


def _version() -> str:
    from . import __version__

    return __version__


def validate_structure(path: str | Path, *, chain: str | None = None) -> list[tuple[str, bool, str]]:
    """Check an input structure without orienting it.

    Returns (check_name, ok, detail) triples so the caller decides how to report.
    """
    checks: list[tuple[str, bool, str]] = []
    path = Path(path)

    if not path.is_file():
        return [("input:exists", False, f"no such file: {path}")]
    checks.append(("input:exists", True, str(path)))

    try:
        from .geometry import load_structure
    except ImportError as exc:
        checks.append(
            (
                "extras:compute",
                False,
                f"the compute extras are not installed ({exc}); "
                f"pip install 'memorient[compute]'",
            )
        )
        return checks
    checks.append(("extras:compute", True, "numpy, scipy and biopython are importable"))

    try:
        structure = load_structure(str(path), chain=chain)
    except Exception as exc:  # noqa: BLE001 - any parse failure is a failed check
        checks.append(("input:parses", False, f"{type(exc).__name__}: {exc}"))
        return checks
    checks.append(("input:parses", True, "coordinates were read"))

    count = len(getattr(structure, "residues", ()) or ())
    if count == 0:
        checks.append(("input:residues", False, "no residues were read from the file"))
    else:
        checks.append(("input:residues", True, f"{count} residue(s)"))

    return checks
