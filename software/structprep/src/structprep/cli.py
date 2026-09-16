from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .components import CLASSES, DEFAULT_DROP, VERSION as TABLE_VERSION
from .core import InputError, default_policy, prepare, validate_policy, write_outputs


def describe() -> dict:
    return {
        "module_id": "structure_preparation", "component_id": "structprep",
        "package": "structprep", "version": "0.1.0",
        "display_name": "Structure Preparation",
        "one_line": "Which atoms should downstream analysis see, and what was removed to get there?",
        "commands": ["describe", "validate", "classes", "run"],
        "inputs": [
            {"name": "structure", "format": "PDB or mmCIF", "required": True},
            {"name": "policy", "format": "JSON", "required": False,
             "note": "declarative; defaults drop solvent, cryoprotectant and buffer"},
        ],
        "outputs": [
            {"name": "PREPARED.pdb or PREPARED.cif", "contract": "prepared_coordinates",
             "format": "pdb or mmcif",
             "note": "matches the parent format; atom records are re-emitted as written"},
            {"name": "PREPARATION.json", "contract": "preparation_record", "format": "json"},
            {"name": "REMOVED.tsv", "contract": "removed_atom_table", "format": "tsv"},
            {"name": "RUN_MANIFEST.json", "contract": "structure_analysis_run_manifest",
             "format": "json"},
        ],
        "deterministic_output": True,
        "component_table_version": TABLE_VERSION,
        "limitations": [
            "Removing a component does not establish that it was not functional.",
            "No protonation, charge assignment or minimisation is performed; the output "
            "is not docking-ready until a step declaring force field and pH has run.",
            "Coordinates are never altered, only selected and relabelled.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser("structprep", description="Structure preparation")
    ap.add_argument("--version", action="version", version="structprep 0.1.0")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("describe", help="print the machine-readable IO contract")
    sub.add_parser("classes", help="print the component classification table")

    for name in ("validate", "run"):
        s = sub.add_parser(name, help="check inputs without writing"
                           if name == "validate" else "prepare the structure")
        s.add_argument("--structure", required=True)
        s.add_argument("--policy", help="JSON policy file")
        s.add_argument("--keep", action="append", default=[],
                       help="residue name always retained. Repeatable.")
        s.add_argument("--drop", action="append", default=[],
                       help="residue name always removed, overriding contact protection.")
        s.add_argument("--chains", help="comma-separated chains to retain")
        s.add_argument("--altloc", help="alternate location to keep (default A)")
        s.add_argument("--no-flatten", action="store_true",
                       help="do not give assembly MODEL copies unique chain ids")
        s.add_argument("--allow-contacting", action="store_true",
                       help="remove classed components even when they touch the polymer")
        if name == "run":
            s.add_argument("--out", required=True)

    args = ap.parse_args(argv)

    if args.cmd == "describe":
        print(json.dumps(describe(), indent=1)); return 0
    if args.cmd == "classes":
        print(json.dumps({"version": TABLE_VERSION, "default_drop": list(DEFAULT_DROP),
                          "classes": {k: sorted(v) for k, v in CLASSES.items()}}, indent=1))
        return 0

    policy = default_policy()
    if args.policy:
        policy.update(json.loads(Path(args.policy).read_text()))
    if args.keep: policy["keep"] = list(policy["keep"]) + args.keep
    if args.drop: policy["drop"] = list(policy["drop"]) + args.drop
    if args.chains: policy["chains"] = [c.strip() for c in args.chains.split(",")]
    if args.altloc: policy["altloc"] = args.altloc
    if args.no_flatten: policy["flatten_models"] = False
    if args.allow_contacting: policy["keep_contacting"] = False

    try:
        policy = validate_policy(policy)
        result = prepare(args.structure, policy)
    except InputError as exc:
        print(f"structprep: {exc}", file=sys.stderr); return 2

    if args.cmd == "validate":
        print(json.dumps({"ok": True, "counts": result["counts"],
                          "refusals": len(result["refusals"]),
                          "warnings": result["warnings"]}, indent=1))
        return 0

    paths = write_outputs(result, args.out)
    print(json.dumps({"outputs": paths, "counts": result["counts"],
                      "refusals": len(result["refusals"]),
                      "warnings": result["warnings"]}, indent=1))
    # Non-zero when something was retained against the policy: the caller asked
    # for a removal that did not happen, and should know without reading a file.
    return 1 if result["refusals"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
