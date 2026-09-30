from __future__ import annotations
import argparse
import json
import sys


def describe():
    return {"module_id": "bio_orient", "package": "yauvi-bio-orient", "version": "0.1.0.dev0",
            "display_name": "Bio-Orient", "commands": ["describe", "validate", "run"],
            "inputs": ["structure", "patch_declaration", "chemical_reference (optional)", "topology_evidence (optional)", "supporting_evidence (optional)"],
            "outputs": ["BIO_ORIENT.json", "BIO_ORIENT_SIDES.tsv", "BIO_ORIENT_EDGES.tsv", "BIO_ORIENT_LAYER.json", "RUN_MANIFEST.json"],
            "scientific_readiness": "experimental_nonblocking", "claim_ceiling": "Geometry and recorded evidence only; functional availability remains separately supported or unknown."}


def main(argv=None):
    parser = argparse.ArgumentParser(prog="bio-orient")
    sub = parser.add_subparsers(dest="command", required=True)
    description = sub.add_parser("describe")
    description.add_argument("--structure", help="include exact scope and a patch declaration starter")
    description.add_argument("--model", type=int, default=0)
    description.add_argument("--assembly-id", default="asu")
    for name in ("validate", "run"):
        p = sub.add_parser(name)
        p.add_argument("--structure", required=True)
        p.add_argument("--patch-declaration", required=True)
        p.add_argument("--model", type=int, default=0)
        p.add_argument("--assembly-id", default="asu")
        p.add_argument("--chemical-reference")
        p.add_argument("--topology-evidence")
        p.add_argument("--evidence", action="append", default=[])
        if name == "run":
            p.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "describe":
            document = describe()
            if args.structure:
                from structqc.coordinate_scope import load_scope
                from .core import binding
                scope = load_scope(args.structure, args.model, args.assembly_id)
                document["binding"] = binding(scope)
                document["components"] = scope["components"]
                document["patch_declaration_starter"] = {"schema_version": "1.0", "binding": binding(scope), "sides": []}
            print(json.dumps(document, indent=2))
            return 0
        from .core import analyze, write_outputs
        result = analyze(args.structure, args.patch_declaration, model_index=args.model, assembly_id=args.assembly_id,
                         chemical_reference=args.chemical_reference, topology_path=args.topology_evidence, evidence_paths=args.evidence)
        if args.command == "run":
            write_outputs(result, args.out)
        print(json.dumps({"valid": True, "sides": len(result["sides"]), "edges": len(result["edges"]), "chirality": result["chirality"]["state"], "scientific_readiness": "experimental_nonblocking"}))
        return 0
    except (ValueError, OSError, KeyError, TypeError, ImportError) as exc:
        print(f"bio-orient: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
