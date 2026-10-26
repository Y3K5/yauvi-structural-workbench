"""``memorient`` command line.

Subcommands
-----------
* ``memorient contexts``            — list every membrane context (model / method / metrics)
* ``memorient describe``            — the module's machine-readable interface
* ``memorient describe <context>``  — full description of one context
* ``memorient validate <pdb>``      — check an input structure without orienting it
* ``memorient fetch``               — what raw files this module needs, and where from
* ``memorient orient <pdb>``        — orient a structure, print the per-residue label table,
                                      and write oriented PDB + 3Dmol JSON + PyMOL script

``contexts`` and ``describe`` are stdlib-only. ``orient`` imports the numeric layer
(numpy/biopython) lazily, so the first two work in a bare ``pip install .``.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path
import platform
import sys
from typing import List, Optional

from .contexts import get_context, list_contexts


def _load_topology_evidence(path: str | None, structure_path: str) -> dict | None:
    if not path:
        return None
    topology_path = Path(path)
    try:
        document = json.loads(topology_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read topology evidence: {exc}") from exc
    if not isinstance(document, dict):
        raise ValueError("topology evidence must be a JSON object")
    observed = hashlib.sha256(Path(structure_path).read_bytes()).hexdigest()
    if document.get("coordinate_sha256") != observed:
        raise ValueError("topology evidence coordinate checksum does not match the structure")
    source = document.get("source", {})
    if not isinstance(source, dict) or not str(source.get("id", "")).strip() or not str(source.get("citation", "")).strip():
        raise ValueError("topology evidence requires source.id and source.citation")
    return document


def _cmd_contexts(args: argparse.Namespace) -> int:
    ctxs = list_contexts()
    if args.json:
        print(json.dumps([c.to_dict() for c in ctxs], indent=2))
        return 0
    name_w = max(len(c.name) for c in ctxs)
    print(f"{'CONTEXT':<{name_w}}  {'MEMBRANE MODEL':<16}  {'ORIENTATION':<16}  METRICS")
    for c in ctxs:
        metrics = ", ".join(m for m in c.metrics)
        print(f"{c.name:<{name_w}}  {c.membrane_model:<16}  {c.orientation_method:<16}  {metrics}")
    return 0


def _cmd_describe(args: argparse.Namespace) -> int:
    # No context named: print the module's own interface. Every module in the
    # platform answers `describe` this way, so they can be discovered uniformly.
    # `memorient describe` with no argument was an error before, so this is
    # additive and the per-context form is unchanged.
    if not args.context:
        from .module_contract import describe as describe_module

        print(json.dumps(describe_module(), indent=2, sort_keys=True))
        return 0

    try:
        c = get_context(args.context)
    except KeyError as e:
        print(str(e), file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(c.to_dict(), indent=2))
        return 0
    tp = c.thickness_prior
    tp_str = f"{tp.mean:.1f} +/- {tp.sd:.1f} A (half-thickness)" if tp else "n/a (no bilayer)"
    print(f"# {c.name}")
    print(c.description)
    print()
    print(f"  membrane model      : {c.membrane_model}")
    print(f"  orientation method  : {c.orientation_method}")
    print(f"  thickness prior     : {tp_str}")
    print(f"  has membrane sides  : {c.has_membrane_sides}")
    print(f"  LPS shielding       : {c.lps_shielding}")
    print(f"  active metrics      : {', '.join(c.metrics)}")
    return 0


def _cmd_orient(args: argparse.Namespace) -> int:
    # Lazy import so `contexts`/`describe` work without numpy/biopython installed.
    try:
        from .geometry import load_structure
        from .orientor import orient_structure
        from .viz import display_oriented, write_3dmol_html, write_pymol_script
    except ImportError as e:  # pragma: no cover - dependency guard
        print(
            f"`memorient orient` needs the compute extras: pip install 'memorient[compute]'\n({e})",
            file=sys.stderr,
        )
        return 3

    try:
        ctx = get_context(args.context)
    except KeyError as e:
        print(str(e), file=sys.stderr)
        return 2

    try:
        topology = _load_topology_evidence(args.topology_evidence, args.pdb)
        structure = load_structure(args.pdb, chain=args.chain)
        result = orient_structure(structure, ctx, topology_evidence=topology)
    except (OSError, ValueError) as exc:
        print(f"memorient: {exc}", file=sys.stderr)
        return 2

    # Per-residue label table (compact). Full result is written to --out-json.
    rows = result.residue_table()
    if args.max_rows and len(rows) > args.max_rows:
        shown = rows[: args.max_rows]
        tail = f"  ... ({len(rows) - args.max_rows} more)"
    else:
        shown, tail = rows, ""
    hdr = [
        "resid", "insertion_code", "resname", "chain", "zone", "facing",
        "accessibility", "extracellular", "rsa", "sidedness_state", "biological_compartment", "legacy_accessibility_basis",
    ]
    print("\t".join(hdr))
    for r in shown:
        print("\t".join(str(r.get(k, "")) for k in hdr))
    if tail:
        print(tail)

    summary = result.summary()
    print("\n# summary", file=sys.stderr)
    for k, v in summary.items():
        print(f"  {k}: {v}", file=sys.stderr)

    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(result.to_dict(), fh, indent=2)
        print(f"wrote {args.out_json}", file=sys.stderr)
    if args.out_pdb:
        result.write_pdb(args.out_pdb)
        print(f"wrote {args.out_pdb}", file=sys.stderr)
    if args.out_viz:
        disp = display_oriented(result)
        with open(args.out_viz, "w") as fh:
            json.dump(disp, fh, indent=2)
        print(f"wrote {args.out_viz}", file=sys.stderr)
    if args.out_pymol:
        write_pymol_script(result, args.out_pymol)
        print(f"wrote {args.out_pymol}", file=sys.stderr)
    if args.out_html:
        write_3dmol_html(result, args.out_html)
        print(f"wrote {args.out_html}", file=sys.stderr)
    return 1 if result.scope_id == "alpha_helical" and result.scientific_state != "placement_evaluated" else 0


def _cmd_run(args: argparse.Namespace) -> int:
    """Common module-contract wrapper around the scientific ``orient`` command."""
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    forwarded = argparse.Namespace(
        pdb=args.structure,
        context=args.context,
        chain=args.chain,
        out_json=str(out / "MEMBRANE_ORIENTATION.json"),
        out_pdb=str(out / "ORIENTED_STRUCTURE.pdb"),
        out_viz=str(out / "MEMBRANE_LAYER.json"),
        out_pymol=str(out / "MEMBRANE_ORIENTATION.pml"),
        out_html=None,
        max_rows=0,
        topology_evidence=args.topology_evidence,
    )
    table = io.StringIO()
    with contextlib.redirect_stdout(table):
        code = _cmd_orient(forwarded)
    sys.stdout.write(table.getvalue())
    if code not in (0, 1):
        return code
    (out / "RESIDUE_ORIENTATION.tsv").write_text(table.getvalue(), encoding="utf-8")
    outputs = [
        "MEMBRANE_ORIENTATION.json",
        "ORIENTED_STRUCTURE.pdb",
        "MEMBRANE_LAYER.json",
        "MEMBRANE_ORIENTATION.pml",
        "RESIDUE_ORIENTATION.tsv",
    ]
    source = Path(args.structure)
    manifest = {
        "schema_version": "1.0",
        "module_id": "membrane_orientation",
        "version": "0.3.0",
        "input_sha256": {
            "structure": hashlib.sha256(source.read_bytes()).hexdigest(),
            **({"topology_evidence": hashlib.sha256(Path(args.topology_evidence).read_bytes()).hexdigest()}
               if args.topology_evidence else {}),
        },
        "parameters": {"context": args.context, "chain": args.chain or "all"},
        "runtime_versions": {"python": platform.python_version()},
        "outputs": outputs,
        "missing_evidence": [
            "external_opm_ppm_benchmark_adoption",
            *(["checksum_bound_topology_evidence"] if args.context in {"eukaryotic_pm", "tm_receptor"} and not args.topology_evidence else []),
        ],
        "limitations": [
            "Modeled orientation and coordinate accessibility do not prove native intact-cell exposure.",
            *(["Alpha-helical orientation is experimental and is not part of the Mark 1 qualified scope."]
              if args.context in {"eukaryotic_pm", "tm_receptor"} else []),
        ],
    }
    (out / "RUN_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return code


def _cmd_validate(args: argparse.Namespace) -> int:
    from .module_contract import validate_structure

    checks = validate_structure(args.pdb, chain=args.chain)
    failed = 0
    for name, ok, detail in checks:
        if ok:
            print(f"  [ok]      {name}" + (f"  {detail}" if detail else ""))
        else:
            failed += 1
            print(f"  [FAILED]  {name}  {detail}")
    return 0 if not failed else 1


def _cmd_fetch(args: argparse.Namespace) -> int:
    from .module_contract import SOURCES_MANIFEST

    if args.plan_only or not SOURCES_MANIFEST.is_file():
        if not SOURCES_MANIFEST.is_file():
            print(f"no source manifest shipped with this install", file=sys.stderr)
            return 2
        print(SOURCES_MANIFEST.read_text(encoding="utf-8"))
        return 0
    try:
        from yauvi_sources.cli import main as fetch_main
    except ImportError:
        print(
            "acquisition is provided by yauvi-sources, which is not installed.\n"
            "  pip install 'yauvi-sources[fetch]'\n"
            f"  then: yauvi-fetch plan --for memorient --manifest {SOURCES_MANIFEST}\n"
            "Or run `memorient fetch --plan` to print the declared sources.",
            file=sys.stderr,
        )
        return 1
    return fetch_main(["plan", "--for", "memorient", "--manifest", str(SOURCES_MANIFEST)])


def _cmd_surface(args: argparse.Namespace) -> int:
    """Build explicitly supplied geometry; never infer an organism membrane."""
    try:
        from .surfaces import Surface, place_receptor
        raw = Path(args.spec).read_bytes()
        spec = json.loads(raw)
        if spec.get("schema_version") != "memorient.surface.v1":
            raise ValueError("expected schema_version memorient.surface.v1")
        surface = Surface(**spec["surface"])
        result = {
            "schema_version": "memorient.surface_result.v1",
            "input_sha256": hashlib.sha256(raw).hexdigest(),
            "units": "angstrom", "angle_units": "radian",
            "scientific_state": "geometric_model_only",
            "specification": spec,
            "surface": surface.sample(**spec.get("sampling", {})),
            "receptors": [place_receptor(surface=surface, **item) for item in spec.get("receptors", [])],
            "limitations": ["No organism-specific reconstruction or lipid packing", "Geometric normal does not identify biological compartments", "No receptor-receptor or receptor-membrane clash checking"],
        }
        destination = Path(args.out)
        with destination.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    except (OSError, ValueError, KeyError, TypeError, ImportError) as exc:
        print(f"memorient surface: {exc}", file=sys.stderr)
        return 2
    return 0


def _cmd_viewer(args):
    from .surface_viewer import serve
    from .starter_catalogue import configure
    configure(args.cache)
    serve(args.port,manifest=args.manifest)
    return 0


def _cmd_place_pdb(args):
    from .pdb_surface import build_scene
    try:
        spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
        spec["pdb"] = Path(args.pdb).read_bytes().decode("utf-8")
        result = build_scene(spec)
        with Path(args.out).open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"memorient place-pdb: {exc}", file=sys.stderr)
        return 2
    return 0


def _cmd_lipid_plan(args):
    from .lipid_plan import plan_lipids
    try:
        result = plan_lipids(json.loads(Path(args.spec).read_text(encoding="utf-8")))
        with Path(args.out).open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"memorient lipid-plan: {exc}", file=sys.stderr)
        return 2
    return 0


def _cmd_build_ecoli_om(args):
    from .explicit_outer_membrane import save_patch
    try:
        receipt = save_patch(args.out)
        print(json.dumps({k:receipt[k] for k in ("state", "full_atoms", "view_atoms", "retained_lipids", "checks")}, indent=2))
    except (ValueError, OSError) as exc:
        print(f"memorient build-ecoli-om: {exc}", file=sys.stderr)
        return 2
    return 0


def _cmd_md_prepare(args):
    from .md_preparation import prepare_and_minimize
    try:
        result = prepare_and_minimize(args.profile, args.out, args.gmx, args.threads)
        print(json.dumps(result, indent=2))
        return 0 if result['state'] == 'minimized' else 3
    except (ValueError, OSError, KeyError, RuntimeError) as exc:
        print(f"memorient md-prepare: {exc}", file=sys.stderr)
        return 2


def _cmd_result_bundle(args):
    from .run_bundle import import_legacy,load_manifest,artifact_path
    import zipfile
    try:
        path=import_legacy(args.run,args.display,args.out)
        doc=load_manifest(path)
        with zipfile.ZipFile(Path(args.out)/'MEMBRANE_RESULT.zip','x',compression=zipfile.ZIP_DEFLATED) as archive:
            archive.write(path,'MEMBRANE_RUN.json')
            for ref in doc['artifacts'].values():archive.write(artifact_path(path.parent,ref['path']),ref['path'])
        print('Accepted local result bundle created. Contains molecular data; review before any sharing.')
        return 0
    except (OSError,ValueError,KeyError,TypeError) as exc:
        print(f'memorient result-bundle: {exc}',file=sys.stderr)
        return 2


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="memorient", description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="command", required=True)
    bundle=sub.add_parser('result-bundle',help='copy a verified legacy result to a portable local bundle; no simulation')
    for name in ('run','display','out'):bundle.add_argument('--'+name,required=True)
    bundle.set_defaults(func=_cmd_result_bundle)
    pmd = sub.add_parser("md-prepare", help="solvate and neutralize reviewed inputs, minimize locally; stop before equilibration")
    pmd.add_argument("--profile", required=True, help="source-bound local JSON profile")
    pmd.add_argument("--out", required=True, help="new output directory")
    pmd.add_argument("--gmx", default="gmx")
    pmd.add_argument("--threads", type=int, default=4)
    pmd.set_defaults(func=_cmd_md_prepare)
    pecoli = sub.add_parser("build-ecoli-om", help="assemble the pinned finite K-12 atomic membrane patch; no MD")
    pecoli.add_argument("--out", required=True, help="new output directory; refuses overwrite")
    pecoli.set_defaults(func=_cmd_build_ecoli_om)

    plipid = sub.add_parser("lipid-plan", help="check explicit-lipid recipe inputs; no packing")
    plipid.add_argument("--spec", required=True)
    plipid.add_argument("--out", required=True)
    plipid.set_defaults(func=_cmd_lipid_plan)

    pview = sub.add_parser("viewer", help="serve the private loopback membrane studio")
    pview.add_argument("--port", type=int, default=8956)
    pview.add_argument('--cache',help='explicit local checksum-addressed catalogue coordinate cache')
    pview.add_argument('--manifest',help='portable MEMBRANE_RUN.json to inspect')
    pview.set_defaults(func=_cmd_viewer)
    ppdb = sub.add_parser("place-pdb", help="place a local PDB on a supplied surface")
    ppdb.add_argument("--pdb", required=True)
    ppdb.add_argument("--spec", required=True, help="JSON surface and placement objects")
    ppdb.add_argument("--out", required=True, help="new scene JSON including placed PDB")
    ppdb.set_defaults(func=_cmd_place_pdb)

    ps = sub.add_parser("surface", help="experimental geometric membrane and receptor placement")
    ps.add_argument("--spec", required=True, help="local memorient.surface.v1 JSON")
    ps.add_argument("--out", required=True, help="new output JSON; refuses overwrite")
    ps.set_defaults(func=_cmd_surface)

    pc = sub.add_parser("contexts", help="list membrane contexts")
    pc.add_argument("--json", action="store_true", help="emit JSON")
    pc.set_defaults(func=_cmd_contexts)

    pd = sub.add_parser(
        "describe",
        help="the module interface, or one membrane context if named",
    )
    pd.add_argument(
        "context", nargs="?", default=None, help="context name (see `memorient contexts`)"
    )
    pd.add_argument("--json", action="store_true", help="emit JSON")
    pd.set_defaults(func=_cmd_describe)

    po = sub.add_parser("orient", help="orient a structure and label residues")
    po.add_argument("pdb", help="path to a PDB or mmCIF file")
    po.add_argument(
        "--context", "-c", default="gram_negative_om",
        help="membrane context (default: gram_negative_om)",
    )
    po.add_argument("--chain", default=None, help="restrict to one chain id")
    po.add_argument("--topology-evidence", default=None, help="checksum-bound transmembrane-span JSON")
    po.add_argument("--out-json", default=None, help="write full result JSON here")
    po.add_argument("--out-pdb", default=None, help="write oriented PDB here")
    po.add_argument("--out-viz", default=None, help="write 3Dmol display JSON here")
    po.add_argument("--out-pymol", default=None, help="write PyMOL .pml script here")
    po.add_argument("--out-html", default=None, help="write self-contained 3Dmol.js HTML viewer here")
    po.add_argument("--max-rows", type=int, default=40, help="max residue rows to print (0 = all)")
    po.set_defaults(func=_cmd_orient)

    pr = sub.add_parser("run", help="common-contract run with deterministic output names")
    pr.add_argument("--structure", required=True, help="path to a PDB or mmCIF file")
    pr.add_argument("--context", "-c", default="gram_negative_om")
    pr.add_argument("--chain", default=None)
    pr.add_argument("--topology-evidence", default=None)
    pr.add_argument("--out", required=True)
    pr.set_defaults(func=_cmd_run)

    pv = sub.add_parser("validate", help="check an input structure without orienting it")
    pv.add_argument("pdb", help="path to a PDB or mmCIF file")
    pv.add_argument("--chain", default=None, help="restrict to one chain id")
    pv.set_defaults(func=_cmd_validate)

    pf = sub.add_parser("fetch", help="what raw files this module needs, and where from")
    pf.add_argument(
        "--plan", dest="plan_only", action="store_true", help="print the declared sources only"
    )
    pf.set_defaults(func=_cmd_fetch)

    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
