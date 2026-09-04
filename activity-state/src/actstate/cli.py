"""`actstate` — activity-state classification from raw files.

    actstate run --in <dir-or-files> --out <dir>
    actstate describe                       machine-readable IO contract
    actstate validate --in <dir>            check inputs without running
    actstate fetch --plan                   what raw files are needed, and where from

`run` takes an annotation table, optionally a FASTA and a directory of
structures, and writes a result document plus a flat table. It needs no
workspace, no project, and no campaign.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Sequence

from . import __version__
from .core import SITE_CLUSTER_MAX_ANGSTROM, LABELS, SIGNAL_STATES, assess
from .io import (
    InputError,
    attach_sequences,
    build_document,
    find_structure,
    read_annotation_table,
    read_fasta,
    read_sidecar,
    write_json,
    write_table,
)
from .structure import StructureError, read_structure

EXIT_OK = 0
EXIT_BLOCKED = 1
EXIT_USAGE = 2

RESULT_NAME = "ACTIVITY_STATE.json"
TABLE_NAME = "ACTIVITY_STATE.tsv"

# Conventional file names inside an input directory, so `--in <dir>` works
# without naming each file.
DEFAULT_ANNOTATION_NAMES = ("annotations.tsv", "uniprot.tsv", "annotation.tsv")
DEFAULT_FASTA_NAMES = ("sequences.fasta", "proteins.fasta", "sequences.faa")
DEFAULT_STRUCTURE_DIRS = ("structures", "pdb")


def _describe() -> dict:
    """The module's declared interface. Kept in step with catalogs/modules/."""
    return {
        # Two names, deliberately both reported. `module_id` is how the platform
        # addresses this module (catalogs/modules/activity_state.yaml); `package`
        # is how pip and `yauvi-fetch --for` address it.
        "module_id": "activity_state",
        "package": "actstate",
        "display_name": "Activity State Classification",
        "version": __version__,
        "one_line": (
            "Does the available evidence support this protein being in a functionally "
            "competent state?"
        ),
        "commands": ["run", "describe", "validate", "fetch"],
        "inputs": [
            {
                "name": "annotation_table",
                "format": "TSV/CSV with UniProt columns",
                "required": True,
                "columns_used": [
                    "accession", "sequence", "ft_act_site", "ft_binding",
                    "ft_site", "cc_cofactor", "ec", "xref_interpro", "xref_pfam",
                ],
            },
            {"name": "fasta", "format": "FASTA", "required": False},
            {"name": "structures", "format": "directory of PDB/mmCIF", "required": False},
            {"name": "fold_state", "format": "JSON keyed by accession", "required": False},
            {
                "name": "reference_comparison",
                "format": "JSON keyed by accession",
                "required": False,
            },
            {
                "name": "expected_residues",
                "format": "JSON keyed by accession, then by catalytic position",
                "required": False,
                "note": (
                    "The residue an experimentally validated reference carries at each "
                    "annotated position. Required to reach active_site_disrupted."
                ),
            },
        ],
        "outputs": [
            {"name": RESULT_NAME, "contract": "activity_state_summary", "format": "JSON"},
            {"name": TABLE_NAME, "contract": "activity_state_table", "format": "TSV"},
        ],
        "labels": list(LABELS),
        "signals": [
            {"name": "completeness", "needs": "annotation + sequence"},
            {"name": "geometry", "needs": "annotation + structure"},
            {"name": "occupancy", "needs": "annotation + structure"},
            {"name": "conformation", "needs": "structural aligner + curated references"},
            {"name": "assembly", "needs": "fold_state output"},
        ],
        "signal_states": list(SIGNAL_STATES),
        "runtimes": ["yauvi-python"],
        "optional_runtimes": ["foldseek"],
        "reproducible": True,
        "deterministic_output": True,
        "limitations": [
            "Signals are reported separately and must never be merged into a single "
            "activity score.",
            "An entry with no ACT_SITE annotation is indeterminate, never inactive; "
            "absence of annotation is not evidence of absence of function.",
            "A predicted model can never on its own yield active_state_supported.",
            "Active-site clustering is evaluated on a representative side-chain atom "
            "per residue, not on full side-chain rotamers; it detects a dispersed site, "
            "not a subtly misaligned one.",
            "The conformation signal requires a curated set of references of known "
            "state. No such set ships with this module, so the signal reports "
            "unavailable until one is supplied.",
            "Catalytic competence is judged by residue identity, which does not "
            "detect a site disabled by a change outside the annotated positions.",
            "active_site_disrupted requires an expected residue for the position from "
            "an experimentally validated reference. Without one, a residue outside the "
            "competence set is reported as contradicting evidence and caps the label at "
            "indeterminate; membership in a broad residue set is not a position-specific "
            "chemistry test.",
        ],
    }


def cmd_describe(args) -> int:
    print(json.dumps(_describe(), indent=2, sort_keys=True))
    return EXIT_OK


def _first_existing(directory: Path, names: Sequence[str]) -> Path | None:
    return next((directory / n for n in names if (directory / n).is_file()), None)


def _resolve_inputs(args) -> dict:
    """Work out which raw files to use, from --in plus any explicit overrides."""
    resolved: dict = {"annotation": None, "fasta": None, "structures": None}
    base = Path(args.input) if args.input else None

    if base and base.is_dir():
        resolved["annotation"] = _first_existing(base, DEFAULT_ANNOTATION_NAMES)
        resolved["fasta"] = _first_existing(base, DEFAULT_FASTA_NAMES)
        for name in DEFAULT_STRUCTURE_DIRS:
            if (base / name).is_dir():
                resolved["structures"] = base / name
                break
    elif base and base.is_file():
        resolved["annotation"] = base
    elif base:
        raise InputError(f"--in path does not exist: {base}")

    if args.annotation:
        resolved["annotation"] = Path(args.annotation)
    if args.fasta:
        resolved["fasta"] = Path(args.fasta)
    if args.structures:
        resolved["structures"] = Path(args.structures)

    if resolved["annotation"] is None:
        raise InputError(
            "no annotation table found. Pass --annotation <file>, or --in <dir> "
            f"containing one of: {', '.join(DEFAULT_ANNOTATION_NAMES)}"
        )
    return resolved


def cmd_validate(args) -> int:
    inputs = _resolve_inputs(args)
    records = read_annotation_table(inputs["annotation"])
    print(f"annotation table: {inputs['annotation']}  ({len(records)} protein(s))")

    with_sites = [r for r in records if r.features().has_catalytic_annotation]
    print(f"  entries with an annotated catalytic site: {len(with_sites)}")

    if inputs["fasta"]:
        sequences = read_fasta(inputs["fasta"])
        records = attach_sequences(records, sequences)
        print(f"fasta: {inputs['fasta']}  ({len(sequences)} sequence(s))")
    missing_sequence = [r.accession for r in records if not r.sequence]
    if missing_sequence:
        print(f"  entries with no sequence: {len(missing_sequence)}")

    if inputs["structures"]:
        found = sum(1 for r in records if find_structure(inputs["structures"], r.accession))
        print(f"structures: {inputs['structures']}  ({found}/{len(records)} matched)")
    else:
        print("structures: none supplied — geometry and occupancy will be unevaluated")

    unparsed = sum(len(r.features().unparsed) for r in records)
    if unparsed:
        print(f"  feature positions that could not be parsed: {unparsed}")

    if not with_sites:
        print()
        print(
            "No entry carries an ACT_SITE annotation. Every protein would be reported "
            "as indeterminate, which is correct but uninformative — check that the "
            "annotation export includes the ft_act_site field."
        )
        return EXIT_BLOCKED
    return EXIT_OK


def cmd_run(args) -> int:
    inputs = _resolve_inputs(args)
    out_dir = Path(args.output)

    records = read_annotation_table(inputs["annotation"])
    if inputs["fasta"]:
        records = attach_sequences(records, read_fasta(inputs["fasta"]))

    fold_states = read_sidecar(args.fold_state)
    comparisons = read_sidecar(args.reference_comparison)
    expected_residues = read_sidecar(args.expected_residues)

    structure_errors: list[str] = []
    assessments = []
    for record in records:
        structure = None
        if inputs["structures"]:
            path = find_structure(inputs["structures"], record.accession)
            if path is not None:
                try:
                    structure = read_structure(path, identifier=record.accession)
                except StructureError as exc:
                    # A structure that cannot be read is recorded as unread. It
                    # never silently becomes "no structure supplied".
                    structure_errors.append(f"{record.accession}: {exc}")
        assessments.append(
            assess(
                record,
                structure=structure,
                chain=args.chain,
                reference_comparison=comparisons.get(record.accession),
                fold_state=fold_states.get(record.accession),
                expected_residues=expected_residues.get(record.accession),
                max_separation=args.max_separation,
            )
        )

    config = {
        "max_separation_angstrom": args.max_separation,
        "chain": args.chain or "all",
        "fasta_supplied": bool(inputs["fasta"]),
        "structures_supplied": bool(inputs["structures"]),
        "fold_state_supplied": bool(fold_states),
        "reference_comparison_supplied": bool(comparisons),
        "expected_residues_supplied": bool(expected_residues),
        # A mistyped expectation does nothing, and doing nothing silently is how
        # a curator concludes their entry was applied. Named here and printed.
        "rejected_expectations": sorted(
            f"{a.accession}: {reason}"
            for a in assessments
            for reason in (a.signal("completeness").values.get("rejected_expectations", ())
                           if a.signal("completeness") else ())
        ),
        "unreadable_structures": sorted(structure_errors),
    }
    document = build_document(assessments, config=config)
    result_path = write_json(out_dir / RESULT_NAME, document)
    table_path = write_table(out_dir / TABLE_NAME, assessments)

    print(f"{len(assessments)} protein(s) assessed")
    for label, count in document["summary"]["labels"].items():
        print(f"  {label:<24} {count}")
    if structure_errors:
        print(f"\n{len(structure_errors)} structure(s) could not be read:")
        for message in structure_errors:
            print(f"  - {message}")
    if config["rejected_expectations"]:
        print(f"\n{len(config['rejected_expectations'])} expected-residue entr(ies) were "
              "rejected and had no effect:")
        for message in config["rejected_expectations"]:
            print(f"  - {message}")
    print(f"\nwrote {result_path}")
    print(f"wrote {table_path}")
    return EXIT_OK


def cmd_fetch(args) -> int:
    """Defer to yauvi-fetch, which owns acquisition, and say so if it is absent."""
    manifest = Path(__file__).with_name("sources.yaml")
    if args.plan_only:
        print(manifest.read_text(encoding="utf-8"))
        return EXIT_OK
    try:
        from yauvi_sources.cli import main as fetch_main  # noqa: PLC0415
    except (ImportError, ModuleNotFoundError):
        print(
            "acquisition is provided by yauvi-sources, which is not installed.\n"
            "  pip install 'yauvi-sources[fetch]'\n"
            f"  then: yauvi-fetch plan --for actstate --manifest {manifest}\n"
            "Or run `actstate fetch --plan` to print the declared sources.",
            file=sys.stderr,
        )
        return EXIT_BLOCKED
    return fetch_main(["plan", "--for", "actstate", "--manifest", str(manifest)])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="actstate",
        description="Classify whether evidence supports a protein being in a working state.",
    )
    parser.add_argument("--version", action="version", version=f"actstate {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_input_args(p):
        p.add_argument("--in", dest="input", help="input directory, or an annotation table")
        p.add_argument("--annotation", help="UniProt-style annotation TSV/CSV")
        p.add_argument("--fasta", help="sequences, if not in the annotation table")
        p.add_argument("--structures", help="directory of PDB/mmCIF files")

    p_run = sub.add_parser("run", help="assess every protein and write the results")
    add_input_args(p_run)
    p_run.add_argument("--out", dest="output", required=True, help="output directory")
    p_run.add_argument("--chain", help="restrict geometry to one chain")
    p_run.add_argument("--fold-state", help="JSON of fold_state records, keyed by accession")
    p_run.add_argument(
        "--expected-residues",
        help=(
            "JSON of expected catalytic residues, keyed by accession then position. "
            "Required to reach active_site_disrupted."
        ),
    )
    p_run.add_argument(
        "--reference-comparison",
        help="JSON of reference-state comparisons, keyed by accession",
    )
    p_run.add_argument(
        "--max-separation",
        type=float,
        default=SITE_CLUSTER_MAX_ANGSTROM,
        help=f"active-site cluster bound in angstrom (default {SITE_CLUSTER_MAX_ANGSTROM})",
    )
    p_run.set_defaults(func=cmd_run)

    p_validate = sub.add_parser("validate", help="check inputs without running")
    add_input_args(p_validate)
    p_validate.set_defaults(func=cmd_validate)

    p_describe = sub.add_parser("describe", help="print the machine-readable IO contract")
    p_describe.set_defaults(func=cmd_describe)

    p_fetch = sub.add_parser("fetch", help="plan the raw files this module needs")
    p_fetch.add_argument(
        "--plan", dest="plan_only", action="store_true", help="print the declared sources only"
    )
    p_fetch.set_defaults(func=cmd_fetch)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (InputError, StructureError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
