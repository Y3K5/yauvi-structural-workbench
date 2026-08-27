"""Reading raw inputs and writing results.

Inputs are the files a user actually has: a UniProt annotation TSV, a FASTA, and
a directory of structures. Nothing here needs a workspace, a project, or a
campaign — which is what makes the module runnable on its own.

Outputs are byte-deterministic. The result document carries no timestamp and no
absolute path, and every mapping is written with sorted keys, so two runs over
the same inputs produce identical bytes and a digest is a meaningful thing to
record. Run metadata that *is* time-varying goes in a separate file.
"""
from __future__ import annotations

from dataclasses import asdict
import csv
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .core import ActivityAssessment, ProteinRecord

SCHEMA_VERSION = "1.0"

# Column aliases, so an export made with UniProt's own field names or with its
# display headers both work without the user having to rename anything.
COLUMN_ALIASES: Mapping[str, tuple[str, ...]] = {
    "accession": ("accession", "Entry", "entry", "protein_id", "id"),
    "sequence": ("sequence", "Sequence"),
    "act_site_raw": ("ft_act_site", "Active site", "active_site"),
    "binding_raw": ("ft_binding", "Binding site", "binding_site"),
    "site_raw": ("ft_site", "Site", "site"),
    "cofactor_raw": ("cc_cofactor", "Cofactor", "cofactor"),
    "ec_number": ("ec", "EC number", "ec_number"),
    "interpro": ("xref_interpro", "InterPro", "interpro"),
    "pfam": ("xref_pfam", "Pfam", "pfam"),
    "protein_name": ("protein_name", "Protein names", "Protein name"),
}


class InputError(RuntimeError):
    pass


def _resolve_columns(header: Sequence[str]) -> dict[str, str]:
    """Map our field names onto whatever the file calls them."""
    present = {name.strip(): name for name in header}
    resolved: dict[str, str] = {}
    for field, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in present:
                resolved[field] = present[alias]
                break
    if "accession" not in resolved:
        raise InputError(
            "annotation table has no recognisable accession column; expected one of: "
            + ", ".join(COLUMN_ALIASES["accession"])
        )
    return resolved


def read_annotation_table(path: str | Path) -> list[ProteinRecord]:
    """Read a UniProt-style TSV (or CSV) into records."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"annotation table not found: {path}")
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise InputError(f"annotation table is empty: {path}")
        columns = _resolve_columns(reader.fieldnames)
        records = []
        for row in reader:
            accession = (row.get(columns["accession"]) or "").strip()
            if not accession:
                continue
            records.append(
                ProteinRecord(
                    accession=accession,
                    **{
                        field: (row.get(column) or "").strip()
                        for field, column in columns.items()
                        if field != "accession"
                    },
                )
            )
    if not records:
        raise InputError(f"annotation table declares no proteins: {path}")
    return records


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read sequences, keyed by accession.

    Handles both a bare accession header and the UniProt `sp|ACC|NAME` form,
    which is what `subproteo` carries proteins by.
    """
    path = Path(path)
    if not path.is_file():
        raise InputError(f"FASTA not found: {path}")

    sequences: dict[str, list[str]] = {}
    key: str | None = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            header = line[1:].strip()
            parts = header.split("|")
            key = parts[1] if len(parts) >= 3 else (header.split() or [""])[0]
            sequences.setdefault(key, [])
        elif key is not None:
            sequences[key].append(line.strip())
    if not sequences:
        raise InputError(f"no sequences found in {path}")
    return {key: "".join(chunks) for key, chunks in sequences.items()}


def attach_sequences(
    records: Iterable[ProteinRecord], sequences: Mapping[str, str]
) -> list[ProteinRecord]:
    """Fill in sequences from a FASTA for records that lack one."""
    out = []
    for record in records:
        if record.sequence or record.accession not in sequences:
            out.append(record)
            continue
        out.append(
            ProteinRecord(**{**asdict(record), "sequence": sequences[record.accession]})
        )
    return out


def find_structure(directory: str | Path, accession: str) -> Path | None:
    """Locate a structure for an accession, by the conventional file names."""
    directory = Path(directory)
    if not directory.is_dir():
        return None
    candidates = (
        f"AF-{accession}-F1.pdb",
        f"AF-{accession}-F1-model_v4.pdb",
        f"{accession}.pdb",
        f"{accession}.cif",
        f"{accession}.mmcif",
    )
    for name in candidates:
        path = directory / name
        if path.is_file():
            return path
    matches = sorted(directory.glob(f"*{accession}*"))
    return next((p for p in matches if p.suffix.lower() in (".pdb", ".cif", ".mmcif")), None)


# -- results --------------------------------------------------------------


def _signal_payload(signal) -> dict:
    return {
        "name": signal.name,
        "state": signal.state,
        "detail": signal.detail,
        "values": dict(signal.values),
    }


def assessment_payload(assessment: ActivityAssessment) -> dict:
    return {
        "accession": assessment.accession,
        "label": assessment.label,
        "rationale": assessment.rationale,
        "catalytic_positions": list(assessment.catalytic_positions),
        "declared_cofactors": list(assessment.declared_cofactors),
        "unparsed_features": list(assessment.unparsed_features),
        "signals": [_signal_payload(s) for s in assessment.signals],
    }


def build_document(
    assessments: Sequence[ActivityAssessment], *, config: Mapping[str, object]
) -> dict:
    """The result document. Deterministic: no timestamps, no absolute paths."""
    counts: dict[str, int] = {}
    for assessment in assessments:
        counts[assessment.label] = counts.get(assessment.label, 0) + 1
    return {
        "schema_version": SCHEMA_VERSION,
        "module_id": "activity_state",
        "config": dict(config),
        "summary": {
            "proteins": len(assessments),
            "labels": dict(sorted(counts.items())),
        },
        "assessments": [
            assessment_payload(a) for a in sorted(assessments, key=lambda a: a.accession)
        ],
    }


def write_json(path: str | Path, document: Mapping[str, object]) -> Path:
    """Write canonically: sorted keys, fixed separators, trailing newline."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(document, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")
    return path


def write_table(path: str | Path, assessments: Sequence[ActivityAssessment]) -> Path:
    """A flat TSV of one row per protein, for reading alongside other channels."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    signal_names = ("completeness", "geometry", "occupancy", "conformation", "assembly")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["accession", "label", *(f"{name}_state" for name in signal_names),
             "catalytic_positions", "declared_cofactors"]
        )
        for assessment in sorted(assessments, key=lambda a: a.accession):
            states = []
            for name in signal_names:
                signal = assessment.signal(name)
                states.append(signal.state if signal else "unavailable")
            writer.writerow(
                [
                    assessment.accession,
                    assessment.label,
                    *states,
                    ";".join(str(p) for p in assessment.catalytic_positions),
                    ";".join(assessment.declared_cofactors),
                ]
            )
    return path


def read_sidecar(path: str | Path | None) -> dict[str, Mapping[str, object]]:
    """Read an optional per-accession JSON map (fold_state, or reference comparisons)."""
    if not path:
        return {}
    path = Path(path)
    if not path.is_file():
        raise InputError(f"sidecar file not found: {path}")
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except ValueError as exc:
        raise InputError(f"sidecar file is not valid JSON: {path}: {exc}") from exc
    if not isinstance(document, Mapping):
        raise InputError(f"sidecar must be a JSON object keyed by accession: {path}")
    return {str(k): v for k, v in document.items() if isinstance(v, Mapping)}
