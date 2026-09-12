"""Deterministic structural editing with a derivation record.

structprep never edits in place. It writes new coordinates plus a record naming
the parent checksum and every atom removed, so `structqc` can be re-run on the
child and the provenance chain survives:

    structqc(raw) -> structprep -> structqc(prepared) -> downstream modules

No chemistry happens here. Protonation, tautomers, charges and minimisation are
force-field- and pH-dependent, and bundling them would let a model-dependent
assumption travel under the authority of a deterministic step.
"""
from __future__ import annotations

import hashlib
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from .components import DEFAULT_DROP, VERSION as TABLE_VERSION, classify


class InputError(RuntimeError):
    pass


class RefusedError(RuntimeError):
    """Raised when an edit would remove something that looks structural."""


CONTACT_CUTOFF_A = 4.0


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------- parsing
def _atom(line: str) -> dict[str, Any]:
    return {
        "record": line[:6].strip(),
        "name": line[12:16].strip(),
        "altloc": line[16].strip(),
        "resname": line[17:20].strip(),
        "chain": line[21],
        "resseq": line[22:26].strip(),
        "icode": line[26].strip(),
        "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
        "occupancy": float(line[54:60]) if line[54:60].strip() else 1.0,
        "line": line,
    }


CIF_NULLS = (".", "?")


def _cif_value(raw: str) -> str:
    return "" if raw in CIF_NULLS else raw


def _lines_with_offsets(text: str) -> list[tuple[int, str]]:
    out, pos = [], 0
    for line in text.splitlines(keepends=True):
        out.append((pos, line))
        pos += len(line)
    return out


def _locate_atom_site(text: str) -> tuple[list[str], int, int]:
    """Return (tags, data_start, data_end) for the one atom_site loop.

    Refuses rather than guesses. A file with no atom_site loop, or more than one,
    is not something this module will edit: picking the wrong loop would rewrite
    coordinates that were never examined.
    """
    lines = _lines_with_offsets(text)
    found = None
    i = 0
    while i < len(lines):
        if lines[i][1].strip() == "loop_":
            tags, j = [], i + 1
            while j < len(lines) and lines[j][1].lstrip().startswith("_"):
                tags.append(lines[j][1].strip().split()[0])
                j += 1
            if tags and tags[0].startswith("_atom_site."):
                if found is not None:
                    raise InputError(
                        "more than one _atom_site loop; refusing to guess which "
                        "one holds the coordinates")
                start = lines[j][0] if j < len(lines) else len(text)
                k = j
                while k < len(lines):
                    s = lines[k][1].lstrip()
                    if s.startswith(("loop_", "data_", "save_", "stop_", "#", "_")):
                        break
                    k += 1
                end = lines[k][0] if k < len(lines) else len(text)
                found = (tags, start, end)
            i = j
            continue
        i += 1
    if found is None:
        raise InputError(
            "no _atom_site loop found; this does not look like a coordinate mmCIF")
    return found


def _cif_tokens(text: str, start: int, end: int) -> list[tuple[str, int, int]]:
    """Tokenise a data region into (value, outer_start, outer_end).

    Offsets span the token exactly as written, quotes included, so a value can be
    spliced back into the source text without reformatting any of its neighbours.
    That is what lets this module edit an mmCIF without ever re-rendering a
    coordinate, the same guarantee the PDB path gives by slicing fixed columns.
    """
    toks: list[tuple[str, int, int]] = []
    i = start
    while i < end:
        c = text[i]
        if c in " \t\r\n":
            i += 1
            continue
        if c == "#":
            nl = text.find("\n", i, end)
            i = end if nl < 0 else nl + 1
            continue
        if c == ";" and (i == start or text[i - 1] == "\n"):
            raise InputError(
                "multi-line ';' value inside the atom_site loop; refusing rather "
                "than guess where the row ends")
        if c in "'\"":
            j = i + 1
            while j < end:
                if text[j] == c and (j + 1 >= end or text[j + 1] in " \t\r\n"):
                    break
                j += 1
            if j >= end:
                raise InputError("unterminated quoted value in the atom_site loop")
            toks.append((text[i + 1:j], i, j + 1))
            i = j + 1
            continue
        j = i
        while j < end and text[j] not in " \t\r\n":
            j += 1
        toks.append((text[i:j], i, j))
        i = j
    return toks


def _read_cif(path: str | Path) -> tuple[list[dict], list[str], dict]:
    """Parse an mmCIF atom_site loop, keeping each row's source text."""
    text = Path(path).read_text(errors="replace")
    tags, ds, de = _locate_atom_site(text)
    idx = {t.split(".", 1)[1]: k for k, t in enumerate(tags)}
    toks = _cif_tokens(text, ds, de)
    ncol = len(tags)
    if ncol == 0 or len(toks) % ncol:
        raise InputError(
            f"atom_site loop holds {len(toks)} values across {ncol} columns; the "
            "row count is not whole, so the file will not be edited")
    for needed in ("Cartn_x", "Cartn_y", "Cartn_z"):
        if needed not in idx:
            raise InputError(f"atom_site loop has no _atom_site.{needed}")

    atoms: list[dict] = []
    for r in range(0, len(toks), ncol):
        row = toks[r:r + ncol]
        row_start, row_end = row[0][1], row[-1][2]

        def get(*keys: str, default: str = "") -> str:
            for key in keys:
                k = idx.get(key)
                if k is not None:
                    value = _cif_value(row[k][0])
                    if value:
                        return value
            return default

        try:
            x, y, z = (float(row[idx[c]][0]) for c in ("Cartn_x", "Cartn_y", "Cartn_z"))
        except ValueError as exc:
            raise InputError(
                f"non-numeric coordinate in atom_site row {r // ncol + 1}") from exc
        try:
            occupancy = float(get("occupancy") or 1.0)
        except ValueError:
            occupancy = 1.0
        try:
            model = int(float(get("pdbx_PDB_model_num") or 0))
        except ValueError:
            model = 0

        # Emit the whole physical line when the row occupies one, so trailing
        # whitespace and column alignment survive exactly as deposited. A row
        # wrapped across lines falls back to its token span.
        if "\n" in text[row_start:row_end]:
            emit_start, emit_end = row_start, row_end
        else:
            emit_start = text.rfind("\n", 0, row_start) + 1
            nl = text.find("\n", row_end)
            emit_end = len(text) if nl < 0 else nl

        # Only auth_asym_id is relabelled on flattening — it is the mmCIF analogue
        # of the PDB chain column. label_asym_id is an entity label, not a chain
        # name, and the two legitimately differ: in 8DBK the ATP of auth chain A
        # carries label_asym_id G. Overwriting it would rewrite the entity
        # assignment while claiming only to have relabelled a chain.
        auth = idx.get("auth_asym_id")
        auth_span = ((row[auth][1] - emit_start, row[auth][2] - emit_start)
                     if auth is not None else None)

        atoms.append({
            "record": (get("group_PDB", default="ATOM")).upper(),
            "name": get("auth_atom_id", "label_atom_id"),
            "altloc": get("label_alt_id", "auth_alt_id"),
            "resname": get("auth_comp_id", "label_comp_id"),
            "chain": get("auth_asym_id", "label_asym_id"),
            "resseq": get("auth_seq_id", "label_seq_id"),
            "icode": get("pdbx_PDB_ins_code"),
            "x": x, "y": y, "z": z,
            "occupancy": occupancy,
            "model": model,
            "line": text[emit_start:emit_end],
            "_auth_span": auth_span,
        })

    if not atoms:
        raise InputError(f"no atom_site rows found in {path}")
    return atoms, [], {"format": "mmcif", "text": text,
                       "data_start": ds, "data_end": de}


def _read_any(path: str | Path) -> tuple[list[dict], list[str], dict]:
    head = Path(path).read_text(errors="replace")[:4096].splitlines()
    if (str(path).lower().endswith((".cif", ".mmcif"))
            or any(l.startswith("data_") for l in head[:50])):
        return _read_cif(path)
    return _read_pdb(path)


def read_structure(path: str | Path) -> tuple[list[dict], list[str]]:
    """Return (atoms, header_lines), for callers that predate mmCIF support."""
    atoms, header, _ = _read_any(path)
    return atoms, header


def _read_pdb(path: str | Path) -> tuple[list[dict], list[str], dict]:
    """Return (atoms, header_lines, layout). Models are tracked so assemblies flatten."""
    text = Path(path).read_text(errors="replace").splitlines()
    atoms, header, model = [], [], 0
    for line in text:
        if line.startswith("MODEL"):
            try:
                model = int(line.split()[1])
            except (IndexError, ValueError):
                model += 1
            continue
        if line.startswith(("ENDMDL", "END")):
            continue
        if line.startswith(("ATOM", "HETATM")):
            a = _atom(line)
            a["model"] = model
            atoms.append(a)
        elif line.startswith(("HEADER", "TITLE", "CRYST1", "REMARK")):
            header.append(line)
    if not atoms:
        raise InputError(f"no ATOM or HETATM records found in {path}")
    return atoms, header, {"format": "pdb"}


# ----------------------------------------------------------------- geometry
def _contacts_polymer(hetatms: list[dict], polymer: list[dict], cutoff: float) -> set[tuple]:
    """Which het components sit within `cutoff` of any polymer atom.

    Brute force is fine at this scale and keeps the result exactly reproducible;
    a spatial index would introduce a tolerance nobody declared.
    """
    if not hetatms or not polymer:
        return set()
    near = set()
    c2 = cutoff * cutoff
    for h in hetatms:
        key = (h["model"], h["chain"], h["resname"], h["resseq"], h["icode"])
        if key in near:
            continue
        for p in polymer:
            dx = h["x"] - p["x"]
            if dx * dx > c2:
                continue
            dy, dz = h["y"] - p["y"], h["z"] - p["z"]
            if dx * dx + dy * dy + dz * dz <= c2:
                near.add(key)
                break
    return near


def _bridging(hetatms: list[dict], polymer: list[dict], cutoff: float) -> set[tuple]:
    """Het components touching two or more polymer chains — interface waters."""
    touched = defaultdict(set)
    c2 = cutoff * cutoff
    for h in hetatms:
        key = (h["model"], h["chain"], h["resname"], h["resseq"], h["icode"])
        for p in polymer:
            dx = h["x"] - p["x"]
            if dx * dx > c2:
                continue
            dy, dz = h["y"] - p["y"], h["z"] - p["z"]
            if dx * dx + dy * dy + dz * dz <= c2:
                touched[key].add(p["chain"])
    return {k for k, chains in touched.items() if len(chains) > 1}


# ----------------------------------------------------------------- policy
def default_policy() -> dict[str, Any]:
    return {
        "drop_classes": list(DEFAULT_DROP),
        "keep": [],              # residue names always retained
        "drop": [],              # residue names always removed, overrides keep_contacting
        "chains": None,          # None = all
        "flatten_models": True,  # assembly files reuse chain ids across MODELs
        "altloc": "A",           # which alternate location to retain
        "min_occupancy": 0.0,
        "keep_contacting": True, # refuse to drop het within CONTACT_CUTOFF_A of polymer
        "contact_cutoff_A": CONTACT_CUTOFF_A,
    }


def validate_policy(policy: dict) -> dict:
    p = {**default_policy(), **(policy or {})}
    from .components import CLASSES
    bad = set(p["drop_classes"]) - set(CLASSES)
    if bad:
        raise InputError(f"unknown component classes: {sorted(bad)}; "
                         f"known: {sorted(CLASSES)}")
    if p["altloc"] and len(str(p["altloc"])) != 1:
        raise InputError("altloc must be a single character")
    return p


# ----------------------------------------------------------------- preparation
def prepare(structure: str | Path, policy: dict | None = None) -> dict[str, Any]:
    p = validate_policy(policy or {})
    atoms, header, layout = _read_any(structure)
    parent = sha256(structure)

    polymer = [a for a in atoms if a["record"] == "ATOM"]
    het = [a for a in atoms if a["record"] == "HETATM"]

    contacting = (_contacts_polymer(het, polymer, p["contact_cutoff_A"])
                  if p["keep_contacting"] else set())
    bridging = _bridging(het, polymer, p["contact_cutoff_A"])

    keep_names = {n.upper() for n in p["keep"]}
    drop_names = {n.upper() for n in p["drop"]}
    drop_classes = set(p["drop_classes"])

    kept, removed, refusals, warnings = [], [], [], []
    for a in atoms:
        key = (a["model"], a["chain"], a["resname"], a["resseq"], a["icode"])
        cls = classify(a["resname"]) if a["record"] == "HETATM" else "polymer"
        reason = None

        if p["chains"] and a["chain"] not in p["chains"]:
            reason = "chain not selected"
        elif p["altloc"] and a["altloc"] and a["altloc"].upper() != p["altloc"].upper():
            reason = f"altloc {a['altloc']} (kept {p['altloc']})"
        elif a["occupancy"] < p["min_occupancy"]:
            reason = f"occupancy {a['occupancy']:.2f} below {p['min_occupancy']}"
        elif a["record"] == "HETATM":
            name = a["resname"].upper()
            if name in drop_names:
                reason = "named in drop list"
            elif name in keep_names:
                reason = None
            elif cls in drop_classes:
                # The refusal that matters: a component in contact with the
                # polymer may be structural — a catalytic haem, a bound metal —
                # and residue code alone cannot tell you. Solvent is exempt:
                # essentially every crystal water touches the protein, so the
                # rule would refuse everything and mean nothing. A water that
                # bridges two chains still earns a warning below.
                if (key in contacting and p["keep_contacting"]
                        and cls != "solvent"):
                    refusals.append({
                        "component": name, "chain": a["chain"], "resseq": a["resseq"],
                        "class": cls,
                        "reason": f"within {p['contact_cutoff_A']} A of polymer; "
                                  f"retained. Name it in `drop` to remove it.",
                    })
                else:
                    reason = f"class {cls}"
                    if key in bridging:
                        warnings.append(
                            f"{name} {a['chain']}{a['resseq']} bridges two chains "
                            f"and was removed; interface geometry may change")

        (removed if reason else kept).append(
            {**a, "removed_because": reason} if reason else a)

    if p["altloc"]:
        alts = {a["altloc"] for a in atoms if a["altloc"]}
        for other in sorted(alts - {p["altloc"].upper()}):
            warnings.append(f"alternate location {other} existed and was discarded")

    return {
        "parent_sha256": parent, "parent_path": str(structure),
        "policy": p, "component_table_version": TABLE_VERSION,
        "kept": kept, "removed": removed,
        "refusals": [dict(t) for t in {tuple(sorted(r.items())) for r in refusals}],
        "warnings": sorted(set(warnings)),
        "header": header,
        "layout": layout,
        "counts": {
            "atoms_in": len(atoms), "atoms_kept": len(kept), "atoms_removed": len(removed),
            "removed_by_class": dict(Counter(
                classify(a["resname"]) if a["record"] == "HETATM" else "polymer"
                for a in removed)),
            "models_in": len({a["model"] for a in atoms}),
        },
    }


# ----------------------------------------------------------------- output
def _chain_ids(*, multichar: bool):
    """Chain labels for flattened assemblies.

    PDB has one column for the chain id, so the single-character alphabet is a
    hard ceiling there — 62 copies. mmCIF has no such limit, so the generator
    continues into two-character ids rather than refusing a large assembly.
    """
    alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    yield from alpha
    if not multichar:
        return
    for first in alpha:
        for second in alpha:
            yield first + second


def _renumber_models(kept: list[dict], flatten: bool,
                     layout: dict | None = None) -> tuple[list[str], dict]:
    """Assembly files reuse chain ids across MODEL records. Flatten to unique
    ids so a downstream tool reading chains alone sees the whole assembly and
    not just the first copy — the failure that reports a crystal-packing
    interface as if it were biological."""
    fmt = (layout or {}).get("format", "pdb")
    remap, out = {}, []
    multi = len({a["model"] for a in kept}) > 1
    ids = _chain_ids(multichar=fmt == "mmcif")
    for a in kept:
        chain = a["chain"]
        if flatten and multi:
            key = (a["model"], a["chain"])
            if key not in remap:
                try:
                    remap[key] = next(ids)
                except StopIteration:
                    raise InputError(
                        "assembly has more chains than available ids") from None
            chain = remap[key]
        if fmt == "mmcif":
            span = a.get("_auth_span")
            if chain == a["chain"] or span is None:
                out.append(a["line"])          # untouched, byte for byte
            else:
                start, end = span
                out.append(a["line"][:start] + chain + a["line"][end:])
        else:
            out.append(a["line"][:21] + chain + a["line"][22:])
    return out, {f"model{m}:{c}": v for (m, c), v in remap.items()}


def write_outputs(result: dict, out_dir: str | Path) -> dict[str, str]:
    import json
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    layout = result.get("layout") or {"format": "pdb"}
    lines, remap = _renumber_models(
        result["kept"], result["policy"]["flatten_models"], layout)
    if layout["format"] == "mmcif":
        # Everything outside the atom_site rows — every other category, the header
        # blocks, the trailing categories — is carried through byte for byte.
        source = layout["text"]
        prepared = out / "PREPARED.cif"
        prepared.write_text(source[:layout["data_start"]] + "\n".join(lines)
                            + "\n" + source[layout["data_end"]:])
    else:
        prepared = out / "PREPARED.pdb"
        prepared.write_text("\n".join(result["header"] + lines) + "\nEND\n")

    rows = ["\t".join(("record", "resname", "chain", "resseq", "atom", "class", "reason"))]
    for a in result["removed"]:
        rows.append("\t".join((
            a["record"], a["resname"], a["chain"], a["resseq"], a["name"],
            classify(a["resname"]) if a["record"] == "HETATM" else "polymer",
            a["removed_because"] or "")))
    (out / "REMOVED.tsv").write_text("\n".join(rows) + "\n")

    record = {
        "schema_version": "1.1",
        "module_id": "structure_preparation",
        "coordinate_format": layout["format"],
        "derived_from": {
            "path": result["parent_path"],
            "sha256": result["parent_sha256"],
        },
        "prepared_sha256": sha256(prepared),
        "component_table_version": result["component_table_version"],
        "policy": result["policy"],
        "counts": result["counts"],
        "chain_remap": remap,
        "refusals": result["refusals"],
        "warnings": result["warnings"],
        "limitations": [
            "Removing a component does not establish that it was not functional.",
            "No chemistry was applied: no protonation, tautomer assignment, charge "
            "assignment or minimisation. This structure is not docking-ready until a "
            "step that declares its force field and pH has run.",
            "Coordinates are unchanged; only selection and chain labelling differ.",
            "Atom records are re-emitted as written in the parent file; no coordinate "
            "was parsed and re-rendered, so no value can drift through this step.",
            "On mmCIF input only _atom_site.auth_asym_id is relabelled when an assembly "
            "is flattened. label_asym_id is left as deposited, so entity labels and "
            "chain names can differ in the prepared file exactly as they did in the parent.",
        ],
    }
    (out / "PREPARATION.json").write_text(json.dumps(record, indent=1) + "\n")
    (out / "RUN_MANIFEST.json").write_text(json.dumps({
        "module_id": "structure_preparation",
        "inputs": {"structure": result["parent_path"],
                   "structure_sha256": result["parent_sha256"]},
        "outputs": {name: sha256(out / name) for name in
                    (prepared.name, "REMOVED.tsv", "PREPARATION.json")},
        "policy": result["policy"],
    }, indent=1) + "\n")
    return {"prepared": str(prepared), "record": str(out / "PREPARATION.json")}
