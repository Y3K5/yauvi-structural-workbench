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


def read_structure(path: str | Path) -> tuple[list[dict], list[str]]:
    """Return (atoms, header_lines). Models are tracked so assemblies flatten."""
    text = Path(path).read_text(errors="replace").splitlines()
    if any(l.startswith("data_") for l in text[:5]):
        raise InputError("mmCIF is not supported yet; convert to PDB first")
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
    return atoms, header


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
    atoms, header = read_structure(structure)
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
        "counts": {
            "atoms_in": len(atoms), "atoms_kept": len(kept), "atoms_removed": len(removed),
            "removed_by_class": dict(Counter(
                classify(a["resname"]) if a["record"] == "HETATM" else "polymer"
                for a in removed)),
            "models_in": len({a["model"] for a in atoms}),
        },
    }


# ----------------------------------------------------------------- output
def _renumber_models(kept: list[dict], flatten: bool) -> tuple[list[str], dict]:
    """Assembly files reuse chain ids across MODEL records. Flatten to unique
    ids so a downstream tool reading chains alone sees the whole assembly and
    not just the first copy — the failure that reports a crystal-packing
    interface as if it were biological."""
    alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    remap, out = {}, []
    multi = len({a["model"] for a in kept}) > 1
    for a in kept:
        chain = a["chain"]
        if flatten and multi:
            key = (a["model"], a["chain"])
            if key not in remap:
                if len(remap) >= len(alpha):
                    raise InputError("assembly has more chains than available ids")
                remap[key] = alpha[len(remap)]
            chain = remap[key]
        out.append(a["line"][:21] + chain + a["line"][22:])
    return out, {f"model{m}:{c}": v for (m, c), v in remap.items()}


def write_outputs(result: dict, out_dir: str | Path) -> dict[str, str]:
    import json
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    lines, remap = _renumber_models(result["kept"], result["policy"]["flatten_models"])
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
        "schema_version": "1.0",
        "module_id": "structure_preparation",
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
        ],
    }
    (out / "PREPARATION.json").write_text(json.dumps(record, indent=1) + "\n")
    (out / "RUN_MANIFEST.json").write_text(json.dumps({
        "module_id": "structure_preparation",
        "inputs": {"structure": result["parent_path"],
                   "structure_sha256": result["parent_sha256"]},
        "outputs": {name: sha256(out / name) for name in
                    ("PREPARED.pdb", "REMOVED.tsv", "PREPARATION.json")},
        "policy": result["policy"],
    }, indent=1) + "\n")
    return {"prepared": str(prepared), "record": str(out / "PREPARATION.json")}
