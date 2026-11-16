"""Files for one domain-aware orientation run (see ``domain_orient``).

Given a local PDB file, this writes four plain files into one folder:

* ``DOMAIN_ORIENTATION.json`` - the method's result, the input file's SHA-256,
  the settings used, and what the result does and does not show.
* ``ORIENTED_STRUCTURE.pdb`` - the first model's atoms moved into the estimated
  membrane frame: membrane centre at the origin, membrane normal along +Z,
  +Z = the side the method calls outside (extracellular).
* ``MEMBRANE_LAYER.json`` - a small viewer descriptor (slab bounds and residue
  zone colours) in the same shape the Workbench viewer already reads.
* ``RESIDUE_ZONES.tsv`` - one row per C-alpha residue with its height (Z) and
  zone: membrane core, interface, positive (+Z) side or negative (-Z) side.

Software output only: the pose is an estimate from C-alpha geometry and
hydropathy, not a measured membrane position or a biological finding.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from . import domain_orient as _do

SCHEMA_VERSION = "1.0"
INTERFACE_BAND_A = 5.0
CLAIM = ("Estimated membrane pose from C-alpha geometry and hydropathy only; not a measured membrane position, "
         "and not evidence of native sidedness, exposure or function.")
ZONE_COLORS = {"membrane_core": "#d4a017", "interface": "#6d5bd0",
               "positive_side": "#2f8f5b", "negative_side": "#3a6fb0"}
SPANNING_METHODS = ("beta_barrel_core", "beta_barrel_core_open", "helical_core")


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _ca_records(path, chain=None):
    """(chain, resnum, resname) for exactly the residues ``domain_orient.read_pdb_ca`` keeps, in its order."""
    rows, seen = [], set()
    for ln in open(path):
        if ln.startswith("ENDMDL"):
            break
        if not ln.startswith(("ATOM", "HETATM")) or ln[12:16].strip() != "CA":
            continue
        if chain and ln[21] != chain:
            continue
        rn = ln[17:20]; key = (ln[21], ln[22:27])
        if rn not in _do.THREE or key in seen:
            continue
        if ln[16] not in " A":
            continue
        seen.add(key); rows.append((ln[21], int(ln[22:26]), rn, ln[26]))
    return rows


def membrane_frame(normal) -> np.ndarray:
    """Rotation whose rows are x, y, z with z = the membrane normal (right-handed)."""
    n = np.asarray(normal, float); n = n / np.linalg.norm(n)
    e1 = np.cross(n, [1.0, 0.0, 0.0] if abs(n[0]) < .9 else [0.0, 1.0, 0.0]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return np.array([e1, e2, n])


def zone_for(z: float, half: float) -> str:
    if abs(z) <= half:
        return "membrane_core"
    if abs(z) <= half + INTERFACE_BAND_A:
        return "interface"
    return "positive_side" if z > 0 else "negative_side"


def _oriented_pdb(path, rotation, origin) -> str:
    out = ["REMARK   1 MEMORIENT DOMAIN-AWARE ESTIMATE: MEMBRANE CENTRE AT ORIGIN, NORMAL ALONG +Z",
           "REMARK   1 +Z IS THE SIDE THE METHOD CALLS OUTSIDE; AN ESTIMATE, NOT A MEASUREMENT"]
    for ln in open(path):
        if ln.startswith("ENDMDL"):
            break
        if ln.startswith(("ATOM", "HETATM")) and len(ln) >= 54:
            xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            x, y, z = rotation @ (xyz - origin)
            out.append(f"{ln[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{ln[54:].rstrip()}")
        elif ln.startswith("TER"):
            out.append(ln.rstrip())
    out.append("END")
    return "\n".join(out) + "\n"


def run(pdb, out_dir, chain=None, half=None, skip_n=0) -> dict:
    """Orient ``pdb`` and write the four files into ``out_dir`` (created, must not exist yet)."""
    pdb = Path(pdb); out = Path(out_dir)
    ca, seq, num = _do.read_pdb_ca(str(pdb), chain)
    if len(ca) < 3:
        raise ValueError("No usable C-alpha atoms" + (f" in chain {chain!r}" if chain else ""))
    records = _ca_records(pdb, chain)
    if len(records) != len(ca):
        raise ValueError("Residue bookkeeping mismatch while reading C-alpha atoms")
    result = _do.orient(ca, seq, num, half=half, skip_n=skip_n)
    core_idx = result.pop("core_idx", None) or []
    out.mkdir(parents=True, exist_ok=False)
    found = result.get("method") not in (None, "none")
    document = {
        "schema_version": SCHEMA_VERSION, "module_id": "memorient.domain_orient",
        "input": {"file_name": pdb.name, "sha256": _sha256(pdb), "chain": chain or "", "c_alpha_residues": int(len(ca))},
        "parameters": {"chain": chain or "", "half": half, "skip_n": int(skip_n)},
        "result": result,
        "core_residues": [{"chain": records[i][0], "resnum": records[i][1]} for i in core_idx],
        "frame": {"axis": "+z", "membrane_center": [0.0, 0.0, 0.0], "positive_z": "outside (as estimated)"},
        "claim_ceiling": CLAIM, "scientific_qualification": "not_established",
        "missing_evidence": ["independent reference orientation not compared in this run",
                             "experimental topology or sidedness evidence not supplied"],
        "warnings": ([] if found else ["No membrane core (beta barrel or apolar transmembrane helices) was found; no pose was written."])
                    + (["Side assignment confidence is low."] if found and result.get("confidence") == "low" else []),
    }
    files = ["DOMAIN_ORIENTATION.json"]
    if found:
        R = membrane_frame(result["normal"]); origin = np.asarray(result["origin"], float); h = float(result["half"])
        (out / "ORIENTED_STRUCTURE.pdb").write_text(_oriented_pdb(pdb, R, origin), encoding="utf-8")
        z = (ca - origin) @ R[2]
        zones = [zone_for(float(v), h) for v in z]
        with open(out / "RESIDUE_ZONES.tsv", "w", encoding="utf-8") as fh:
            fh.write("chain\tresnum\tinsertion_code\tresname\tz_A\tzone\tmembrane_core_residue\n")
            core = set(core_idx)
            for i, (rec, v, zone) in enumerate(zip(records, z, zones)):
                fh.write(f"{rec[0]}\t{rec[1]}\t{rec[3].strip()}\t{rec[2]}\t{v:.2f}\t{zone}\t{'yes' if i in core else 'no'}\n")
        spanning = result["method"] in SPANNING_METHODS
        layer = {
            "context": "domain_aware_estimate", "method": result["method"],
            "label": "barrel" if result["method"].startswith("beta_barrel") else "helical" if result["method"] == "helical_core" else "anchored",
            "membrane_spanning": spanning,
            "orientation": {"axis": "+z", "membrane_center": [0.0, 0.0, 0.0],
                            "side_assignment": {"state": "estimated_" + result.get("confidence", "low") + "_confidence",
                                                "positive_compartment": "outside (estimated)",
                                                "negative_compartment": "inside (estimated)"}},
            "membrane_slab": {"half_thickness": h, "core_lower_z": -h, "core_upper_z": h, "side_a_z": h, "side_b_z": -h,
                              "asymmetric": False} if spanning else None,
            "residue_colors": [{"resid": rec[1], "chain": rec[0], "zone": zone, "accessibility": zone,
                                "extracellular": zone == "positive_side", "color": ZONE_COLORS[zone]}
                               for rec, zone in zip(records, zones)],
            "surface_set": [], "color_legend": ZONE_COLORS,
        }
        (out / "MEMBRANE_LAYER.json").write_text(json.dumps(layer, indent=1) + "\n", encoding="utf-8")
        document["zone_counts"] = {k: zones.count(k) for k in ZONE_COLORS}
        files += ["ORIENTED_STRUCTURE.pdb", "MEMBRANE_LAYER.json", "RESIDUE_ZONES.tsv"]
    document["outputs"] = files
    (out / "DOMAIN_ORIENTATION.json").write_text(json.dumps(document, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    return document
