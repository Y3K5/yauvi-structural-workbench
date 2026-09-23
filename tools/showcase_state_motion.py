"""Derived, inspectable endpoint geometry for the ABL showcase; no trajectory claim."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.vectors import calc_angle, calc_dihedral


def fit(moving, fixed):
    """Proper row-vector Kabsch transform, with reflections forbidden."""
    mc, fc = moving.mean(axis=0), fixed.mean(axis=0)
    u, _, vt = np.linalg.svd((moving - mc).T @ (fixed - fc))
    if np.linalg.det(u @ vt) < 0:
        u[:, -1] *= -1
    rotation = u @ vt
    return rotation, fc - mc @ rotation


def angles(residues, position):
    r = residues[position]
    out = {"phi": None, "psi": None, "chi1": None, "n_ca_c": None}
    def measure(names, key, angle=False):
        try:
            atoms = [residues[p][a] for p, a in names]
            if any(a.is_disordered() for a in atoms):
                return
            value = (calc_angle if angle else calc_dihedral)(*[a.get_vector() for a in atoms])
            if math.isfinite(value):
                out[key] = round(math.degrees(value), 2)
        except KeyError:
            pass
    # Never calculate a torsion across an absent residue or broken peptide bond.
    if position - 1 in residues and "C" in residues[position - 1] and "N" in r:
        if 1.0 < np.linalg.norm(residues[position - 1]["C"].coord - r["N"].coord) < 1.9:
            measure([(position - 1, "C"), (position, "N"), (position, "CA"), (position, "C")], "phi")
    if position + 1 in residues and "C" in r and "N" in residues[position + 1]:
        if 1.0 < np.linalg.norm(r["C"].coord - residues[position + 1]["N"].coord) < 1.9:
            measure([(position, "N"), (position, "CA"), (position, "C"), (position + 1, "N")], "psi")
    # Restrict chi1 to DFG residues with the unambiguous N-CA-CB-CG definition.
    if r.resname in {"ASP", "PHE"}:
        measure([(position, a) for a in ["N", "CA", "CB", "CG"]], "chi1")
    measure([(position, a) for a in ["N", "CA", "C"]], "n_ca_c", angle=True)
    return out


def pdb_text(residues, rotation, translation):
    lines, serial, previous = [], 0, None
    for pos, res in sorted(residues.items()):
        if previous is not None and pos != previous + 1:
            lines.append("TER")
        for atom in res:
            if atom.is_disordered():
                continue
            serial += 1
            x, y, z = atom.coord.astype(float) @ rotation + translation
            lines.append(f"ATOM  {serial:5d} {atom.name:>4} {res.resname:>3} A{pos:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{atom.bfactor:6.2f}          {atom.element:>2}")
        previous = pos
    return "\n".join(lines + ["TER", "END", ""])


def build_motion(panel: Path):
    structures, provenance, exclusions = [], [], []
    for entry, state in [("2GQG", "active"), ("2HYY", "inactive")]:
        source = panel / f"sources/wwpdb/{entry}.cif"
        mapping = panel / f"references/abl/alignment_map.{entry}.json"
        run = panel / f"results/execution-abl/v2-stateatlas-{state}-held_out-{entry}/RUN_MANIFEST.json"
        hashes = json.loads(run.read_text())["input_sha256"]
        for path, key in [(source, "structure"), (mapping, "alignment_map")]:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if digest != hashes[key]:
                raise ValueError(f"State motion input drift: {path.relative_to(panel)}")
            provenance.append({"path": str(path.relative_to(panel)), "sha256": digest})
        doc = json.loads(mapping.read_text())
        if doc["coordinate_system"] != "uniprot" or doc["domain"] != {"uniprot_start": 242, "uniprot_end": 495}:
            raise ValueError("Unexpected state motion coordinate system or domain")
        chain = MMCIFParser(QUIET=True).get_structure(entry, source)[0]["A"]
        mapped = {}
        for row in doc["query"]:
            pos = row["uniprot_position"]
            if row["mapping_state"] != "exact" or row.get("chain_id") != "A":
                continue
            key = (" ", int(row["auth_seq_id"]), row.get("insertion_code") or " ")
            if key in chain and "CA" in chain[key] and not chain[key]["CA"].is_disordered():
                mapped[pos] = chain[key]
        structures.append(mapped)
    a, b = structures
    positions = sorted(p for p in a.keys() & b.keys() if a[p].resname == b[p].resname)
    exclusions = sorted(set(range(242, 496)) - set(positions))
    if len(positions) / 254 < .9:
        raise ValueError("State motion shared exact mapping covers less than 90% of domain")
    a = {p: a[p] for p in positions}; b = {p: b[p] for p in positions}
    active = np.array([a[p]["CA"].coord for p in positions], dtype=float)
    inactive = np.array([b[p]["CA"].coord for p in positions], dtype=float)
    rotation, translation = fit(inactive, active)
    placed = inactive @ rotation + translation
    displacement = np.linalg.norm(active - placed, axis=1)
    residues = []
    for i, pos in enumerate(positions):
        residues.append({"position": pos, "resname": a[pos].resname,
                         "active": active[i].round(5).tolist(), "inactive": placed[i].round(5).tolist(),
                         "displacement": round(float(displacement[i]), 4),
                         "active_angles": angles(a, pos), "inactive_angles": angles(b, pos)})
    regions = [{"id": "activation", "label": "Activation loop / DFG · 381–402", "start": 381, "end": 402},
               {"id": "dfg", "label": "DFG motif · Asp381–Phe382–Gly383", "start": 381, "end": 383},
               {"id": "alpha", "label": "αC region · 280–299", "start": 280, "end": 299},
               {"id": "ploop", "label": "P-loop region · 248–255", "start": 248, "end": 255},
               {"id": "all", "label": "Whole mapped kinase domain · 242–495", "start": 242, "end": 495}]
    for region in regions:
        selected = [r["displacement"] for r in residues if region["start"] <= r["position"] <= region["end"]]
        region["mean_displacement"] = round(float(np.mean(selected)), 4)
        region["mapped_count"] = len(selected)
    return {"schema_version": 1, "active_entry": "2GQG", "inactive_entry": "2HYY", "chain": "A",
            "coordinate_system": "UniProt P00519; exact frozen query maps; first coordinate model",
            "method": "Unweighted CA Kabsch fit of inactive onto active over all shared exact same-residue positions within 242–495; no flexible-region exclusion or outlier rejection.",
            "claim": "New derived endpoint geometry for teaching; not a StateAtlas classification rerun, molecular-dynamics trajectory or transition pathway.",
            "source_inputs": provenance, "mapped_count": len(positions), "excluded_positions": exclusions,
            "exclusion_note": "Missing, nonstandard, ambiguous or sequence-different positions are omitted. Tyr393 is phosphorylated (PTR) in 2GQG and is not bridged by the motion trace. Disordered atoms are excluded from angle measurements and exported models.",
            "rmsd": round(float(np.sqrt(np.mean(displacement ** 2))), 4),
            "rotation": rotation.tolist(), "translation": translation.tolist(),
            "regions": regions, "residues": residues,
            "active_pdb": pdb_text(a, np.eye(3), np.zeros(3)),
            "inactive_pdb": pdb_text(b, rotation, translation)}
