"""memorient.geometry — structure IO, rigid frames, and PCA frame canonicalization.

The one data object that flows through the whole pipeline is :class:`Structure`: a light
container of per-residue arrays (CA coordinates, residue identity, side-chain centroid
vectors, per-residue confidence) plus the raw atom records needed to write an oriented PDB
back out. Everything downstream (SASA, membrane fit, labeling, viz) reads a ``Structure``.

Frame canonicalization
----------------------
A predicted structure arrives in an arbitrary coordinate frame. Before any orientation
metric is computed, :func:`canonicalize` rotates the structure into a *deterministic* PCA
frame so every downstream decision depends only on intrinsic geometry, not on how the file
happened to be saved. The construction guarantees a **proper rotation** (det = +1, never a
reflection):

* principal axes = eigenvectors of the CA covariance, ordered by descending eigenvalue;
* the sign of the first two axes is fixed by the sign of the coordinate's third moment
  (skew) along that axis — a frame-intrinsic tiebreak;
* the third axis is replaced by ``cross(axis0, axis1)`` so the basis is right-handed.

This is what makes ``rotation_validate`` reach Jaccard 1.0: re-orienting the same structure
from any input rotation yields the same canonical frame (up to the degenerate case of equal
eigenvalues, which is reported).
"""

from __future__ import annotations

import io
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Residues we treat as standard amino acids (three-letter). Non-standard/hetero are dropped
# from the per-residue arrays but kept in the atom records so the PDB still round-trips.
_AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE",  # selenomethionine, treated as MET
}
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


@dataclass
class Structure:
    """Per-residue view of a protein chain (or set of chains) for orientation.

    Arrays are parallel and indexed by residue order in the file:

    * ``ca`` (N,3)          — C-alpha coordinates (Angstrom)
    * ``resids`` (N,)       — author residue numbers (int)
    * ``icodes`` (N,)       — author insertion codes (empty string when absent)
    * ``resnames`` (N,)     — three-letter residue names
    * ``chains`` (N,)       — chain ids
    * ``sc_vec`` (N,3)      — unit vector CA -> side-chain centroid (0 for GLY / missing)
    * ``plddt`` (N,)        — per-residue confidence (CA B-factor; pLDDT for AF models), or NaN
    * ``atoms``             — list of raw atom dicts for PDB serialization (all atoms, all residues)
    """

    ca: np.ndarray
    resids: np.ndarray
    resnames: np.ndarray
    chains: np.ndarray
    sc_vec: np.ndarray
    plddt: np.ndarray
    atoms: List[dict] = field(default_factory=list)
    source: str = ""
    icodes: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        self.ca = np.asarray(self.ca, dtype=float).reshape(-1, 3)
        n = len(self.ca)
        self.resids = np.asarray(self.resids).reshape(n)
        self.resnames = np.asarray(self.resnames, dtype=object).reshape(n)
        self.chains = np.asarray(self.chains, dtype=object).reshape(n)
        self.sc_vec = np.asarray(self.sc_vec, dtype=float).reshape(n, 3)
        self.plddt = np.asarray(self.plddt, dtype=float).reshape(n)
        if self.icodes is None:
            self.icodes = np.full(n, "", dtype=object)
        else:
            self.icodes = np.asarray(self.icodes, dtype=object).reshape(n)

    def __len__(self) -> int:
        return len(self.ca)

    @property
    def sequence(self) -> str:
        return "".join(_THREE_TO_ONE.get(str(r), "X") for r in self.resnames)

    def copy(self) -> "Structure":
        return Structure(
            ca=self.ca.copy(),
            resids=self.resids.copy(),
            resnames=self.resnames.copy(),
            chains=self.chains.copy(),
            sc_vec=self.sc_vec.copy(),
            plddt=self.plddt.copy(),
            atoms=[dict(a) for a in self.atoms],
            source=self.source,
            icodes=self.icodes.copy(),
        )

    def transformed(self, R: np.ndarray, t: Optional[np.ndarray] = None) -> "Structure":
        """Return a copy with ``x -> R @ x + t`` applied to CA, atoms, and side-chain vectors.

        ``R`` must be a proper rotation (3x3). ``t`` defaults to zero. Side-chain *vectors*
        rotate but do not translate.
        """
        R = np.asarray(R, dtype=float).reshape(3, 3)
        t = np.zeros(3) if t is None else np.asarray(t, dtype=float).reshape(3)
        out = self.copy()
        out.ca = self.ca @ R.T + t
        # rotate side-chain unit vectors (no translation)
        out.sc_vec = self.sc_vec @ R.T
        for a in out.atoms:
            xyz = np.array([a["x"], a["y"], a["z"]]) @ R.T + t
            a["x"], a["y"], a["z"] = float(xyz[0]), float(xyz[1]), float(xyz[2])
        return out


# --------------------------------------------------------------------------------------
# Rigid-frame helpers
# --------------------------------------------------------------------------------------


def rotation_matrix_to_z(v: Sequence[float]) -> np.ndarray:
    """Proper rotation mapping unit vector ``v`` onto +Z (``R @ v = [0,0,1]``).

    Uses Rodrigues' formula; handles the antiparallel case (v == -Z) by a 180 deg flip about
    X. Returns a right-handed rotation matrix (det = +1).
    """
    v = np.asarray(v, dtype=float).reshape(3)
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("cannot rotate a zero vector")
    v = v / n
    z = np.array([0.0, 0.0, 1.0])
    c = float(np.dot(v, z))
    if c > 1.0 - 1e-12:
        return np.eye(3)
    if c < -1.0 + 1e-12:
        # antiparallel: 180 deg about any axis perpendicular to z, e.g. X
        return np.diag([1.0, -1.0, -1.0])
    axis = np.cross(v, z)
    s = np.linalg.norm(axis)
    axis = axis / s
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    # angle between v and z
    theta = np.arctan2(s, c)
    R = np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)
    return R


def principal_axes(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (centroid, eigvecs, eigvals) of the coordinate covariance.

    ``eigvecs`` columns are ordered by *descending* eigenvalue. ``eigvals`` are the matching
    variances. No sign convention is imposed here (see :func:`canonical_rotation`).
    """
    coords = np.asarray(coords, dtype=float).reshape(-1, 3)
    centroid = coords.mean(axis=0)
    X = coords - centroid
    cov = (X.T @ X) / max(len(X) - 1, 1)
    vals, vecs = np.linalg.eigh(cov)  # ascending
    order = np.argsort(vals)[::-1]
    return centroid, vecs[:, order], vals[order]


def ordered_intrinsic_rotation(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Return a rigid frame determined by the ordered point cloud itself.

    The frame is a *computational gauge*, not a structural-axis estimate.  Its purpose is to
    make a subsequent coordinate-grid search independent of the input file's rotation.  The
    first basis vector points from the centroid to the first farthest point; the second points
    toward the first point with the largest component perpendicular to the first.  Both
    selections use only rotation-invariant distances and stable input order.  The third vector
    is their cross product, so the result is a proper rotation.

    A point cloud that is coincident or exactly collinear cannot define a full 3-D frame.  This
    function reports that condition with ``ValueError`` rather than silently choosing a lab-
    frame direction.  Near-collinearity is returned in ``info`` for callers that need to record
    the conditioning of the gauge.
    """
    points = np.asarray(coords, dtype=float)
    if points.ndim != 2 or points.shape[1:] != (3,) or len(points) < 2:
        raise ValueError("at least two 3-D coordinates are required to define an intrinsic frame")
    if not np.isfinite(points).all():
        raise ValueError("coordinates must be finite")

    centroid = points.mean(axis=0)
    centred = points - centroid
    squared_radius = np.einsum("ij,ij->i", centred, centred)
    radius2 = float(squared_radius.max())
    scale2 = max(radius2, 1.0)
    tie_tol = 128.0 * np.finfo(float).eps * scale2
    if radius2 <= tie_tol:
        raise ValueError("coincident coordinates cannot define an intrinsic frame")

    # Use the first point in a numerical tie. Input order is preserved by rigid transforms,
    # whereas a lab-coordinate component fallback is not rotation-equivariant.
    primary_index = int(np.flatnonzero(squared_radius >= radius2 - tie_tol)[0])
    e0 = centred[primary_index] / np.sqrt(squared_radius[primary_index])

    perpendicular = centred - (centred @ e0)[:, None] * e0[None, :]
    squared_perp = np.einsum("ij,ij->i", perpendicular, perpendicular)
    perp2 = float(squared_perp.max())
    if perp2 <= tie_tol:
        raise ValueError("collinear coordinates cannot define rotation about their line")
    secondary_index = int(np.flatnonzero(squared_perp >= perp2 - tie_tol)[0])
    e1 = perpendicular[secondary_index] / np.sqrt(squared_perp[secondary_index])
    e2 = np.cross(e0, e1)
    e2 /= np.linalg.norm(e2)
    # Recompute e1 to remove the last round-off component along e0.
    e1 = np.cross(e2, e0)

    rotation = np.vstack([e0, e1, e2])
    info = {
        "primary_index": primary_index,
        "secondary_index": secondary_index,
        "perpendicular_ratio": float(np.sqrt(perp2 / radius2)),
        "near_collinear": bool(perp2 / radius2 < 1e-12),
    }
    return rotation, centroid, info


def canonical_axis_sign(axis: np.ndarray, centred: np.ndarray) -> np.ndarray:
    """Choose an intrinsic sign for an axis without consulting lab coordinates.

    Third-moment skew supplies the usual sign.  If the distribution is symmetric along the
    axis, the first point with a non-zero projection supplies the sign.  That fallback depends
    on point order and intrinsic projections, so it transforms with the structure.  A largest-
    Cartesian-component fallback would instead change when the input file is rotated.
    """
    axis = np.asarray(axis, dtype=float).reshape(3)
    norm = float(np.linalg.norm(axis))
    if norm <= np.finfo(float).eps or not np.isfinite(norm):
        raise ValueError("axis must be a finite non-zero vector")
    axis = axis / norm
    points = np.asarray(centred, dtype=float).reshape(-1, 3)
    if not np.isfinite(points).all():
        raise ValueError("centred coordinates must be finite")
    projection = points @ axis
    m3 = float(np.mean(projection ** 3)) if len(projection) else 0.0
    scale = max(float(np.max(np.abs(projection))) if len(projection) else 0.0, 1.0)
    moment_tol = max(1e-9, 128.0 * np.finfo(float).eps * scale ** 3)
    if abs(m3) > moment_tol:
        return -axis if m3 < 0 else axis
    nonzero = np.flatnonzero(np.abs(projection) > 128.0 * np.finfo(float).eps * scale)
    if len(nonzero):
        return -axis if projection[int(nonzero[0])] < 0 else axis
    # The point cloud contains no information about this axis's sign. Preserve the caller's
    # representative; downstream code must continue to treat the plane as undirected.
    return axis


def canonical_rotation(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Deterministic PCA rotation into the canonical frame.

    Returns ``(R, centroid, info)`` such that ``(coords - centroid) @ R.T`` is the canonical
    frame. ``R`` is a proper rotation (det = +1). ``info`` reports the eigenvalues and whether
    the frame is degenerate (near-equal eigenvalues make the axis choice ill-conditioned).
    """
    centroid, vecs, vals = principal_axes(coords)
    X = np.asarray(coords, dtype=float).reshape(-1, 3) - centroid

    # Sign each of the first two axes using only intrinsic projections. Near-degenerate PCA
    # subspaces remain explicitly reported below; this does not pretend to resolve them.
    axes = [vecs[:, i].copy() for i in range(3)]
    for i in (0, 1):
        axes[i] = canonical_axis_sign(axes[i], X)
    # third axis = cross product -> guaranteed right-handed
    axes[2] = np.cross(axes[0], axes[1])

    R = np.vstack(axes)  # rows are the new basis vectors; y = R @ (x - c)
    # numerical guard: ensure proper rotation
    if np.linalg.det(R) < 0:
        R[2] = -R[2]

    # degeneracy report
    vsort = np.sort(vals)[::-1]
    denom = vsort[0] if vsort[0] > 0 else 1.0
    gap01 = (vsort[0] - vsort[1]) / denom
    gap12 = (vsort[1] - vsort[2]) / denom
    info = {
        "eigenvalues": vsort.tolist(),
        "gap01": float(gap01),
        "gap12": float(gap12),
        "degenerate": bool(gap01 < 0.02 or gap12 < 0.02),
    }
    return R, centroid, info


def canonicalize(structure: Structure) -> Tuple[Structure, dict]:
    """Return a copy of ``structure`` rotated into the canonical PCA frame, plus info.

    Deterministic and (up to eigenvalue degeneracy) frame-independent: the canonical frame of
    ``structure`` and of any rigidly-rotated copy of it are identical.
    """
    R, centroid, info = canonical_rotation(structure.ca)
    out = structure.transformed(R, t=-(R @ centroid))
    return out, info


# --------------------------------------------------------------------------------------
# IO
# --------------------------------------------------------------------------------------


def load_structure(path: str, chain: Optional[str] = None, model: int = 0) -> Structure:
    """Parse a PDB or mmCIF file into a :class:`Structure`.

    Parameters
    ----------
    path : str
        ``.pdb``/``.ent`` or ``.cif``/``.mmcif`` file.
    chain : str, optional
        Restrict to a single chain id. Default: all chains of the chosen model.
    model : int
        Model index (0-based) for multi-model files.
    """
    from Bio.PDB import MMCIFParser, PDBParser

    lower = path.lower()
    if lower.endswith((".cif", ".mmcif")):
        parser = MMCIFParser(QUIET=True)
        biostruct = parser.get_structure("s", path)
    else:
        parser = PDBParser(QUIET=True)
        try:
            biostruct = parser.get_structure("s", path)
        except (UnboundLocalError, KeyError, ValueError):
            # Some PDBs (e.g. OPM-oriented files) carry malformed REMARK/BIOMT header records
            # that crash BioPython's header parser. The coordinates are fine, so re-parse from
            # ATOM/HETATM lines only, bypassing the header.
            with open(path) as fh:
                body = "".join(
                    ln for ln in fh if ln.startswith(("ATOM", "HETATM", "TER", "MODEL", "ENDMDL", "END"))
                )
            biostruct = parser.get_structure("s", io.StringIO(body))
    models = list(biostruct.get_models())
    bio_model = models[model]
    return _structure_from_biomodel(bio_model, chain=chain, source=path)


def structure_from_string(text: str, fmt: str = "pdb", chain: Optional[str] = None) -> Structure:
    """Parse a PDB/mmCIF string (handy for tests and in-memory structures)."""
    from Bio.PDB import MMCIFParser, PDBParser

    if fmt == "cif":
        parser = MMCIFParser(QUIET=True)
        handle = io.StringIO(text)
    else:
        parser = PDBParser(QUIET=True)
        handle = io.StringIO(text)
    biostruct = parser.get_structure("s", handle)
    bio_model = list(biostruct.get_models())[0]
    return _structure_from_biomodel(bio_model, chain=chain, source="<string>")


def _structure_from_biomodel(bio_model, chain: Optional[str], source: str) -> Structure:
    ca_list: List[np.ndarray] = []
    resids: List[int] = []
    resnames: List[str] = []
    chains: List[str] = []
    sc_vecs: List[np.ndarray] = []
    plddts: List[float] = []
    icodes: List[str] = []
    atoms: List[dict] = []

    serial = 0
    for ch in bio_model:
        cid = ch.id
        if chain is not None and cid != chain:
            continue
        for res in ch:
            hetflag, resseq, icode = res.id
            resname = res.resname.strip()
            # record all atoms for serialization
            res_atoms = []
            for atom in res:
                serial += 1
                rec = {
                    "serial": serial,
                    "name": atom.get_name(),
                    "resname": resname,
                    "chain": cid,
                    "resseq": int(resseq),
                    "icode": icode.strip(),
                    "x": float(atom.coord[0]),
                    "y": float(atom.coord[1]),
                    "z": float(atom.coord[2]),
                    "element": atom.element.strip() if atom.element else atom.get_name()[0],
                    "bfactor": float(atom.get_bfactor()),
                    "hetero": hetflag.strip() != "",
                }
                atoms.append(rec)
                res_atoms.append(rec)
            # per-residue arrays only for standard amino acids with a CA
            if resname not in _AA3:
                continue
            ca_atom = None
            for atom in res:
                if atom.get_name() == "CA":
                    ca_atom = atom
                    break
            if ca_atom is None:
                continue
            ca = np.array(ca_atom.coord, dtype=float)
            ca_list.append(ca)
            resids.append(int(resseq))
            icodes.append(icode.strip())
            resnames.append("MET" if resname == "MSE" else resname)
            chains.append(cid)
            plddts.append(float(ca_atom.get_bfactor()))
            # side-chain centroid vector (heavy sidechain atoms; CB fallback; 0 for GLY)
            sc_atoms = [
                a.coord for a in res
                if a.get_name() not in ("N", "CA", "C", "O")
                and a.element.strip() not in ("H", "D")
            ]
            if sc_atoms:
                centroid = np.mean(np.asarray(sc_atoms, dtype=float), axis=0)
                v = centroid - ca
                nrm = np.linalg.norm(v)
                sc_vecs.append(v / nrm if nrm > 1e-6 else np.zeros(3))
            else:
                sc_vecs.append(np.zeros(3))

    if not ca_list:
        raise ValueError(f"no standard amino-acid residues with CA found in {source}")

    return Structure(
        ca=np.asarray(ca_list),
        resids=np.asarray(resids),
        resnames=np.asarray(resnames, dtype=object),
        chains=np.asarray(chains, dtype=object),
        sc_vec=np.asarray(sc_vecs),
        plddt=np.asarray(plddts),
        atoms=atoms,
        source=source,
        icodes=np.asarray(icodes, dtype=object),
    )


def to_pdb_string(structure: Structure) -> str:
    """Serialize the structure's atom records to a PDB string (current coordinates)."""
    lines: List[str] = []
    for a in structure.atoms:
        record = "HETATM" if a.get("hetero") else "ATOM"
        name = a["name"]
        # PDB atom-name justification: 4-char field, element-dependent
        if len(name) >= 4:
            nm = name[:4]
        elif len(a.get("element", "")) == 1 and len(name) <= 3:
            nm = " " + name.ljust(3)
        else:
            nm = name.ljust(4)
        line = (
            f"{record:<6}{a['serial'] % 100000:>5} {nm:<4} "
            f"{a['resname']:>3} {str(a['chain'])[:1]:>1}"
            f"{a['resseq'] % 10000:>4}{(a.get('icode') or ''):<1}   "
            f"{a['x']:>8.3f}{a['y']:>8.3f}{a['z']:>8.3f}"
            f"{1.00:>6.2f}{a.get('bfactor', 0.0):>6.2f}          "
            f"{a.get('element', ''):>2}"
        )
        lines.append(line)
    lines.append("END")
    return "\n".join(lines) + "\n"


def write_pdb(structure: Structure, path: str) -> None:
    with open(path, "w") as fh:
        fh.write(to_pdb_string(structure))
