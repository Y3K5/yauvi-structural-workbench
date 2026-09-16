"""Synthetic structure builders for tests — ideal barrel, TM helix, soluble blob.

These produce :class:`memorient.geometry.Structure` objects with realistic-enough geometry
and residue identity that the numeric metrics fire in the right direction, without needing a
network fetch. They are the ground truth for the rotation-invariance and routing tests: we
*know* where the membrane is because we placed it.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np

from memorient.geometry import Structure

# hydrophobic / polar / charged pools for placing residues by role
HYDROPHOBIC = ["LEU", "ILE", "VAL", "PHE", "ALA", "MET"]
AROMATIC_GIRDLE = ["TRP", "TYR"]
POLAR = ["SER", "THR", "ASN", "GLN"]
POS = ["LYS", "ARG"]
NEG = ["ASP", "GLU"]


def _atoms_for_residue(serial0: int, resname: str, chain: str, resseq: int,
                       ca: np.ndarray, sc_dir: Optional[np.ndarray], bf: float) -> Tuple[list, np.ndarray]:
    """Minimal atom set: N, CA, C, O, and (unless GLY) a CB placed along ``sc_dir``.

    Returns (atom_records, sc_vec_unit). ``sc_dir`` is the intended side-chain direction
    (unit); CB is placed 1.5 A along it from CA.
    """
    atoms = []
    s = serial0

    def add(name, xyz, elem):
        nonlocal s
        s += 1
        atoms.append({
            "serial": s, "name": name, "resname": resname, "chain": chain,
            "resseq": resseq, "icode": "", "x": float(xyz[0]), "y": float(xyz[1]),
            "z": float(xyz[2]), "element": elem, "bfactor": float(bf), "hetero": False,
        })

    # place backbone atoms as a small rigid cluster around CA (geometry need only be plausible)
    add("N", ca + np.array([-0.6, 0.8, 0.0]), "N")
    add("CA", ca, "C")
    add("C", ca + np.array([0.6, -0.8, 0.0]), "C")
    add("O", ca + np.array([1.2, -1.4, 0.3]), "O")
    sc_vec = np.zeros(3)
    if resname != "GLY" and sc_dir is not None:
        d = np.asarray(sc_dir, dtype=float)
        n = np.linalg.norm(d)
        if n > 1e-6:
            d = d / n
            sc_vec = d
            add("CB", ca + 1.5 * d, "C")
            # a second side-chain atom further out so the centroid vector is well-defined
            add("CG", ca + 3.0 * d, "C")
    return atoms, sc_vec


def make_barrel(n_strands: int = 10, strand_len: int = 10, radius: float = 12.0,
                half_thick: float = 13.0, ec_loop_len: int = 8, peri_loop_len: int = 2,
                seed: int = 0) -> Structure:
    """An ideal antiparallel beta-barrel spanning the membrane along Z.

    * Strands run parallel to Z between ``-half_thick`` and ``+half_thick``.
    * Side chains on odd strand-residues point *outward* (lipid) and are hydrophobic;
      even ones point *inward* (pore) and are polar — the barrel signature.
    * Aromatics (Trp/Tyr) are planted at both membrane interfaces (the girdle).
    * Long extracellular loops at +Z, short periplasmic turns at -Z (the side signal).
    The extracellular side is +Z by construction.
    """
    rng = np.random.default_rng(seed)
    ca: List[np.ndarray] = []
    resnames: List[str] = []
    resids: List[int] = []
    chains: List[str] = []
    sc_vecs: List[np.ndarray] = []
    atoms: List[dict] = []
    serial = 0
    resid = 1

    def emit(pos, resname, sc_dir, bf=80.0):
        nonlocal serial, resid
        recs, sv = _atoms_for_residue(serial, resname, "A", resid, pos, sc_dir, bf)
        serial = recs[-1]["serial"] if recs else serial
        atoms.extend(recs)
        ca.append(np.asarray(pos, dtype=float))
        resnames.append(resname)
        resids.append(resid)
        chains.append("A")
        sc_vecs.append(sv)
        resid += 1

    zc = np.linspace(-half_thick, half_thick, strand_len)
    for k in range(n_strands):
        ang = 2 * np.pi * k / n_strands
        radial = np.array([np.cos(ang), np.sin(ang), 0.0])
        going_up = (k % 2 == 0)
        zseq = zc if going_up else zc[::-1]
        for i, z in enumerate(zseq):
            pos = np.array([radius * np.cos(ang), radius * np.sin(ang), z])
            outward = (i % 2 == 0)
            # interface aromatic girdle near |z| ~ half_thick
            at_interface = abs(z) > half_thick - 3.0
            if at_interface and rng.random() < 0.6:
                resname = AROMATIC_GIRDLE[int(rng.integers(0, 2))]
                sc_dir = radial
            elif outward:
                resname = HYDROPHOBIC[int(rng.integers(0, len(HYDROPHOBIC)))]
                sc_dir = radial  # lipid-facing
            else:
                resname = POLAR[int(rng.integers(0, len(POLAR)))]
                sc_dir = -radial  # pore-facing
            emit(pos, resname, sc_dir)
        # connecting loop after this strand
        top = going_up  # ended at +Z if going up
        loop_len = ec_loop_len if top else peri_loop_len
        z_end = half_thick if top else -half_thick
        next_ang = 2 * np.pi * (k + 1) / n_strands
        for j in range(loop_len):
            frac = (j + 1) / (loop_len + 1)
            bulge = 6.0 * np.sin(np.pi * frac)  # loop pokes away from the barrel
            # extracellular loops rise well clear of the headgroup interface (real OMP loops
            # extend 10-30 A); periplasmic turns stay short.
            rise = (16.0 * np.sin(np.pi * frac) if top else 2.0 * np.sin(np.pi * frac))
            zpos = z_end + (rise if top else -rise)
            a = ang + (next_ang - ang) * frac
            r = radius + bulge
            pos = np.array([r * np.cos(a), r * np.sin(a), zpos])
            resname = POLAR[int(rng.integers(0, len(POLAR)))]
            emit(pos, resname, np.array([np.cos(a), np.sin(a), 0.0]))

    return Structure(
        ca=np.asarray(ca), resids=np.asarray(resids),
        resnames=np.asarray(resnames, dtype=object), chains=np.asarray(chains, dtype=object),
        sc_vec=np.asarray(sc_vecs), plddt=np.full(len(ca), 80.0), atoms=atoms,
        source="<synthetic barrel>",
    )


def make_tm_helix(n_res: int = 30, half_thick: float = 15.0, rise: float = 1.5,
                  turn: float = 100.0, extra_ecto: int = 10, extra_cyto: int = 10,
                  seed: int = 1) -> Structure:
    """A single-pass alpha-helix crossing the membrane along Z, with flanking soluble tails.

    * Core (|z| < half_thick) residues are hydrophobic (the TM belt).
    * The cytoplasmic tail (-Z) is enriched in Lys/Arg (positive-inside).
    * The ectodomain tail (+Z) is polar/negative.
    Extracellular side is +Z by construction.
    """
    rng = np.random.default_rng(seed)
    ca: List[np.ndarray] = []
    resnames: List[str] = []
    resids: List[int] = []
    chains: List[str] = []
    sc_vecs: List[np.ndarray] = []
    atoms: List[dict] = []
    serial = 0
    resid = 1

    def emit(pos, resname, sc_dir, bf=75.0):
        nonlocal serial, resid
        recs, sv = _atoms_for_residue(serial, resname, "A", resid, pos, sc_dir, bf)
        serial = recs[-1]["serial"] if recs else serial
        atoms.extend(recs)
        ca.append(np.asarray(pos, dtype=float)); resnames.append(resname)
        resids.append(resid); chains.append("A"); sc_vecs.append(sv); resid += 1

    r_helix = 2.3
    # z spans from -(half_thick+extra_cyto*rise) to +(half_thick+extra_ecto*rise)
    n_core = int(2 * half_thick / rise) + 1
    total = extra_cyto + n_core + extra_ecto
    z0 = -half_thick - extra_cyto * rise
    for i in range(total):
        z = z0 + i * rise
        ang = np.deg2rad(turn) * i
        pos = np.array([r_helix * np.cos(ang), r_helix * np.sin(ang), z])
        radial = np.array([np.cos(ang), np.sin(ang), 0.0])
        if abs(z) <= half_thick:
            resname = HYDROPHOBIC[int(rng.integers(0, len(HYDROPHOBIC)))]
        elif z < -half_thick:  # cytoplasmic
            resname = POS[int(rng.integers(0, len(POS)))] if rng.random() < 0.6 else POLAR[int(rng.integers(0, len(POLAR)))]
        else:  # extracellular ectodomain
            resname = NEG[int(rng.integers(0, len(NEG)))] if rng.random() < 0.4 else POLAR[int(rng.integers(0, len(POLAR)))]
        emit(pos, resname, radial)

    return Structure(
        ca=np.asarray(ca), resids=np.asarray(resids),
        resnames=np.asarray(resnames, dtype=object), chains=np.asarray(chains, dtype=object),
        sc_vec=np.asarray(sc_vecs), plddt=np.full(len(ca), 75.0), atoms=atoms,
        source="<synthetic tm helix>",
    )


def make_soluble_blob(n_res: int = 120, radius: float = 16.0, seed: int = 2) -> Structure:
    """A compact globular blob with no membrane character (roughly spherical)."""
    rng = np.random.default_rng(seed)
    pts = rng.normal(size=(n_res, 3))
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    rad = radius * rng.random(n_res) ** (1.0 / 3.0)
    coords = pts * rad[:, None]
    ca: List[np.ndarray] = []
    resnames: List[str] = []
    resids: List[int] = []
    chains: List[str] = []
    sc_vecs: List[np.ndarray] = []
    atoms: List[dict] = []
    serial = 0
    pool = HYDROPHOBIC + POLAR + POS + NEG
    for i in range(n_res):
        pos = coords[i]
        outward = pos / (np.linalg.norm(pos) + 1e-9)
        resname = pool[int(rng.integers(0, len(pool)))]
        recs, sv = _atoms_for_residue(serial, resname, "A", i + 1, pos, outward, 70.0)
        serial = recs[-1]["serial"]
        atoms.extend(recs)
        ca.append(pos); resnames.append(resname); resids.append(i + 1)
        chains.append("A"); sc_vecs.append(sv)
    return Structure(
        ca=np.asarray(ca), resids=np.asarray(resids),
        resnames=np.asarray(resnames, dtype=object), chains=np.asarray(chains, dtype=object),
        sc_vec=np.asarray(sc_vecs), plddt=np.full(n_res, 70.0), atoms=atoms,
        source="<synthetic soluble>",
    )


def make_ellipsoid(n_res: int = 150, axes: tuple = (24.0, 15.0, 8.0), seed: int = 3) -> Structure:
    """A blob with three DISTINCT principal axes (non-degenerate) for strict frame tests.

    Because the semi-axes differ, the PCA frame is uniquely defined and canonicalization is
    exactly frame-independent — unlike the symmetric barrel/helix/sphere where an in-plane
    rotation is a true symmetry.
    """
    rng = np.random.default_rng(seed)
    pts = rng.normal(size=(n_res, 3))
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    rad = rng.random(n_res) ** (1.0 / 3.0)
    coords = pts * rad[:, None] * np.asarray(axes)[None, :]
    ca: List[np.ndarray] = []
    resnames: List[str] = []
    resids: List[int] = []
    chains: List[str] = []
    sc_vecs: List[np.ndarray] = []
    atoms: List[dict] = []
    serial = 0
    pool = HYDROPHOBIC + POLAR + POS + NEG
    for i in range(n_res):
        pos = coords[i]
        outward = pos / (np.linalg.norm(pos) + 1e-9)
        resname = pool[int(rng.integers(0, len(pool)))]
        recs, sv = _atoms_for_residue(serial, resname, "A", i + 1, pos, outward, 70.0)
        serial = recs[-1]["serial"]
        atoms.extend(recs)
        ca.append(pos); resnames.append(resname); resids.append(i + 1)
        chains.append("A"); sc_vecs.append(sv)
    return Structure(
        ca=np.asarray(ca), resids=np.asarray(resids),
        resnames=np.asarray(resnames, dtype=object), chains=np.asarray(chains, dtype=object),
        sc_vec=np.asarray(sc_vecs), plddt=np.full(n_res, 70.0), atoms=atoms,
        source="<synthetic ellipsoid>",
    )


def random_rotation(seed: int = 0) -> np.ndarray:
    """A uniformly random proper rotation matrix (QR method, det forced +1)."""
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(3, 3))
    Q, R = np.linalg.qr(A)
    Q *= np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q
