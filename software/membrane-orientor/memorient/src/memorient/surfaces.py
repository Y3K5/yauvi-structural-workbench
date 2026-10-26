"""Experimental geometric membranes and rigid placement, in Angstroms.

No lipid packing, biological sidedness inference, or energy minimization.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


def _vector(value, name):
    result = np.asarray(value, dtype=float)
    if result.shape != (3,) or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain three finite numbers")
    return result


def _unit(value, name):
    result = _vector(value, name)
    length = np.linalg.norm(result)
    if length < 1e-12:
        raise ValueError(f"{name} must be nonzero")
    return result / length


@dataclass(frozen=True)
class Surface:
    """Axis-aligned midsurface; dimensions are explicit, never species defaults.

    plane: size=(width,height); sphere: size=(radius,); cylinder:
    size=(radius,length), axis z; ellipsoid: size=(a,b,c).
    Parameters u,v: plane x/y Angstroms; sphere/ellipsoid longitude and
    latitude radians; cylinder azimuth radians and axial z Angstroms.
    """
    kind: str
    size: tuple
    thickness: float
    center: tuple = (0., 0., 0.)

    def __post_init__(self):
        counts = {"plane": 2, "sphere": 1, "cylinder": 2, "ellipsoid": 3}
        if self.kind not in counts:
            raise ValueError("unsupported surface kind")
        size = np.asarray(self.size, dtype=float)
        if size.shape != (counts[self.kind],) or not np.isfinite(size).all() or np.any(size <= 0):
            raise ValueError("surface dimensions must be finite and positive")
        if not math.isfinite(self.thickness) or self.thickness <= 0:
            raise ValueError("thickness must be finite and positive")
        _vector(self.center, "center")
        # Conservative bound prevents inward offset folding at high curvature.
        if self.kind in {"sphere", "cylinder"}:
            curvature_radius = size[0]
        elif self.kind == "ellipsoid":
            curvature_radius = min(size)**2 / max(size)
        else:
            curvature_radius = float("inf")
        if self.thickness / 2 >= curvature_radius:
            raise ValueError("half thickness exceeds safe local curvature radius")

    def frame(self, u, v):
        """Return midsurface position and geometric outward normal."""
        if not np.isfinite([u, v]).all():
            raise ValueError("surface parameters must be finite")
        if self.kind == "plane":
            if abs(u) > self.size[0]/2 or abs(v) > self.size[1]/2:
                raise ValueError("point outside plane patch")
            point, normal = np.array([u, v, 0.]), np.array([0., 0., 1.])
        elif self.kind == "cylinder":
            r, length = self.size
            if abs(v) > length/2:
                raise ValueError("point outside open cylinder")
            normal = np.array([math.cos(u), math.sin(u), 0.])
            point = r * normal + [0., 0., v]
        else:
            if abs(v) > math.pi/2:
                raise ValueError("latitude outside [-pi/2, pi/2]")
            axes = np.repeat(self.size[0], 3) if self.kind == "sphere" else np.asarray(self.size)
            direction = np.array([math.cos(v)*math.cos(u), math.cos(v)*math.sin(u), math.sin(v)])
            point = axes * direction
            normal = _unit(direction / axes, "normal")
        return point + _vector(self.center, "center"), normal

    def sample(self, nu=32, nv=17):
        """Point cloud of midsurface and offset interfaces; no molecular atoms."""
        if type(nu) is not int or type(nv) is not int or min(nu, nv) < 2 or nu*nv > 100000:
            raise ValueError("sample grid requires integers >=2 and at most 100000 points")
        if self.kind == "plane":
            us = np.linspace(-self.size[0]/2, self.size[0]/2, nu)
            vs = np.linspace(-self.size[1]/2, self.size[1]/2, nv)
        else:
            us = np.linspace(-math.pi, math.pi, nu, endpoint=False)
            vs = np.linspace(-self.size[1]/2, self.size[1]/2, nv) if self.kind == "cylinder" else np.linspace(-math.pi/2, math.pi/2, nv)
        pairs = [self.frame(u, v) for u in us for v in vs]
        points, normals = np.array([x[0] for x in pairs]), np.array([x[1] for x in pairs])
        return {"mid": points.tolist(), "positive_interface": (points + normals*self.thickness/2).tolist(),
                "negative_interface": (points - normals*self.thickness/2).tolist(), "normals": normals.tolist()}


def place_receptor(coordinates, surface, *, anchor, axis, u, v, side=1, depth=0., spin=0.):
    """Rigidly align a supplied receptor axis to a local normal.

    axis sign, anchor, insertion depth and side are user assertions, not
    inferred biological evidence. Row vectors transform as xyz @ R.T + t.
    """
    xyz = np.asarray(coordinates, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError("coordinates must be a nonempty finite Nx3 array")
    if side not in (-1, 1) or not np.isfinite([depth, spin]).all():
        raise ValueError("side must be +/-1 and depth/spin finite")
    anchor = _vector(anchor, "anchor")
    source = _unit(axis, "axis")
    point, normal = surface.frame(u, v)
    target = side * normal
    cross = np.cross(source, target)
    cosine = np.clip(source @ target, -1., 1.)
    if cosine < -1 + 1e-12:
        perpendicular = _unit(np.cross(source, np.eye(3)[np.argmin(abs(source))]), "perpendicular")
        rotation = 2*np.outer(perpendicular, perpendicular) - np.eye(3)
    else:
        skew = np.array([[0., -cross[2], cross[1]], [cross[2], 0., -cross[0]], [-cross[1], cross[0], 0.]])
        rotation = np.eye(3) + skew + skew @ skew / (1 + cosine)
    skew = np.array([[0., -target[2], target[1]], [target[2], 0., -target[0]], [-target[1], target[0], 0.]])
    twist = np.eye(3) + math.sin(spin)*skew + (1-math.cos(spin))*(skew @ skew)
    rotation = twist @ rotation
    translation = point + depth*target - rotation @ anchor
    placed = xyz @ rotation.T + translation
    return {"coordinates": placed.tolist(), "rotation": rotation.tolist(), "translation": translation.tolist(),
            "surface_point": point.tolist(), "target_axis": target.tolist(),
            "scientific_state": "geometric_placement_only", "sidedness": "user_asserted",
            "limitations": ["No lipid atoms or force field", "No collision or energy assessment", "Local tangent placement does not establish curved-membrane compatibility", "No native exposure or biological sidedness validation"]}
