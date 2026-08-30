"""memorient.orientor — the unified, context-aware entry point.

:func:`orient_structure` is the one call the rest of the platform makes. It:

1. **canonicalizes** the input into a deterministic PCA frame (so the answer never depends on
   how the coordinates happened to arrive);
2. routes to the orientation method the :class:`~memorient.contexts.MembraneContext` declares —

   * ``barrel_normal``  → :func:`~memorient.barrel.fit_membrane` + barrel/surface classifier,
     extracellular side from loop architecture (:mod:`memorient.labeler`);
   * ``tm_helix_belt``  → the experimental ``tm_helix_axis_v2`` path: exact declared
     transmembrane spans define per-helix Cα axes and an unsigned consensus normal; membrane
     centre and thickness are then optimized along that normal. Declared topology is preferred
     for sidedness. Positive-inside evidence is used only with sufficient mapped flanks and a
     predeclared charge asymmetry; otherwise sidedness remains unresolved;
   * ``anchor_relative`` → principal axis for the direction, N-terminal membrane-proximal
     anchor for the sign (no bilayer to fit);
   * ``sasa_only``      → no membrane frame at all; solvent accessibility decides exposure;
3. computes SASA, projects membrane zones **only when the context has a bilayer**, and labels
   every residue;
4. gates ``host_antibody_accessible`` through an **injected localization call** (the P1 seam)
   that can veto geometry — a residue can be geometrically extracellular yet biologically
   shielded (a periplasmic assembly, an LPS-buried surface). The default is pass-through;
5. runs :func:`five_fold_validate` — re-orients under random rotations and separately checks
   unsigned normal drift, non-vacuous embedded-residue agreement, and supported non-empty
   extracellular-residue agreement. This proves *self-consistency*, which is necessary but
   not sufficient for correctness; an external benchmark checks correctness.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Tuple

import numpy as np

from .barrel import (
    KD,
    MembraneClass,
    MembraneFit,
    classify_membrane_protein,
    fit_membrane,
    fit_membrane_on_normal,
)
from .contexts import MembraneContext, OrientationMethod
from .geometry import Structure, canonical_rotation, rotation_matrix_to_z
from .labeler import LabelSet, ResidueLabel, SideCall, call_extracellular_side, label_residues
from .membrane import (
    ACC_ANTIBODY,
    ACC_BURIED,
    MembraneProjection,
    RSA_EXPOSED,
    context_metrics,
    project_membrane,
)
from .sasa import compute_sasa


# --------------------------------------------------------------------------------------
# P1 seam: injected localization gate (replaces the hard-coded KNOWN_LOCALIZATION dict)
# --------------------------------------------------------------------------------------


@dataclass
class LocalizationCall:
    """Biological localization that can veto geometry, independent of structure.

    A geometrically-extracellular residue is only *publishable* as an antibody target if the
    protein is actually surface-exposed in the intact cell. A periplasmic flagellar sheath
    subunit is geometrically a barrel-like object with an "outer" face, yet it is shielded in
    vivo. This is the seam a localization predictor (PSORTb / SignalP / DeepTMHMM, fed by a
    ``subproteo``-style upstream) plugs into. The default is pass-through: geometry decides.
    """

    localization: str = "unknown"       # e.g. outer_membrane / periplasmic / secreted / cytoplasmic
    surface_exposed: bool = True        # does biology place any of it on the cell surface?
    source: str = "default_passthrough"
    confidence: float = 0.0

    @property
    def vetoes_surface(self) -> bool:
        return not self.surface_exposed


DEFAULT_LOCALIZATION = LocalizationCall()


# --------------------------------------------------------------------------------------
# Result container
# --------------------------------------------------------------------------------------


@dataclass
class OrientationResult:
    context: str
    method: str
    label: str                          # barrel / surface / tm_helix / anchored / soluble
    confidence: float
    structure: Structure                # oriented: membrane centred at origin, normal = +Z, EC = +Z
    rsa: np.ndarray
    labels: LabelSet
    host_antibody_accessible: bool
    localization: LocalizationCall
    fit: Optional[MembraneFit] = None
    classification: Optional[MembraneClass] = None
    side: Optional[SideCall] = None
    projection: Optional[MembraneProjection] = None
    metrics: Dict[str, float] = field(default_factory=dict)
    validation: Dict[str, object] = field(default_factory=dict)
    scientific_state: str = "placement_evaluated"
    scope_id: str = "unresolved"
    scientific_readiness: str = "prototype"
    input_normal: Optional[np.ndarray] = None
    topology_evidence: Dict[str, object] = field(default_factory=dict)

    # -- views the CLI / viz consume ----------------------------------------------------

    def residue_table(self) -> List[dict]:
        return self.labels.to_rows()

    def extracellular_resids(self) -> List[int]:
        return self.labels.extracellular_resids()

    def summary(self) -> Dict[str, object]:
        s = {
            "context": self.context,
            "method": self.method,
            "label": self.label,
            "confidence": round(self.confidence, 3),
            "n_residues": len(self.structure),
            "n_extracellular": len(self.extracellular_resids()),
            "n_surface_set": len(self.labels.surface_set),
            "host_antibody_accessible": self.host_antibody_accessible,
            "localization": self.localization.localization,
            "localization_source": self.localization.source,
            "scientific_state": self.scientific_state,
            "scope_id": self.scope_id,
            "scientific_readiness": self.scientific_readiness,
        }
        if self.fit is not None:
            s["half_thickness"] = round(self.fit.half_thickness, 2)
            s["delta_kd"] = round(self.fit.delta_kd, 2)
            s["n_embedded"] = self.fit.n_embedded
        if self.side is not None:
            s["ec_sign_confidence"] = round(self.side.confidence, 3)
        if self.validation:
            s["rotation_invariant"] = self.validation.get("passed")
            s["mean_jaccard"] = self.validation.get("mean_jaccard")
        for k, v in self.metrics.items():
            s[f"metric.{k}"] = round(float(v), 3) if v == v else None  # NaN -> None
        return s

    def to_dict(self) -> Dict[str, object]:
        surface_residue_keys = [
            {
                "chain_id": label.chain,
                "auth_seq_id": label.resid,
                "insertion_code": label.insertion_code,
            }
            for label in self.labels.labels
            if label.accessibility == ACC_ANTIBODY
        ]
        extracellular_residue_keys = [
            {
                "chain_id": label.chain,
                "auth_seq_id": label.resid,
                "insertion_code": label.insertion_code,
            }
            for label in self.labels.labels if label.extracellular
        ]
        return {
            "summary": self.summary(),
            "residues": self.residue_table(),
            "surface_set": self.labels.surface_set,
            "surface_residue_keys": surface_residue_keys,
            "extracellular_resids": self.extracellular_resids(),
            "extracellular_residue_keys": extracellular_residue_keys,
            "votes": (self.side.votes if self.side else {}),
            "metrics": {k: (float(v) if v == v else None) for k, v in self.metrics.items()},
            "validation": self.validation,
            "scientific_scope": {
                "scope_id": self.scope_id,
                "state": self.scientific_state,
                "readiness": self.scientific_readiness,
                "supported_subject_class": (
                    "beta-barrel membrane proteins" if self.scope_id == "beta_barrel"
                    else "alpha-helical membrane proteins" if self.scope_id == "alpha_helical"
                    else self.scope_id
                ),
                "known_limitations": ([
                    "Experimental method; not part of the Mark 1 qualified scope.",
                    "Membrane orientation is not native intact-cell exposure.",
                ] if self.scope_id == "alpha_helical" else [
                    "Membrane orientation is not native intact-cell exposure.",
                ]),
            },
            "topology_evidence": self.topology_evidence,
        }

    def write_pdb(self, path: str) -> None:
        from .geometry import write_pdb
        write_pdb(self.structure, path)


# --------------------------------------------------------------------------------------
# Exact topology mapping and sidedness for α-helical TM proteins
# --------------------------------------------------------------------------------------


def _residue_index(structure) -> Dict[Tuple[str, int, str], int]:
    index: Dict[Tuple[str, int, str], int] = {}
    for position, (chain, resid, icode) in enumerate(
        zip(structure.chains, structure.resids, structure.icodes)
    ):
        key = (str(chain), int(resid), str(icode))
        if key in index:
            raise ValueError(
                f"ambiguous coordinate residue mapping for {key[0]}:{key[1]}{key[2]}"
            )
        index[key] = position
    return index


def _topology_span_indices(structure, topology: Mapping[str, Any]) -> List[np.ndarray]:
    """Resolve exact declared TM spans to coordinate indices.

    Explicit residue lists are preferred because author numbering can contain gaps.  Range
    spans are accepted only when every integer position is present exactly once.
    """
    index = _residue_index(structure)
    resolved: List[np.ndarray] = []
    spans = topology.get("spans", [])
    if not isinstance(spans, list) or not spans:
        raise ValueError("topology evidence has no transmembrane spans")
    used: set[Tuple[str, int, str]] = set()
    for span_number, span in enumerate(spans, 1):
        if not isinstance(span, Mapping):
            raise ValueError(f"topology span {span_number} is not an object")
        keys: List[Tuple[str, int, str]] = []
        if isinstance(span.get("residues"), list) and span["residues"]:
            for residue in span["residues"]:
                keys.append((
                    str(residue.get("chain_id", "")),
                    int(residue["auth_seq_id"]),
                    str(residue.get("insertion_code", "")),
                ))
        else:
            chain = str(span.get("chain_id", ""))
            start = int(span.get("start_auth_seq_id", 0))
            end = int(span.get("end_auth_seq_id", -1))
            if not chain or start > end:
                raise ValueError(f"topology span {span_number} has an invalid chain or range")
            keys = [(chain, resid, "") for resid in range(start, end + 1)]
        if len(keys) < 6:
            raise ValueError(f"topology span {span_number} maps fewer than six residues")
        missing = [
            f"{chain}:{resid}{icode}" for chain, resid, icode in keys
            if (chain, resid, icode) not in index
        ]
        if missing:
            raise ValueError(f"topology span {span_number} has missing coordinate residues: {', '.join(missing[:8])}")
        if used.intersection(keys):
            raise ValueError(f"topology span {span_number} overlaps another declared span")
        used.update(keys)
        resolved.append(np.asarray([index[key] for key in keys], dtype=int))
    return resolved


def _fit_tm_helix_axis(structure, ctx: MembraneContext, topology: Mapping[str, Any]) -> MembraneFit:
    """Estimate an unsigned membrane normal from mapped transmembrane helix axes."""
    spans = _topology_span_indices(structure, topology)
    scatter = np.zeros((3, 3), dtype=float)
    for indices in spans:
        coords = structure.ca[indices]
        centered = coords - coords.mean(axis=0)
        _values, vectors = np.linalg.eigh(centered.T @ centered)
        axis = vectors[:, -1]
        axis /= np.linalg.norm(axis) + 1e-12
        scatter += len(indices) * np.outer(axis, axis)
    values, vectors = np.linalg.eigh(scatter)
    if len(values) < 3 or values[-1] <= 0:
        raise ValueError("declared transmembrane spans do not define a stable helix axis")
    normal = vectors[:, -1]
    return fit_membrane_on_normal(structure, ctx, normal)


def _tm_side_call(
    structure,
    fit: MembraneFit,
    ctx: MembraneContext,
    topology: Optional[Mapping[str, Any]] = None,
) -> tuple[SideCall, str]:
    """Extracellular-side call for a single-/multi-pass α-helical protein.

    A TM helix has no long-loop asymmetry to exploit. Declared exact topology is preferred;
    positive-inside evidence is used only when both mapped flanks are sufficiently populated
    and their predeclared charge-asymmetry threshold is met. An unsupported tie stays
    ``sides_unresolved``.
    """
    proj_along = (structure.ca - fit.centroid) @ fit.normal
    z = proj_along - fit.center
    pos = np.array([1.0 if str(r) in ("LYS", "ARG") else 0.0 for r in structure.resnames])
    plus = z > fit.half_thickness
    minus = z < -fit.half_thickness
    f_plus = pos[plus].mean() if plus.sum() else 0.0
    f_minus = pos[minus].mean() if minus.sum() else 0.0

    votes: Dict[str, int] = {}
    scores = {
        "poscharge_plus": float(f_plus), "poscharge_minus": float(f_minus),
        "flank_residues_plus": float(plus.sum()), "flank_residues_minus": float(minus.sum()),
    }
    sidedness = topology.get("sidedness", {}) if isinstance(topology, Mapping) else {}
    marker = sidedness.get("extracellular_residue") if isinstance(sidedness, Mapping) else None
    if isinstance(marker, Mapping):
        key = (
            str(marker.get("chain_id", "")),
            int(marker.get("auth_seq_id", 0)),
            str(marker.get("insertion_code", "")),
        )
        index = _residue_index(structure)
        if key not in index:
            raise ValueError(
                f"declared extracellular topology residue is unmapped: {key[0]}:{key[1]}{key[2]}"
            )
        marker_depth = float(z[index[key]])
        scores["declared_extracellular_depth"] = marker_depth
        votes["declared_topology"] = 1 if marker_depth > fit.half_thickness else -1 if marker_depth < -fit.half_thickness else 0
        if votes["declared_topology"]:
            sign = votes["declared_topology"]
            return SideCall(sign, votes, scores, 1.0, 1.0), "placement_evaluated"

    charge_gap = abs(f_plus - f_minus)
    enough_flanks = plus.sum() >= 4 and minus.sum() >= 4
    if enough_flanks and charge_gap >= 0.10:
        votes["positive_inside"] = 1 if f_plus < f_minus else -1
    else:
        votes["positive_inside"] = 0

    # secondary: which face carries more polar/charged (ecto) mass, and terminus location
    term_z = np.array([z[0], z[-1]])
    off = np.abs(term_z) > fit.half_thickness
    if off.any():
        term_side = float(np.sign(term_z[off].mean()))
        scores["terminus_side"] = term_side
        # no universal rule for helices, so terminus is only a light tiebreaker
        votes["terminus"] = 0
    else:
        scores["terminus_side"] = 0.0
        votes["terminus"] = 0

    if votes["positive_inside"] == 0:
        # +1 is only an unsigned coordinate-frame convention here.  The returned state and
        # zero confidence prevent it from becoming an extracellular-side assertion.
        return SideCall(1, votes, scores, 0.0, 0.0), "sides_unresolved"
    ec_sign = votes["positive_inside"]
    conf = min(0.4 + charge_gap * 2.0, 0.9)
    return SideCall(ec_sign, votes, scores, 1.0, float(conf)), "placement_evaluated"


# --------------------------------------------------------------------------------------
# Labels for contexts without a bilayer
# --------------------------------------------------------------------------------------


def _label_no_membrane(structure, rsa: np.ndarray, outward: Optional[np.ndarray] = None) -> LabelSet:
    """Label residues when there is no bilayer: exposure alone defines the surface set.

    ``outward`` (optional, signed distance along an anchor axis) tags residues as
    'extracellular'/distal when positive; without it, every exposed residue is surface.
    """
    labels: List[ResidueLabel] = []
    surface: List[int] = []
    for i in range(len(structure)):
        exposed = rsa[i] >= RSA_EXPOSED
        is_ec = bool(exposed and (outward is None or outward[i] > 0))
        acc = ACC_ANTIBODY if (exposed and is_ec) else ACC_BURIED
        lbl = ResidueLabel(
            resid=int(structure.resids[i]), resname=str(structure.resnames[i]),
            chain=str(structure.chains[i]), zone="", facing="",
            accessibility=acc, extracellular=is_ec, rsa=float(rsa[i]),
            ec_depth=float(outward[i]) if outward is not None else 0.0,
            insertion_code=str(structure.icodes[i]),
        )
        labels.append(lbl)
        if acc == ACC_ANTIBODY:
            surface.append(lbl.resid)
    return LabelSet(labels=labels, surface_set=surface)


# --------------------------------------------------------------------------------------
# Re-frame the structure into the oriented membrane frame (normal +Z, EC on top)
# --------------------------------------------------------------------------------------


def _reframe_to_membrane(structure, fit: MembraneFit, ec_sign: int) -> Structure:
    """Rotate so the (signed) membrane normal maps to +Z, translate slab centre to origin."""
    n = ec_sign * fit.normal
    R = rotation_matrix_to_z(n)
    # membrane centre in lab coords
    center_pt = fit.centroid + fit.center * fit.normal
    return structure.transformed(R, t=-R @ center_pt)


# --------------------------------------------------------------------------------------
# Rotation-invariance validation
# --------------------------------------------------------------------------------------


def _jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def _random_rotation(rng) -> np.ndarray:
    A = rng.normal(size=(3, 3))
    Q, R = np.linalg.qr(A)
    Q *= np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


def orient_structure(
    structure, context: MembraneContext,
    localization: Optional[LocalizationCall] = None,
    validate: bool = True, n_points: int = 240, n_validate_seeds: int = 8,
    topology_evidence: Optional[Mapping[str, Any]] = None,
) -> OrientationResult:
    """Orient and label a structure in its membrane context. The one public entry point."""
    loc = localization if localization is not None else DEFAULT_LOCALIZATION

    # 1. canonical frame (deterministic; downstream is frame-independent regardless)
    canonical_R, canonical_centroid, _info = canonical_rotation(structure.ca)
    canon = structure.transformed(canonical_R, t=-(canonical_R @ canonical_centroid))

    # 2. SASA (always)
    sasa_res = compute_sasa(canon, n_points=n_points)
    rsa = sasa_res["rsa"]

    method = context.orientation_method
    fit = None
    classification = None
    side = None
    projection = None
    metrics: Dict[str, float] = {}
    label = "unknown"
    scientific_state = "placement_evaluated"
    scope_id = "unresolved"
    scientific_readiness = "prototype"
    input_normal = None
    topology_summary: Dict[str, object] = {}

    if method in (OrientationMethod.BARREL_NORMAL, OrientationMethod.TM_HELIX_BELT):
        if method == OrientationMethod.BARREL_NORMAL:
            fit = fit_membrane(canon, context)
            classification = classify_membrane_protein(canon, context, fit)
            label = classification.label
            side = call_extracellular_side(canon, fit, context)
            scope_id = "beta_barrel"
            scientific_readiness = "conditionally_qualified"
        else:  # TM_HELIX_BELT
            scope_id = "alpha_helical"
            scientific_readiness = "prototype"
            label = "tm_helix_experimental"
            if topology_evidence and topology_evidence.get("spans"):
                fit = _fit_tm_helix_axis(canon, context, topology_evidence)
                method = "tm_helix_axis_v2"
                topology_summary = {
                    "state": "checksum_bound_spans_supplied",
                    "span_count": len(topology_evidence.get("spans", [])),
                    "source": topology_evidence.get("source", {}),
                }
            else:
                # Compatibility-only path: it may produce an inspectable plane, but it is
                # never allowed to carry qualified alpha-helical or sidedness claims.
                fit = fit_membrane(canon, context)
                method = "tm_helix_belt_legacy_experimental"
                scientific_state = "insufficient_topology"
                topology_summary = {"state": "missing", "span_count": 0}
            side, side_state = _tm_side_call(canon, fit, context, topology_evidence)
            if scientific_state != "insufficient_topology":
                scientific_state = side_state

        input_normal = canonical_R.T @ fit.normal
        oriented = _reframe_to_membrane(canon, fit, side.ec_sign)
        # recompute SASA + fit-derived quantities in the oriented frame for the labels
        rsa_o = compute_sasa(oriented, n_points=n_points)["rsa"]
        if scope_id == "alpha_helical" and topology_evidence and topology_evidence.get("spans"):
            fit_o = _fit_tm_helix_axis(oriented, context, topology_evidence)
        else:
            fit_o = fit_membrane(oriented, context)
        if scope_id == "beta_barrel":
            side_o = call_extracellular_side(oriented, fit_o, context)
        else:
            side_o, side_state = _tm_side_call(oriented, fit_o, context, topology_evidence)
            if scientific_state != "insufficient_topology":
                scientific_state = side_state
        projection = project_membrane(oriented, fit_o, context, ec_sign=side_o.ec_sign, rsa=rsa_o)
        metrics = context_metrics(oriented, fit_o, context, projection)
        labels = label_residues(oriented, projection, rsa_o, context, fit_o)
        if scope_id == "alpha_helical" and scientific_state != "placement_evaluated":
            # Preserve membrane depth/core evidence while suppressing an unsupported side.
            labels.surface_set = []
            for residue_label in labels.labels:
                residue_label.extracellular = False
                if residue_label.zone != "hydrophobic_core":
                    residue_label.zone = "side_unresolved"
                residue_label.accessibility = "sidedness_unresolved"
        rsa = rsa_o
        fit = fit_o
        side = side_o
        struct_out = oriented

    elif method == OrientationMethod.ANCHOR_RELATIVE:
        scope_id = "anchored_surface"
        # principal axis = outward direction; N-terminal anchor breaks the sign
        from .geometry import principal_axes
        centroid = canon.ca.mean(axis=0)
        _, vecs, _ = principal_axes(canon.ca)
        axis = vecs[:, 0]
        proj = (canon.ca - centroid) @ axis
        # N-terminal residues are membrane-proximal (anchored); outward points away from them
        n_term = proj[: max(3, len(proj) // 20)].mean()
        if n_term > 0:
            axis = -axis
            proj = -proj
        outward = proj - proj.min()  # distal (surface) residues have large positive outward
        # centre the outward score so exposed distal half is "extracellular"
        outward = proj - np.median(proj)
        R = rotation_matrix_to_z(axis)
        struct_out = canon.transformed(R, t=-R @ centroid)
        rsa = compute_sasa(struct_out, n_points=n_points)["rsa"]
        out2 = (struct_out.ca - struct_out.ca.mean(axis=0)) @ np.array([0, 0, 1.0])
        labels = _label_no_membrane(struct_out, rsa, outward=out2)
        label = "anchored"

    else:  # SASA_ONLY
        scope_id = "soluble_surface"
        struct_out = canon
        labels = _label_no_membrane(struct_out, rsa, outward=None)
        label = "soluble"

    # 4. localization gate (P1 seam): biology can veto geometry
    geometry_says_surface = len(labels.surface_set) > 0
    host_antibody_accessible = bool(geometry_says_surface and not loc.vetoes_surface)

    # 5. rotation-invariance validation
    validation: Dict[str, object] = {}
    if validate:
        validation = _run_five_fold(
            structure, context, loc, n_validate_seeds, topology_evidence=topology_evidence,
        )

    conf = 0.5
    if classification is not None:
        conf = classification.confidence
    elif side is not None:
        conf = side.confidence
    elif method == OrientationMethod.SASA_ONLY:
        conf = 1.0

    return OrientationResult(
        context=context.name, method=method, label=label, confidence=float(conf),
        structure=struct_out, rsa=rsa, labels=labels,
        host_antibody_accessible=host_antibody_accessible, localization=loc,
        fit=fit, classification=classification, side=side, projection=projection,
        metrics=metrics, validation=validation, scientific_state=scientific_state,
        scope_id=scope_id, scientific_readiness=scientific_readiness,
        input_normal=input_normal, topology_evidence=topology_summary,
    )


def _evidence_sets(structure, context, loc, n_points, topology_evidence=None):
    result = orient_structure(
        structure, context, localization=loc, validate=False, n_points=n_points,
        topology_evidence=topology_evidence,
    )
    extracellular = {
        (label.chain, label.resid, label.insertion_code)
        for label in result.labels.labels if label.extracellular
    }
    embedded = {
        (label.chain, label.resid, label.insertion_code) for label in result.labels.labels
        if label.zone == "hydrophobic_core"
    }
    if result.fit is None:
        embedded = {
            (label.chain, label.resid, label.insertion_code) for label in result.labels.labels
            if label.accessibility == ACC_ANTIBODY
        }
    return result, embedded, extracellular


def _run_five_fold(
    structure, context, loc, seeds, n_points: int = 160,
    threshold: float = 0.95, normal_threshold_deg: float = 1.0,
    topology_evidence: Optional[Mapping[str, Any]] = None,
) -> Dict[str, object]:
    """Repeat orientation under rotations without allowing empty sets to pass."""
    reference, reference_embedded, reference_extracellular = _evidence_sets(
        structure, context, loc, n_points, topology_evidence,
    )
    rng = np.random.default_rng(0)
    embedded_jaccards: List[Optional[float]] = []
    extracellular_jaccards: List[Optional[float]] = []
    normal_errors: List[Optional[float]] = []
    for _ in range(seeds):
        R = _random_rotation(rng)
        rotated = structure.transformed(R)
        result, embedded, extracellular = _evidence_sets(
            rotated, context, loc, n_points, topology_evidence,
        )
        embedded_jaccards.append(
            None if not reference_embedded and not embedded else _jaccard(reference_embedded, embedded)
        )
        extracellular_jaccards.append(
            None if not reference_extracellular and not extracellular
            else _jaccard(reference_extracellular, extracellular)
        )
        if reference.input_normal is not None and result.input_normal is not None:
            expected = R @ reference.input_normal
            cosine = min(abs(float(np.dot(expected, result.input_normal))), 1.0)
            normal_errors.append(float(np.degrees(np.arccos(cosine))))
        else:
            normal_errors.append(None)
    applicable_embedded = [value for value in embedded_jaccards if value is not None]
    applicable_extracellular = [value for value in extracellular_jaccards if value is not None]
    applicable_normal = [value for value in normal_errors if value is not None]
    placement_passed = bool(
        applicable_embedded
        and all(value >= threshold for value in applicable_embedded)
        and (not applicable_normal or all(value <= normal_threshold_deg for value in applicable_normal))
    )
    sidedness_passed = bool(
        reference.scientific_state == "placement_evaluated"
        and applicable_extracellular
        and all(value >= threshold for value in applicable_extracellular)
    )
    mean_j = float(np.mean(applicable_extracellular)) if applicable_extracellular else None
    return {
        "passed": bool(placement_passed and (sidedness_passed if context.has_membrane_sides else True)),
        "placement_passed": placement_passed,
        "sidedness_passed": sidedness_passed if context.has_membrane_sides else None,
        "jaccards": [None if value is None else round(value, 3) for value in extracellular_jaccards],
        "mean_jaccard": None if mean_j is None else round(mean_j, 3),
        "embedded_jaccards": [None if value is None else round(value, 3) for value in embedded_jaccards],
        "normal_angle_errors_deg": [None if value is None else round(value, 6) for value in normal_errors],
        "threshold": threshold,
        "normal_threshold_deg": normal_threshold_deg,
        "n_reference_embedded": len(reference_embedded),
        "n_reference_extracellular": len(reference_extracellular),
        "extracellular_comparison_state": (
            "evaluated" if applicable_extracellular else "not_applicable_empty_or_unresolved"
        ),
    }


def five_fold_validate(structure, context: MembraneContext,
                       localization: Optional[LocalizationCall] = None,
                       seeds: int = 8, threshold: float = 0.95,
                       n_points: int = 160,
                       topology_evidence: Optional[Mapping[str, Any]] = None,
                       normal_threshold_deg: float = 1.0) -> Dict[str, object]:
    """Public rotation-invariance check for placement, normal, and supported sidedness.

    Embedded-residue agreement must be non-vacuous, unsigned normal drift is
    bounded when a membrane axis exists, and extracellular-set agreement is
    evaluated only when the reference and rotated runs both carry a supported,
    non-empty sided set. Empty/empty sidedness is ``not_applicable`` rather than
    a favorable Jaccard.

    Self-consistency, not correctness — a stable wrong answer still passes. The P4 benchmark
    checks the membrane normal against experimentally-oriented references.

    The rotation count is ``seeds`` and is now 8. It was 5, and the ``five_fold``
    names are historical: they record the original default, not a fixed protocol.
    Five rotations under-sampled the basin structure -- on Qualification v2's
    beta_barrel stratum, cases that showed zero drift over five rotations in one
    environment drifted by more than 8 degrees over the same five in another, so
    the sample was small enough that whether a second basin was visited at all
    depended on the machine. Any caller that needs the historical behaviour can
    pass ``seeds=5`` explicitly.
    """
    loc = localization if localization is not None else DEFAULT_LOCALIZATION
    return _run_five_fold(
        structure, context, loc, seeds, n_points=n_points, threshold=threshold,
        normal_threshold_deg=normal_threshold_deg, topology_evidence=topology_evidence,
    )
