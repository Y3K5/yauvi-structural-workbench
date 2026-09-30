"""Mapped compartment evidence in an unsigned input-coordinate membrane frame."""
from __future__ import annotations
import numpy as np


def resolve_sidedness(atoms, center, normal, half_thickness, context, topology=None):
    from .contexts import get_context
    ctx = get_context(context) if isinstance(context, str) else context
    topology = topology or {}
    compartments = tuple(topology.get("compartments", ctx.compartments))
    if len(compartments) != 2 or len(set(compartments)) != 2 or not all(isinstance(c, str) and c.strip() for c in compartments):
        raise ValueError("Two distinct compartment names are required")
    if ctx.name == "organelle_membrane" and compartments == ctx.compartments:
        raise ValueError("Other organelle membranes require explicit compartment names")
    if ctx.name != "organelle_membrane" and compartments != ctx.compartments:
        raise ValueError("Compartment names conflict with the declared membrane context")
    result = {"state": "unknown", "compartments": list(compartments), "positive_compartment": None,
              "negative_compartment": None, "markers": [], "source": topology.get("source", {}),
              "interpretation_limit": "Mapped compartment assignment is not native accessibility or qualified membrane accuracy."}
    markers = list(topology.get("side_markers", []))
    legacy = topology.get("sidedness", {}).get("extracellular_residue")
    if legacy:
        if "extracellular" not in compartments:
            raise ValueError("Extracellular marker is incompatible with this membrane")
        markers.append({**legacy, "compartment": "extracellular"})
    votes, bases = [], []
    for marker in markers:
        if marker.get("compartment") not in compartments:
            raise ValueError("Marker names a different compartment")
        matches = [a for a in atoms if a["atom"] == "CA" and marker.get("chain_id") in {a["chain_id"], a.get("chain_copy_id")}
                   and a["auth_seq_id"] == marker.get("auth_seq_id") and a.get("insertion_code", "") == marker.get("insertion_code", "")]
        if len(matches) != 1:
            raise ValueError("Sidedness marker must map to one exact coordinate residue")
        depth = float((np.array(matches[0]["xyz"]) - center) @ normal)
        sign = 1 if depth > half_thickness else -1 if depth < -half_thickness else 0
        assignment = marker["compartment"] if sign == 1 else compartments[1 - compartments.index(marker["compartment"])] if sign == -1 else None
        result["markers"].append({**marker, "depth_A": depth, "positive_compartment": assignment})
        if assignment:
            votes.append(assignment)
        bases.append(marker.get("basis", topology.get("source", {}).get("basis", "curated")))
    if markers and (not result["source"].get("id") or not result["source"].get("citation")):
        raise ValueError("Sidedness markers require a source id and citation")
    if any(b not in {"curated", "experimental", "predicted", "derived"} for b in bases):
        raise ValueError("Unsupported sidedness evidence basis")
    if len(set(votes)) > 1:
        result["state"] = "conflicting"
    elif votes and len(votes) == len(markers):
        result.update(state="predicted" if any(b in {"predicted", "derived"} for b in bases) else "supported",
                      positive_compartment=votes[0], negative_compartment=compartments[1 - compartments.index(votes[0])])
    return result


def verify_placement_sidedness(atoms, membrane):
    """Re-map exported markers rather than trusting a stored compartment label."""
    stored = membrane.get("side_assignment", {"state": "unknown"})
    state = membrane.get("sidedness", stored.get("state", "unknown"))
    if state == "unknown":
        if stored.get("state", "unknown") != "unknown":
            raise ValueError("Unsigned placement has a conflicting signed record")
        return {**stored, "state": "unknown", "positive_compartment": None, "negative_compartment": None}
    if state not in {"supported", "predicted", "conflicting"} or not stored.get("markers"):
        raise ValueError("Signed placements require mapped markers")
    resolved = resolve_sidedness(atoms, np.asarray(membrane["center"]), np.asarray(membrane["normal"]),
                                membrane["half_thickness"], membrane["context"],
                                {"compartments": stored.get("compartments"), "side_markers": stored["markers"], "source": stored.get("source", {})})
    for key in ("state", "positive_compartment", "negative_compartment"):
        if resolved.get(key) != stored.get(key):
            raise ValueError("Stored sidedness conflicts with mapped coordinates")
    for old, new in zip(stored["markers"], resolved["markers"]):
        if not np.isclose(old.get("depth_A", np.nan), new["depth_A"], atol=1e-5, rtol=0):
            raise ValueError("Sidedness marker depth is in another coordinate frame")
    return resolved
