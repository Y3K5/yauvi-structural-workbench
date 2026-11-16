"""Workbench adapter for MemOrient's domain-aware orientation estimate.

What goes in: one local PDB structure (and optionally a chain and a number of
leading residues to ignore for the side vote).
What happens: the run calls ``memorient domain-orient`` in a separate process,
with the checksum-verified copy of the structure.
What comes out (in ``outputs/memorient/``): ``DOMAIN_ORIENTATION.json``,
``RESIDUE_ZONES.tsv`` and, when a membrane core is found,
``ORIENTED_STRUCTURE.pdb`` + ``MEMBRANE_LAYER.json`` in the names the Workbench
viewer already reads, so the oriented model and membrane slab show up there.
The report is built by the store's normal report path from those files.
"""
from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path
from typing import Any

TYPE = "domain_membrane_orientation"
TYPES = {TYPE}
TITLE = "Membrane orientation · domain-aware estimate"
CLAIM = ("Estimated membrane pose from C-alpha geometry and hydropathy only; not a measured membrane position, "
         "and not evidence of native sidedness, exposure or function.")
_SOURCE_FILES = ("domain_orient.py", "domain_report.py", "cli.py")


def definitions() -> list[dict[str, Any]]:
    return [{
        "analysis_type": TYPE, "title": TITLE,
        "question": "Where is the membrane-embedded core of this structure, and which way does it sit in the membrane?",
        "module_ids": ["membrane_orientation"], "readiness": "experimental", "claim_ceiling": CLAIM,
        "use_when": "Estimate a membrane frame for barrels with soluble domains or apolar helix bundles from a local PDB.",
        "measures": "Membrane core found (barrel or helices), normal and depth, side votes, residue heights and zones.",
        "receives": ["Oriented model in the viewer", "Membrane slab overlay", "Residue zone table", "Side-vote record"],
        "non_claim": CLAIM,
        "scientific_readiness": {"software_state": "experimental", "external_benchmark": "not_evaluated",
                                 "release_gate": "independent_reference_comparison_required"},
        "inputs": [{
            "role": "structure", "label": "Local PDB structure (first model is used)", "required": True, "multiple": False,
            "extensions": [".pdb", ".ent"], "accepted_extensions": [".pdb", ".ent"], "template_id": None,
            "description": "Coordinates whose C-alpha atoms are used to find a membrane core.",
            "why_needed": "The estimate is computed from these exact coordinates.",
            "absence_effect": "blocked", "accepted_artifact_types": [], "source_ids": [],
            "format_guide": "PDB format; only the first model and standard amino acids are read.",
            "validator_id": "structure", "sensitivity": "sensitive_by_default",
        }, {
            # Same optional coordinate-validation input every structure workflow offers; when supplied,
            # the store's preflight checks it is bound to this exact structure. The estimate does not use it.
            "role": "validation_report", "label": "Optional coordinate validation", "required": False, "multiple": False,
            "extensions": [".json", ".xml"], "accepted_extensions": [".json", ".xml"], "template_id": None,
            "description": "A checksum-bound wwPDB, MolProbity, or Phenix validation result for the same coordinates.",
            "why_needed": "Records coordinate quality next to the estimate; it does not change the orientation.",
            "absence_effect": "unevaluated", "accepted_artifact_types": ["wwpdb.validation"], "source_ids": ["wwpdb_validation"],
            "format_guide": "wwPDB validation XML or a supported normalized local validation JSON.",
            "validator_id": "structure_validation_report", "sensitivity": "sensitive_by_default",
        }],
        "parameters": [
            {"name": "chain", "label": "Chain (blank = all chains)", "type": "text", "default": ""},
            {"name": "skip_n", "label": "Leading residues to ignore for the side vote", "type": "integer", "default": 0},
        ],
    }]


def _memorient_root() -> Path | None:
    spec = importlib.util.find_spec("memorient")
    return Path(spec.origin).parent if spec and spec.origin else None


def source_identity() -> dict[str, str | None]:
    root = _memorient_root()
    if root is None:
        return {"memorient_domain_orient": None}
    digest = hashlib.sha256()
    for name in _SOURCE_FILES:
        digest.update(name.encode() + b"\0" + (root / name).read_bytes() + b"\0")
    return {"memorient_domain_orient": digest.hexdigest()}


def readiness() -> list[dict[str, Any]]:
    try:
        identity = source_identity()["memorient_domain_orient"]
        state, detail = ("experimental", CLAIM) if identity else ("blocked_missing_runtime", "MemOrient is unavailable")
    except OSError as exc:
        identity, state, detail = None, "blocked_source_unavailable", str(exc)
    return [{"contract_id": "tool_readiness_record", "analysis_type": TYPE, "title": TITLE, "state": state,
             "package_source_sha256": identity, "required_runtimes": {}, "optional_runtimes": {},
             "external_benchmark": {"state": "not_evaluated", "collection": "none"},
             "scientific_scopes": [], "limitations": [detail]}]


def _settings(manifest) -> tuple[str, int]:
    params = manifest.get("parameters", {}) or {}
    chain = str(params.get("chain", "") or "")
    skip_n = params.get("skip_n", 0)
    if len(chain) > 1:
        raise ValueError("Chain must be one character or blank")
    if type(skip_n) is not int or not 0 <= skip_n <= 2000:
        raise ValueError("Leading residues to ignore must be a whole number from 0 to 2000")
    return chain, skip_n


def preflight(store, manifest, inputs) -> dict[str, Any]:
    from memorient.domain_orient import read_pdb_ca
    chain, _ = _settings(manifest)
    ca, _seq, _num = read_pdb_ca(str(store._one(inputs, "structure")), chain or None)
    if len(ca) < 3:
        raise ValueError("No usable C-alpha atoms" + (f" in chain {chain!r}" if chain else ""))
    return {"c_alpha_residues": int(len(ca))}


def execute(store, manifest, inputs, run_dir, cancel_event=None, on_process=None, on_progress=None):
    chain, skip_n = _settings(manifest)
    output = run_dir / "outputs" / "memorient"
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, "-m", "memorient.cli", "domain-orient", str(store._one(inputs, "structure")),
               "--out", str(output), "--skip-n", str(skip_n)]
    if chain:
        command.extend(["--chain", chain])
    code = store._execute(command, cwd=run_dir, log_path=run_dir / "logs" / "domain-orient.log",
                          cancel_event=cancel_event, on_process=on_process)
    # 0 = pose written; 1 = no membrane core found (scientifically incomplete); anything else = failed.
    code = code if code in (0, 1) else 2
    return [{"module_id": TYPE, "exit_code": code}], code
