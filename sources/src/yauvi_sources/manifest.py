"""What a module declares it needs.

Each module ships a `sources.yaml` naming the registry entries it consumes and
what each one is for. That file is the module's half of the contract: the
registry says what a source *is*, the manifest says why *this* module wants it.

Keeping the requirement list in the module rather than in the registry is what
lets a module be installed and planned on its own. `yauvi-fetch plan --for
memorient` works with nothing else from the workspace present except the
registry itself.

Manifest shape:

    schema_version: "1.0"
    module_id: subproteo
    requires:
      - source_id: uniprot_proteomes
        role: "target and host proteomes"
        required: true
      - source_id: deg
        role: "essentiality reference"
        required: false
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import yaml


class ManifestError(RuntimeError):
    pass


@dataclass(frozen=True)
class Requirement:
    source_id: str
    role: str = ""
    required: bool = True


@dataclass(frozen=True)
class ModuleManifest:
    module_id: str
    requires: Sequence[Requirement]
    origin: str = ""

    def source_ids(self) -> list[str]:
        return [r.source_id for r in self.requires]

    def required_only(self) -> list[Requirement]:
        return [r for r in self.requires if r.required]


def parse_manifest(document: Mapping[str, object], *, origin: str = "") -> ModuleManifest:
    module_id = str(document.get("module_id") or "").strip()
    if not module_id:
        raise ManifestError(f"module manifest is missing module_id: {origin or '<inline>'}")

    raw_requires = document.get("requires")
    if raw_requires is None:
        raw_requires = []
    if not isinstance(raw_requires, list):
        raise ManifestError(f"'requires' must be a list in {origin or module_id}")

    requirements: list[Requirement] = []
    seen: set[str] = set()
    for item in raw_requires:
        if isinstance(item, str):
            item = {"source_id": item}
        if not isinstance(item, Mapping):
            raise ManifestError(f"malformed requirement in {origin or module_id}: {item!r}")
        source_id = str(item.get("source_id") or "").strip()
        if not source_id:
            raise ManifestError(f"requirement without source_id in {origin or module_id}")
        if source_id in seen:
            raise ManifestError(
                f"{module_id} declares source {source_id!r} twice; a requirement is stated once"
            )
        seen.add(source_id)
        requirements.append(
            Requirement(
                source_id=source_id,
                role=str(item.get("role", "")).strip(),
                required=bool(item.get("required", True)),
            )
        )

    return ModuleManifest(module_id=module_id, requires=tuple(requirements), origin=origin)


def load_manifest(path: str | Path) -> ModuleManifest:
    path = Path(path)
    if not path.is_file():
        raise ManifestError(f"module manifest not found: {path}")
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        raise ManifestError(f"module manifest is not valid YAML: {exc}") from exc
    if not isinstance(document, Mapping):
        raise ManifestError(f"module manifest must be a mapping: {path}")
    return parse_manifest(document, origin=str(path))


# Import names for the modules that ship a manifest, so `--for <id>` can find one
# in an installed environment without knowing where the workspace is.
KNOWN_MODULE_PACKAGES: Mapping[str, str] = {
    "subproteo": "subproteo",
    "subproteo-pipeline": "subproteo_pipeline",
    "memorient": "memorient",
    "sf_csa": "sf_csa",
    "sf-csa": "sf_csa",
    "actstate": "actstate",
    "activity_state": "actstate",
    "structqc": "structqc",
    "structure_quality": "structqc",
    "assembly_context": "assembly_context",
    "assembly-context": "assembly_context",
    "site_context": "site_context",
    "site-context": "site_context",
    "state_atlas": "state_atlas",
    "state-atlas": "state_atlas",
    "conformational_state": "state_atlas",
    "structcons": "structcons",
    "structural_conservation": "structcons",
    "oral_ecosystem": "oral_ecosystem",
    "oral-ecosystem": "oral_ecosystem",
}

# Where each module's manifest lives relative to the workspace root, for the
# common case of running inside a checkout before anything is pip-installed.
WORKSPACE_MANIFEST_PATHS: Mapping[str, str] = {
    "subproteo": "Subtractive Proteomics/src/subproteo/sources.yaml",
    "subproteo-pipeline": "Subtractive Proteomics/src/subproteo_pipeline/sources.yaml",
    "memorient": "Membrane Orientor/memorient/src/memorient/sources.yaml",
    "sf_csa": "sf-csa/src/sf_csa/sources.yaml",
    "sf-csa": "sf-csa/src/sf_csa/sources.yaml",
    "actstate": "activity-state/src/actstate/sources.yaml",
    "activity_state": "activity-state/src/actstate/sources.yaml",
    "structqc": "structqc/src/structqc/sources.yaml",
    "structure_quality": "structqc/src/structqc/sources.yaml",
    "assembly_context": "assembly-context/src/assembly_context/sources.yaml",
    "assembly-context": "assembly-context/src/assembly_context/sources.yaml",
    "site_context": "site-context/src/site_context/sources.yaml",
    "site-context": "site-context/src/site_context/sources.yaml",
    "state_atlas": "state-atlas/src/state_atlas/sources.yaml",
    "state-atlas": "state-atlas/src/state_atlas/sources.yaml",
    "conformational_state": "state-atlas/src/state_atlas/sources.yaml",
    "structcons": "structcons/src/structcons/sources.yaml",
    "structural_conservation": "structcons/src/structcons/sources.yaml",
    "oral_ecosystem": "oral-ecosystem/src/oral_ecosystem/sources.yaml",
    "oral-ecosystem": "oral-ecosystem/src/oral_ecosystem/sources.yaml",
}


def _from_installed_package(module_id: str) -> ModuleManifest | None:
    package = KNOWN_MODULE_PACKAGES.get(module_id)
    if not package:
        return None
    try:
        from importlib.resources import files  # noqa: PLC0415

        resource = files(package) / "sources.yaml"
        if not resource.is_file():
            return None
        document = yaml.safe_load(resource.read_text(encoding="utf-8")) or {}
    except (ImportError, ModuleNotFoundError, FileNotFoundError, OSError):
        return None
    except yaml.YAMLError as exc:
        raise ManifestError(f"installed manifest for {module_id} is not valid YAML: {exc}") from exc
    return parse_manifest(document, origin=f"<installed:{package}>")


def _from_workspace(module_id: str, workspace: Path) -> ModuleManifest | None:
    relative = WORKSPACE_MANIFEST_PATHS.get(module_id)
    if not relative:
        return None
    candidate = workspace / relative
    return load_manifest(candidate) if candidate.is_file() else None


def resolve_manifest(
    module_id: str,
    *,
    explicit_path: str | Path | None = None,
    workspace: str | Path | None = None,
) -> ModuleManifest:
    """Find a module's manifest: explicit path, then installed package, then workspace."""
    if explicit_path:
        return load_manifest(explicit_path)

    manifest = _from_installed_package(module_id)
    if manifest is not None:
        return manifest

    if workspace:
        manifest = _from_workspace(module_id, Path(workspace))
        if manifest is not None:
            return manifest

    known = ", ".join(sorted(set(KNOWN_MODULE_PACKAGES) | set(WORKSPACE_MANIFEST_PATHS)))
    raise ManifestError(
        f"no source manifest found for module {module_id!r}.\n"
        f"  Pass --manifest <path>, install the module, or run from the workspace root.\n"
        f"  Modules with a known manifest location: {known}"
    )
