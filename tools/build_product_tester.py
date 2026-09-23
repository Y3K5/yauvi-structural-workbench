"""Assemble a local tester candidate from an explicit allowlist; never publish."""
from __future__ import annotations
import argparse
import hashlib
import html
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parent))
import content_screen  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NAME = "yauvi-product-tester-2026-09-19"
PACKAGES = {
    "structural-workbench": "yauvi_structural_workbench",
    "platform": "yauvi_platform/structural_workbench",
    "sources": "yauvi_sources", "structqc": "structqc",
    "Membrane Orientor/memorient": "memorient", "state-atlas": "state_atlas",
    "site-context": "site_context", "activity-state": "actstate",
    "assembly-context": "assembly_context", "sf-csa": "sf_csa", "structprep": "structprep",
}
TEMPLATE_FILES = ["START_HERE.html", "start_workbench.py", "Launch Mac.command", "Launch Linux.sh",
                  "Launch Windows.cmd", "VALIDATION.md", "THIRD_PARTY_NOTICES.md"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def screen(files, denylist, denylist_source=""):
    """Inspect the actual resolved payload, including wheel members; report rule names/counts only."""
    return content_screen.screen(files, denylist, denylist_source)


def add_denylist_arguments(parser):
    parser.add_argument("--denylist", help="Private JSON list of terms kept outside the repository; never packaged or reported.")
    parser.add_argument("--no-private-denylist", action="store_true",
                        help="Build without a private denylist; generic path and credential rules still apply.")


def denylist_from(args):
    try:
        return content_screen.load_denylist(args.denylist, allow_missing=args.no_private_denylist)
    except content_screen.DenylistError as exc:
        raise SystemExit(f"Candidate held: {exc}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable, help="Interpreter with local build dependencies")
    parser.add_argument("--out", default=str(ROOT / "build/product-tester"))
    add_denylist_arguments(parser)
    args = parser.parse_args()
    denylist, denylist_source = denylist_from(args)
    out = Path(args.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="yauvi-package-") as temp:
        stage = Path(temp) / "source"; stage.mkdir()
        for relative in ["pyproject.toml", "yauvi-structural-workbench/README.md", "yauvi-structural-workbench/LICENSE"]:
            target = stage / relative; target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / relative, target)
        for project, package in PACKAGES.items():
            relative = Path("software") / project / "src" / package
            shutil.copytree(ROOT / relative, stage / relative,
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".DS_Store", "tests"))
        wheels = Path(temp) / "wheels"
        subprocess.run([args.python, "-m", "build", "--wheel", "--no-isolation", "--outdir", str(wheels)], cwd=stage, check=True)
        wheel = next(wheels.glob("*.whl"))
        subprocess.run([args.python, str(ROOT / "tools/verify_structural_workbench_wheel.py"), str(wheel)], check=True)
        payload = {"application/" + wheel.name: wheel.read_bytes()}
        for name in TEMPLATE_FILES:
            payload[name] = (ROOT / "tools/tester-kit" / name).read_bytes()
        for name in ["index.html", "workbench/index.html", "workbench/showcase.js", "ABL-endpoint-geometry.json"]:
            payload["showcase/" + name] = (ROOT / "build/showcase-3d" / name).read_bytes()
        payload["START_HERE.html"] = payload["START_HERE.html"].decode().replace(
            "__EMBEDDED_SHOWCASE__", html.escape(payload["showcase/index.html"].decode(), quote=True)).encode()
        payload["LICENSE"] = (ROOT / "yauvi-structural-workbench/LICENSE").read_bytes()
        vendor = ROOT / "software/structural-workbench/src/yauvi_structural_workbench/ui/vendor"
        payload["3DMOL-LICENSE.txt"] = (vendor / "3DMOL-LICENSE.txt").read_bytes()
        audit = screen(payload, denylist, denylist_source)
        (out / "CONTENT_SCREEN.json").write_text(json.dumps(audit, indent=2) + "\n")
        if audit["findings"]:
            print(json.dumps(audit, indent=2))
            raise SystemExit("Candidate held: content findings require review; no ZIP created.")
        payload["CONTENT_SCREEN.json"] = (json.dumps(audit, indent=2) + "\n").encode()
        manifest = {"schema_version": 1, "package": NAME, "channel": "development_tester_candidate",
                    "purpose": "Guided usability pilot; conditional own-input scientific execution",
                    "standalone_desktop_application": False,
                    "external_distribution": "Prepared locally; no sending, upload or publication performed",
                    "files": [{"path": n, "bytes": len(v), "sha256": sha(v)} for n, v in sorted(payload.items())]}
        payload["PACKAGE_MANIFEST.json"] = (json.dumps(manifest, indent=2) + "\n").encode()
        destination = out / NAME; destination.mkdir(exist_ok=True)
        archive_path = out / (NAME + ".zip")
        with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for name, data in sorted(payload.items()):
                path = destination / name; path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(data)
                mode = 0o755 if name.endswith((".command", ".sh")) else 0o644
                path.chmod(mode)
                info = zipfile.ZipInfo(NAME + "/" + name, (2026, 9, 19, 0, 0, 0))
                info.external_attr = (0o100000 | mode) << 16
                info.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(info, data)
        # Re-open the ZIP so the manifest check applies to what will actually be shared.
        with zipfile.ZipFile(archive_path) as archive:
            assert archive.testzip() is None
            for entry in manifest["files"]:
                assert sha(archive.read(NAME + "/" + entry["path"])) == entry["sha256"]
        digest = sha(archive_path.read_bytes())
        (out / (NAME + ".zip.sha256")).write_text(digest + "  " + archive_path.name + "\n")
        print(json.dumps({"archive": str(archive_path), "sha256": digest, "bytes": archive_path.stat().st_size,
                          "payload_files": len(payload), "screened_files": audit["entries_screened_including_archive_members"],
                          "content_findings": 0}, indent=2))


if __name__ == "__main__":
    main()
