"""Portable tester launcher. Uses a separate environment and local analysis library."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import venv


def verify(root):
    manifest = json.loads((root / "PACKAGE_MANIFEST.json").read_text())
    for entry in manifest["files"]:
        path = root / entry["path"]
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise RuntimeError("Package integrity check failed for " + entry["path"] + ". Request a fresh copy.")
    return manifest


def main():
    root = Path(__file__).resolve().parent
    verify(root)
    if "--check" in sys.argv:
        print("Package file checksums match.")
        return 0
    if sys.version_info < (3, 10):
        raise RuntimeError("Python 3.10 or newer is required. Python 3.12 is recommended for this preview.")
    runtime = root / ".runtime"
    python = runtime / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    marker = runtime / "installation-complete"
    wheel = next((root / "application").glob("*.whl"))
    wheel_hash = hashlib.sha256(wheel.read_bytes()).hexdigest()
    if not python.exists() or not marker.exists() or marker.read_text() != wheel_hash:
        print("First-time setup installs YAUVI and its Python dependencies into this package's .runtime folder.")
        print("Internet is required to download dependencies from PyPI. Your protein files are not uploaded.")
        print("For installation-free testing, cancel and open START_HERE.html instead.")
        if input("Install the local workbench now? [y/N] ").strip().lower() != "y":
            return 0
        if not python.exists():
            venv.EnvBuilder(with_pip=True).create(runtime)
        subprocess.run([str(python), "-m", "pip", "install", str(wheel)], check=True)
        marker.write_text(wheel_hash)
    workspace = Path.home() / "YAUVI Tester Library"
    workspace.mkdir(exist_ok=True)
    env = os.environ.copy()
    env["YAUVI_SHOWCASE_DIR"] = str(root / "showcase/workbench")
    # Ensure scientific child commands use the same environment as the workbench.
    env["PATH"] = str(python.parent) + os.pathsep + env.get("PATH", "")
    env.pop("PYTHONPATH", None)
    print("Your analyses will be saved in:", workspace)
    print("The workbench opens in your browser at http://127.0.0.1:8962 .")
    print("Keep this window open. Press Ctrl+C here to stop the server.")
    return subprocess.call([str(python), "-m", "yauvi_structural_workbench.cli",
                            "--workspace", str(workspace), "workbench", "open",
                            "--port", "8962", "--label", "YAUVI Tester Preview"], env=env)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, OSError, subprocess.CalledProcessError, ValueError) as error:
        print("Could not start:", error)
        print("Save this message for feedback. You can still use START_HERE.html offline.")
        input("Press Enter to close.")
        raise SystemExit(1)
