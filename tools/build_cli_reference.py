#!/usr/bin/env python3
"""Generate docs/cli-reference.md from each console script's own --help.

Run with the distribution installed and its environment active:

    python tools/build_cli_reference.py

Generating from --help rather than hand-writing the reference is deliberate:
the document cannot drift away from the parsers it documents.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

CLIS: list[tuple[str, str, list[str]]] = [
    ("yauvi", "Structural analysis case store: create, add inputs, validate, run, export.",
     ["analysis create", "analysis add", "analysis validate", "analysis run", "analysis export"]),
    ("structqc", "Coordinate trust: completeness, provenance class, imported validation.",
     ["describe", "validate", "fetch", "run"]),
    ("memorient", "Membrane orientation and sidedness labelling.",
     ["contexts", "describe", "orient", "run", "validate", "fetch"]),
    ("state-atlas", "Conformational-state resemblance against declared references.",
     ["describe", "fetch", "validate", "run"]),
    ("site-context", "Functional-site roles, cofactors, ligands, and pockets.",
     ["describe", "fetch", "validate", "run"]),
    ("actstate", "Activity-state evidence assembly.",
     ["describe", "fetch", "validate", "run"]),
    ("assembly-context", "Biological assembly, stoichiometry, contacts, and burial.",
     ["describe", "fetch", "validate", "run"]),
    ("sf-csa", "Structure- and sequence-based functional comparison.",
     ["describe", "validate", "fetch", "run", "verify", "build-manifests"]),
    ("yauvi-fetch", "Registered public-source acquisition and staging.",
     ["sources", "plan", "get", "stage", "verify", "where"]),
]

PREAMBLE = """# CLI reference

Generated from each command's own `--help`, so it cannot drift from the code.
Regenerate with `python tools/build_cli_reference.py`.

Common conventions across the scientific modules:

- `describe` prints the module contract as JSON: inputs, outputs, and the claim ceiling.
- `validate` checks inputs without producing an evidence record.
- `run` performs the analysis and writes a result bundle to `--out`.
- `fetch` resolves registered public sources. Acquisition never adopts a file
  into an analysis automatically, and stays disabled without an explicit
  reference-fetch flag.
- Exit codes: `0` completed, `1` scientifically incomplete, `2` invalid input or
  configuration. `1` is a real result, not a failure to be retried around.
"""


def capture(args: list[str]) -> str | None:
    exe = shutil.which(args[0])
    if exe is None:
        return None
    try:
        done = subprocess.run([exe, *args[1:]], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    text = (done.stdout or done.stderr).rstrip()
    if not text:
        return None
    # Some parsers print a default path under the caller's home directory.
    # Generated documentation must not embed whoever ran the generator.
    home = str(Path.home())
    if home and home != "/":
        text = text.replace(home, "~")
    return text


def build() -> tuple[str, list[str]]:
    missing: list[str] = []
    out = [PREAMBLE, "\n| Command | Purpose |\n|---|---|"]
    for name, purpose, _ in CLIS:
        out.append(f"| [`{name}`](#{name}) | {purpose} |")
    out.append("")
    for name, purpose, subs in CLIS:
        top = capture([name, "--help"])
        if top is None:
            missing.append(name)
            continue
        out.append(f"\n---\n\n## {name}\n\n{purpose}\n")
        out.append(f"```\n{top}\n```\n")
        for sub in subs:
            body = capture([name, *sub.split(), "--help"])
            if body:
                out.append(f"### `{name} {sub}`\n")
                out.append(f"```\n{body}\n```\n")
    return "\n".join(out) + "\n", missing


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path,
        default=Path("yauvi-structural-workbench/docs/cli-reference.md"),
        help="Destination markdown file.",
    )
    parser.add_argument(
        "--check", action="store_true",
        help="Exit 1 if the generated reference differs from the file on disk.",
    )
    args = parser.parse_args(argv)

    text, missing = build()
    if missing:
        print(
            "not on PATH, so no reference was generated for: " + ", ".join(missing)
            + "\ninstall the distribution first: python -m pip install -e \".[dev]\"",
            file=sys.stderr,
        )
        return 2
    if args.check:
        current = args.out.read_text(encoding="utf-8") if args.out.exists() else ""
        if current != text:
            print(f"{args.out} is stale; regenerate it.", file=sys.stderr)
            return 1
        print(f"{args.out} is current.")
        return 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text, encoding="utf-8")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
