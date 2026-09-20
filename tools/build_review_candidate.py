#!/usr/bin/env python3
"""Build a local review archive from an explicit, checksum-bound file manifest.

Never discovers additional files, copies a whole workspace, publishes, or modifies
source evidence. A failed screen writes only a findings report and no archive.
Optional private denylist terms stay outside the exported bundle and findings.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import stat
import zipfile


PATTERNS = {
    "absolute_home_path": re.compile(rb"/(?:Users|home)/[A-Za-z0-9_.-]+/"),
    "github_token": re.compile(rb"(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})"),
    "provider_secret": re.compile(rb"sk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{40,}"),
    "private_key": re.compile(rb"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
}
FORBIDDEN_PARTS = {".git", ".claude", ".env", ".venv", "__pycache__", "personal"}
MAX_EXPANDED_BYTES = 128 * 1024 * 1024


def safe_relative(name: str) -> PurePosixPath:
    path = PurePosixPath(name)
    if (not name or path.is_absolute() or "\\" in name or
            any(p in {".", ".."} for p in name.split("/")) or
            any(p in FORBIDDEN_PARTS or p.startswith(".env.") for p in path.parts)):
        raise ValueError("unsafe or excluded relative path")
    if len(path.parts) > 3 and path.parts[:2] == ("evidence", "benchmarks") and "sources" in path.parts:
        raise ValueError("acquired benchmark sources are not redistributable")
    return path


def screen(name: str, data: bytes, denylist: list[str], depth: int = 0,
           budget: list[int] | None = None) -> list[dict]:
    budget = budget if budget is not None else [MAX_EXPANDED_BYTES]
    budget[0] -= len(data)
    if budget[0] < 0 or depth > 3:
        return [{"file": name, "rule": "archive_inspection_limit", "count": 1}]
    screened_bytes = name.encode() + b"\n" + data
    findings = [{"file": name, "rule": rule, "count": len(pattern.findall(screened_bytes))}
                for rule, pattern in PATTERNS.items() if pattern.search(screened_bytes)]
    for index, term in enumerate(denylist):
        count = screened_bytes.lower().count(term.encode().lower())
        if term and count:
            findings.append({"file": name, "rule": f"private_term_{index + 1}", "count": count})
    if zipfile.is_zipfile(io.BytesIO(data)):
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            seen = set()
            for item in archive.infolist():
                if item.is_dir():
                    continue
                child = name + "!" + item.filename
                try:
                    safe_relative(item.filename)
                    if item.filename in seen or stat.S_ISLNK(item.external_attr >> 16):
                        raise ValueError("duplicate or symlink archive member")
                    seen.add(item.filename)
                    if item.file_size > budget[0]:
                        raise ValueError("expanded archive too large")
                    content = archive.read(item)
                except (ValueError, RuntimeError, zipfile.BadZipFile):
                    findings.append({"file": child, "rule": "unsafe_archive_member", "count": 1})
                    continue
                findings.extend(screen(child, content, denylist, depth + 1, budget))
    elif name.endswith((".gz", ".tgz", ".tar", ".bz2", ".xz", ".7z")):
        findings.append({"file": name, "rule": "unsupported_archive_requires_review", "count": 1})
    return findings


def build(source: Path, manifest: dict, output: Path, denylist: list[str]) -> dict:
    source = source.resolve()
    output = output.resolve()
    if output == source or source in output.parents:
        raise ValueError("review outputs must live outside the source tree")
    if output.exists() and any(output.iterdir()):
        raise ValueError("output directory must be empty; preserve prior results")
    records = manifest.get("files", [])
    if not records:
        raise ValueError("manifest must name at least one file")
    payloads = {}
    findings = []
    for record in records:
        name = record["path"]
        if name == "REVIEW_MANIFEST.json":
            raise ValueError("reserved archive manifest filename")
        relative = safe_relative(name)
        path = source / relative
        if name in payloads:
            raise ValueError("duplicate manifest entry")
        if not path.is_file() or path.is_symlink() or source not in path.resolve().parents:
            raise ValueError(f"missing, symlinked or escaping file: {name}")
        if any((source.joinpath(*relative.parts[:i])).is_symlink() for i in range(1, len(relative.parts))):
            raise ValueError(f"symlinked parent: {name}")
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != record["sha256"]:
            raise ValueError(f"source changed since manifest review: {name}")
        payloads[name] = data
        findings.extend(screen(name, data, denylist))
    manifest_bytes = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    findings.extend(screen("REVIEW_MANIFEST.json", manifest_bytes, denylist))
    reviewed_findings = []
    seen_exceptions = set()
    for exception in manifest.get("reviewed_exceptions", []):
        name, rule = exception["path"], exception["rule"]
        key = (name, rule)
        if (key in seen_exceptions or rule != "absolute_home_path"
                or name not in payloads or not exception.get("reason")
                or not exception.get("approved_by")):
            raise ValueError("invalid or duplicate reviewed exception")
        seen_exceptions.add(key)
        if hashlib.sha256(payloads[name]).hexdigest() != exception["sha256"]:
            raise ValueError("reviewed exception is stale; exact file changed")
        matches = [f for f in findings if f["file"] == name and f["rule"] == rule]
        if len(matches) != 1 or matches[0]["count"] != exception["count"]:
            raise ValueError("reviewed exception finding count changed")
        reviewed_findings.append({**matches[0], "sha256": exception["sha256"],
                                  "reason": exception["reason"],
                                  "approved_by": exception["approved_by"]})
        findings.remove(matches[0])
    output.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": "1.0", "file_count": len(payloads),
        "screen_passed": not findings, "findings": findings,
        "reviewed_findings": reviewed_findings,
        "authorization": "local_review_only_no_publication_authorization",
        "coverage": "Exact manifest bytes, file paths, and nested ZIP members; no upload performed.",
        "limitation": "Pattern screening supplements human review; it cannot establish absence of all sensitive meaning.",
    }
    if not findings:
        archive_path = output / "yauvi-review-candidate.zip"
        with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for name, data in sorted({**payloads, "REVIEW_MANIFEST.json": manifest_bytes}.items()):
                info = zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                archive.writestr(info, data)
        report["archive"] = archive_path.name
        report["archive_sha256"] = hashlib.sha256(archive_path.read_bytes()).hexdigest()
    (output / "SCREEN_REPORT.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--denylist", type=Path, help="Private JSON list of terms; never exported.")
    args = parser.parse_args()
    try:
        denylist = json.loads(args.denylist.read_text()) if args.denylist else []
        if not isinstance(denylist, list) or any(not isinstance(x, str) or not x for x in denylist):
            raise ValueError("denylist must be a list of nonempty strings")
        report = build(args.source, json.loads(args.manifest.read_text()), args.output_dir, denylist)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(2, f"candidate build refused: {exc}\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["screen_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
