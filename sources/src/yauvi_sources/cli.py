"""`yauvi-fetch` — plan, acquire, and verify a module's raw input files.

    yauvi-fetch plan   --for subproteo          what is needed, what is here
    yauvi-fetch get    --for subproteo          retrieve what policy permits
    yauvi-fetch stage  deg <path>               adopt a manually acquired file
    yauvi-fetch verify [--source-id ID]         re-hash the cache
    yauvi-fetch where  uniprot_proteomes        print the cached path
    yauvi-fetch sources [--channel localization]  list the registry

Exit codes carry meaning, because these commands are run from scripts:
    0  the request was satisfied
    1  a required source is absent and a human must act
    2  usage or configuration error
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from typing import Sequence

from .cache import CacheError, SourceCache, default_cache_dir
from .manifest import ManifestError, resolve_manifest
from .planner import Plan, build_plan, render_plan
from .policy import FetchClass, PolicyError, classify, instructions_for
from .registry import RegistryError, SourceRegistry

EXIT_OK = 0
EXIT_UNSATISFIED = 1
EXIT_USAGE = 2

REGISTRY_ENV = "YAUVI_SOURCES_REGISTRY"
REGISTRY_RELATIVE = Path("catalogs") / "sources.yaml"


# -- locating the registry ------------------------------------------------


def find_registry(explicit: str | None = None, start: Path | None = None) -> Path:
    """Locate catalogs/sources.yaml: flag, then env, then upward from cwd."""
    if explicit:
        return Path(explicit).expanduser()
    from_env = os.environ.get(REGISTRY_ENV)
    if from_env:
        return Path(from_env).expanduser()
    here = (start or Path.cwd()).resolve()
    for directory in (here, *here.parents):
        candidate = directory / REGISTRY_RELATIVE
        if candidate.is_file():
            return candidate
    raise RegistryError(
        f"could not locate {REGISTRY_RELATIVE}. Pass --registry, set {REGISTRY_ENV}, "
        f"or run from inside the workspace."
    )


def workspace_of(registry_path: Path) -> Path:
    """The workspace root implied by a registry at <root>/catalogs/sources.yaml."""
    return registry_path.resolve().parent.parent


def _parse_args_pairs(pairs: Sequence[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for pair in pairs or ():
        if "=" not in pair:
            raise ValueError(f"--arg expects source_id=value, got {pair!r}")
        key, value = pair.split("=", 1)
        out[key.strip()] = value.strip()
    return out


# -- commands -------------------------------------------------------------


def _load_plan(args) -> tuple[Plan, SourceCache, SourceRegistry]:
    registry_path = find_registry(args.registry)
    registry = SourceRegistry.load(registry_path)
    manifest = resolve_manifest(
        args.module,
        explicit_path=args.manifest,
        workspace=workspace_of(registry_path),
    )
    cache = SourceCache(args.cache)
    return build_plan(manifest, registry, cache), cache, registry


def cmd_plan(args) -> int:
    plan, _, _ = _load_plan(args)
    if args.json:
        print(json.dumps(_plan_as_dict(plan), indent=2, sort_keys=True))
    else:
        print(render_plan(plan, verbose=args.verbose))
    return EXIT_OK if plan.satisfied else EXIT_UNSATISFIED


def _plan_as_dict(plan: Plan) -> dict:
    return {
        "module_id": plan.module_id,
        "satisfied": plan.satisfied,
        "items": [
            {
                "source_id": item.source_id,
                "display_name": item.source.display_name,
                "required": item.requirement.required,
                "role": item.requirement.role,
                "fetch_class": item.fetch_class.value,
                "status": item.status,
                "blocks": item.blocks,
                "license_note": item.source.license_note,
                "cached": (
                    {
                        "sha256": item.cached.sha256,
                        "bytes": item.cached.bytes,
                        "retrieved_at": item.cached.retrieved_at,
                        "version": item.cached.version,
                        "origin": item.cached.origin,
                    }
                    if item.cached
                    else None
                ),
            }
            for item in plan.items
        ],
    }


def cmd_get(args) -> int:
    # lazy import: keeps `plan` usable without the network extra installed
    from .fetchers import FETCH_HOSTS, NAMED_FETCHERS, fetch_url, host_reachable

    plan, cache, _ = _load_plan(args)
    fetch_args = _parse_args_pairs(args.arg)

    targets = [i for i in plan.fetchable() if not args.source_id or i.source_id == args.source_id]
    if not targets:
        print("Nothing to fetch: every fetchable source is already cached.")
    failures: list[str] = []

    # One bounded probe before attempting anything. Without it an offline machine
    # stalls on every source in turn; with it the run says so immediately and
    # still prints the manual instructions, which are useful offline.
    if targets and not args.no_probe:
        if not any(host_reachable(host) for host in FETCH_HOSTS):
            print(
                "No network route to any source endpoint "
                f"({', '.join(FETCH_HOSTS)}).\n"
                "Nothing was downloaded. Sources already in the cache remain usable;\n"
                "run `yauvi-fetch plan` to see what is present, or stage files by hand.\n"
                "Use --no-probe to attempt retrieval anyway.",
                file=sys.stderr,
            )
            return EXIT_UNSATISFIED

    for item in targets:
        source_id = item.source_id
        fetcher = NAMED_FETCHERS.get(source_id)
        argument = fetch_args.get(source_id)

        if fetcher is not None:
            if not argument:
                failures.append(
                    f"{source_id}: needs an identifier — pass --arg {source_id}=<value> "
                    f"(config key: {item.source.config_key or 'see registry'})"
                )
                continue
            outcome = fetcher(argument)
        elif item.source.url:
            outcome = fetch_url(item.source.url)
        else:
            failures.append(
                f"{source_id}: the registry entry names no download URL and there is no "
                f"named fetcher for it. Add a `url:` to catalogs/sources.yaml or stage it by hand."
            )
            continue

        if not outcome.ok:
            failures.append(f"{source_id}: {outcome.reason}")
            print(f"  [miss] {source_id}: {outcome.reason}")
            continue

        entry = cache.store(
            source_id,
            outcome.payload,
            filename=outcome.filename,
            origin=outcome.origin,
            version=outcome.version,
        )
        print(
            f"  [ok]   {source_id}: {entry.filename} "
            f"({entry.bytes} bytes, sha256 {entry.sha256[:12]}"
            + (f", release {entry.version}" if entry.version else "")
            + ")"
        )
        _record_provenance(args.run_dir, entry)

    # Anything a human must acquire is reported here, in full, and never attempted.
    manual = plan.manual()
    if manual:
        print()
        print("Not retrieved — these require you to act:")
        for item in manual:
            print()
            print(f"  {item.source_id} — {item.source.display_name}")
            for line in instructions_for(item.source).splitlines():
                print(f"    {line}")

    if failures:
        print()
        print(f"{len(failures)} source(s) not acquired:")
        for failure in failures:
            print(f"  - {failure}")

    still_blocking = manual or failures
    return EXIT_UNSATISFIED if still_blocking else EXIT_OK


def _record_provenance(run_dir: str | None, entry) -> None:
    """Append the acquisition to the platform run ledger when one is in use.

    Soft dependency by design: yauvi-sources must install and work without the
    platform. When no ledger is available the cache manifest remains the record.
    """
    if not run_dir:
        return
    try:
        from yauvi_platform.runstore import AppendOnlyRunStore  # noqa: PLC0415
    except (ImportError, ModuleNotFoundError):
        print(
            f"  [warn] --run-dir given but yauvi_platform is not installed; "
            f"acquisition of {entry.source_id} recorded in the cache manifest only",
            file=sys.stderr,
        )
        return
    AppendOnlyRunStore(run_dir).append("source_acquisition", entry)


def cmd_stage(args) -> int:
    registry_path = find_registry(args.registry)
    registry = SourceRegistry.load(registry_path)
    source = registry.get(args.source_id)  # raises on an undeclared source
    cache = SourceCache(args.cache)
    entry = cache.stage(args.source_id, args.path, note=args.note or "")
    print(
        f"staged {source.display_name} as {entry.filename} "
        f"({entry.bytes} bytes, sha256 {entry.sha256})"
    )
    print(f"  -> {cache.path_for(entry)}")
    _record_provenance(args.run_dir, entry)
    return EXIT_OK


def cmd_verify(args) -> int:
    cache = SourceCache(args.cache)
    checked = 0
    bad = 0
    for entry, ok, detail in cache.verify(args.source_id):
        checked += 1
        if not ok:
            bad += 1
            print(f"  [FAIL] {entry.source_id}/{entry.filename}: {detail}")
        elif args.verbose:
            print(f"  [ok]   {entry.source_id}/{entry.filename}")
    if checked == 0:
        print(f"cache is empty: {cache.root}")
        return EXIT_OK
    print(f"{checked} cached file(s) checked, {bad} failed.")
    return EXIT_UNSATISFIED if bad else EXIT_OK


def cmd_where(args) -> int:
    cache = SourceCache(args.cache)
    entry = cache.latest(args.source_id)
    if entry is None:
        print(f"{args.source_id} is not cached in {cache.root}", file=sys.stderr)
        return EXIT_UNSATISFIED
    print(cache.path_for(entry))
    return EXIT_OK


def cmd_sources(args) -> int:
    registry = SourceRegistry.load(find_registry(args.registry))
    rows = [s for s in registry if not args.channel or s.channel == args.channel]
    rows.sort(key=lambda s: (s.channel, s.source_id))
    if args.json:
        print(
            json.dumps(
                [
                    {
                        "source_id": s.source_id,
                        "display_name": s.display_name,
                        "channel": s.channel,
                        "kind": s.kind,
                        "status": s.status,
                        "access": s.access,
                        "fetch_class": classify(s).value,
                        "license_note": s.license_note,
                    }
                    for s in rows
                ],
                indent=2,
                sort_keys=True,
            )
        )
        return EXIT_OK
    width = max((len(s.source_id) for s in rows), default=10)
    channel = None
    for source in rows:
        if source.channel != channel:
            channel = source.channel
            print(f"\n{channel}")
        print(
            f"  {source.source_id:<{width}}  {classify(source).value:<14} "
            f"{source.status:<20} {source.display_name}"
        )
    print(f"\n{len(rows)} source(s).")
    return EXIT_OK


# -- argument parsing -----------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="yauvi-fetch",
        description="Plan, acquire, and verify the raw input files a module declares.",
    )
    parser.add_argument("--registry", help="path to catalogs/sources.yaml")
    parser.add_argument(
        "--cache",
        help=f"source cache directory (default: {default_cache_dir()})",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def add_module_args(p):
        p.add_argument("--for", dest="module", required=True, help="module id, e.g. subproteo")
        p.add_argument("--manifest", help="explicit path to the module's sources.yaml")

    p_plan = sub.add_parser("plan", help="report what a module needs and what is present")
    add_module_args(p_plan)
    p_plan.add_argument("--json", action="store_true", help="machine-readable output")
    p_plan.add_argument("-v", "--verbose", action="store_true")
    p_plan.set_defaults(func=cmd_plan)

    p_get = sub.add_parser("get", help="retrieve the sources policy permits")
    add_module_args(p_get)
    p_get.add_argument("--source-id", help="fetch only this source")
    p_get.add_argument(
        "--arg",
        action="append",
        default=[],
        metavar="SOURCE_ID=VALUE",
        help="identifier for a source that needs one, e.g. uniprot_proteomes=UP000005640",
    )
    p_get.add_argument("--run-dir", help="append acquisitions to this platform run ledger")
    p_get.add_argument(
        "--no-probe",
        action="store_true",
        help="skip the up-front reachability check and attempt every retrieval",
    )
    p_get.set_defaults(func=cmd_get)

    p_stage = sub.add_parser("stage", help="adopt a file you acquired by hand")
    p_stage.add_argument("source_id")
    p_stage.add_argument("path")
    p_stage.add_argument("--note", help="how it was obtained")
    p_stage.add_argument("--run-dir")
    p_stage.set_defaults(func=cmd_stage)

    p_verify = sub.add_parser("verify", help="re-hash cached files against their manifests")
    p_verify.add_argument("--source-id")
    p_verify.add_argument("-v", "--verbose", action="store_true")
    p_verify.set_defaults(func=cmd_verify)

    p_where = sub.add_parser("where", help="print the cached path for a source")
    p_where.add_argument("source_id")
    p_where.set_defaults(func=cmd_where)

    p_sources = sub.add_parser("sources", help="list the declared registry")
    p_sources.add_argument("--channel", help="filter by channel, e.g. localization")
    p_sources.add_argument("--json", action="store_true")
    p_sources.set_defaults(func=cmd_sources)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (RegistryError, ManifestError, CacheError, PolicyError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
