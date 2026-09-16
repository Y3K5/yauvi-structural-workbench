#!/usr/bin/env python3
"""Assign sf-csa relationship strata from SCOP, not from the curator's memory.

The fold_analogy stratum is the one that decides whether the panel's
false-positive gate means anything, and it is the easiest to get wrong in the
flattering direction: a pair curated as "analogy" that is really a distant
homolog turns a passing gate into a vacuous one. Two proposed families were
withdrawn for exactly this reason after the literature was checked.

So the stratum is derived here from a classification a reviewer can re-query,
using the definition in Russell, Saqi, Sayle, Bates and Sternberg (J Mol Biol
1997): analogues share a SCOP fold but sit in different superfamilies, with
little evidence of common ancestry; homologues share a superfamily.

    same superfamily              -> homologous_superfamily
    same fold, other superfamily  -> fold_analogy
    different fold                -> unrelated

`exact` is not decided here: it is accession identity, settled before any
structural reasoning.

Network is used for curation only. Panel *execution* forbids it; this script is
not called from the runner.

    python3 sf_csa_stratum_from_scop.py 1tim 8tim 256b 1fha ...
"""
from __future__ import annotations

import json
import sys
import urllib.request
from itertools import combinations

PDBE = "https://www.ebi.ac.uk/pdbe/api/mappings/scop/"


def scop(entry: str) -> list[dict]:
    """Every SCOP domain PDBe records for one entry."""
    try:
        with urllib.request.urlopen(PDBE + entry.lower(), timeout=30) as fh:
            data = json.load(fh)
    except Exception as exc:
        print(f"  {entry}: lookup failed ({type(exc).__name__})", file=sys.stderr)
        return []
    out = []
    for _, body in data.items():
        for sunid, info in (body.get("SCOP") or {}).items():
            out.append({
                "sunid": sunid,
                "sccs": info.get("sccs"),
                "fold": (info.get("fold") or {}).get("description"),
                "fold_id": (info.get("fold") or {}).get("sunid"),
                "superfamily": (info.get("superfamily") or {}).get("description"),
                "superfamily_id": (info.get("superfamily") or {}).get("sunid"),
            })
    return out


def primary(domains: list[dict]) -> dict | None:
    """The entry's largest classified domain, by first appearance.

    An entry with several domains cannot be reduced to one fold without a
    choice; the choice is recorded rather than hidden, and any entry whose
    domains disagree on fold is reported so it can be excluded or split.
    """
    if not domains:
        return None
    folds = {d["fold_id"] for d in domains}
    d = dict(domains[0])
    d["multi_fold"] = len(folds) > 1
    d["domain_count"] = len(domains)
    return d


def stratum(a: dict, b: dict) -> str:
    if a["superfamily_id"] == b["superfamily_id"]:
        return "homologous_superfamily"
    if a["fold_id"] == b["fold_id"]:
        return "fold_analogy"
    return "unrelated"


def main(entries: list[str]) -> int:
    if not entries:
        print(__doc__)
        return 2
    resolved: dict[str, dict] = {}
    print(f"{'entry':7} {'sccs':12} {'fold':38} superfamily")
    print("-" * 110)
    for e in entries:
        p = primary(scop(e))
        if p is None:
            print(f"{e:7} {'-':12} {'NOT CLASSIFIED IN SCOP':38} -")
            continue
        resolved[e.lower()] = p
        flag = "  [multi-fold entry]" if p["multi_fold"] else ""
        print(f"{e:7} {str(p['sccs']):12} {str(p['fold'])[:38]:38} {str(p['superfamily'])[:44]}{flag}")

    if len(resolved) < 2:
        print("\nnot enough classified entries to form a pair")
        return 1

    print(f"\n{'pair':17} {'stratum':24} basis")
    print("-" * 110)
    buckets: dict[str, list[str]] = {}
    for (ea, a), (eb, b) in combinations(resolved.items(), 2):
        s = stratum(a, b)
        if s == "homologous_superfamily":
            basis = f"shared superfamily: {a['superfamily']}"
        elif s == "fold_analogy":
            basis = f"shared fold {a['fold']!r}; superfamilies differ: {a['superfamily']} vs {b['superfamily']}"
        else:
            basis = f"different folds: {a['fold']} vs {b['fold']}"
        buckets.setdefault(s, []).append(f"{ea}/{eb}")
        print(f"{ea + '/' + eb:17} {s:24} {basis[:66]}")

    print("\ncounts by stratum (panel needs 4 of each):")
    for s in ("homologous_superfamily", "fold_analogy", "unrelated"):
        got = buckets.get(s, [])
        print(f"  {s:24} {len(got):>4}   {'OK' if len(got) >= 4 else 'SHORT'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
