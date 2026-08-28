#!/usr/bin/env python3
"""Turn an M-CSA catalytic-site entry into a site-context annotations document.

Curating sixteen functional-site cases by hand would mean transcribing residue
numbers and roles out of a web page, which is the same class of mistake that
produced a wrong expectation earlier in this collection. The annotations are
therefore generated from the M-CSA API and carry the release they came from.

Two vocabularies meet here and they are not the same:

* M-CSA describes what a residue *does* ("proton acceptor", "metal ligand",
  "nucleofuge") at fine grain -- two dozen distinct functions.
* site-context accepts five roles: nucleophile, acid_base, charge_relay,
  metal_ligand, unspecified.

The mapping below is lossy on purpose, and anything unrecognised becomes
``unspecified`` rather than being guessed into a role it may not have. It was
validated by regenerating M-CSA entry 1 and comparing against the hand-curated
annotations from the qualification v1 collection. The
panel's stratum is a separate judgement about the entry's *mechanism*, reported
alongside so a curator can select across strata without re-reading the entry.

Usage:
  python tools/harvest_mcsa_annotations.py --entry 1 --out annotations.json
  python tools/harvest_mcsa_annotations.py --survey 40      # classify candidates
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path
from typing import Any

API = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/"

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}

# M-CSA function -> site-context role, in priority order: the mechanistically
# decisive function wins, not whichever happens to be listed first. A residue
# annotated both "hydrogen bond acceptor" and "proton acceptor" is an acid/base
# residue; taking the first listed function would have called it something else.
#
# Interaction-level functions (hydrogen bonding, electrostatic stabilisation)
# deliberately map to "unspecified" rather than "charge_relay". A charge relay
# is a specific catalytic arrangement, and M-CSA does not annotate one; calling
# every polar contact a relay would assert more than the curation supports.
ROLE_PRIORITY = [
    ("nucleophile", {"nucleophile", "nucleofuge", "covalently attached"}),
    ("metal_ligand", {"metal ligand"}),
    # Only proton transfer by *this* residue counts as acid/base. "increase
    # acidity" and "increase basicity" describe modulating a neighbour's pKa,
    # which is a different claim: in M-CSA entry 1, Ser8 raises Asp7's basicity
    # while Asp7 is the residue that accepts the proton. Mapping modulation to
    # acid_base would have labelled Ser8 a catalytic base it is not.
    ("acid_base", {"proton acceptor", "proton donor", "proton relay"}),
]

# Which stratum an entry belongs to, most specific mechanism first. An entry
# with both a nucleophile and a metal is classified by the nucleophile, because
# that is the mechanistically decisive feature.
STRATUM_ORDER = [
    ("nucleophile_or_covalent", {"nucleophile", "nucleofuge", "covalently attached"}),
    ("metal_or_cofactor", {"metal ligand"}),
    ("acid_base", {"proton acceptor", "proton donor", "proton relay"}),
]


def fetch(url: str) -> Any:
    req = urllib.request.Request(url, headers={"User-Agent": "yauvi-qualification/2.0"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.loads(r.read().decode("utf-8"))


def classify(entry: dict) -> str:
    functions = {
        str(role.get("function", "")).lower()
        for residue in entry.get("residues") or []
        for role in residue.get("roles") or []
    }
    for stratum, markers in STRATUM_ORDER:
        if functions & markers:
            return stratum
    return "unclassified"


def reference_structure(entry: dict) -> tuple[str, str] | None:
    """The PDB id and chain the entry's reference residues are numbered against."""
    counts: dict[tuple[str, str], int] = {}
    for residue in entry.get("residues") or []:
        for chain in residue.get("residue_chains") or []:
            if chain.get("is_reference") and chain.get("pdb_id"):
                key = (str(chain["pdb_id"]).upper(), str(chain.get("chain_name") or "A"))
                counts[key] = counts.get(key, 0) + 1
    return max(counts, key=counts.get) if counts else None


def build_annotations(entry: dict) -> dict:
    ref = reference_structure(entry)
    if ref is None:
        raise SystemExit(f"M-CSA {entry.get('mcsa_id')} has no reference PDB chain")
    pdb_id, chain_name = ref
    reference_uniprot = entry.get("reference_uniprot_id")
    sites, skipped = [], []
    for residue in entry.get("residues") or []:
        chains = [c for c in (residue.get("residue_chains") or [])
                  if str(c.get("pdb_id", "")).upper() == pdb_id
                  and str(c.get("chain_name") or "A") == chain_name]
        if not chains:
            continue
        chain = chains[0]

        # site-context resolves an annotation position through the reference
        # sequence map, not through PDB author numbering, so positions must be
        # given in UniProt numbering. The two coincide in some entries -- 1B73
        # among them -- and diverge wherever a signal peptide or expression tag
        # shifts the deposited numbering. Using auth_resid made every curated
        # residue in 1NIA resolve to the wrong position and report a mismatch.
        sequences = [s for s in (residue.get("residue_sequences") or [])
                     if not reference_uniprot or s.get("uniprot_id") == reference_uniprot]
        if not sequences:
            skipped.append({"auth_resid": chain.get("auth_resid"), "code": chain.get("code"),
                            "reason": "no reference-sequence position"})
            continue
        sequence = sequences[0]
        position = sequence.get("resid")

        one = THREE_TO_ONE.get(str(sequence.get("code") or chain.get("code", "")).capitalize())
        if one is None:
            # A modified or non-standard curated residue. Recorded as skipped
            # rather than mapped to a standard letter it is not.
            skipped.append({"auth_resid": chain.get("auth_resid"), "code": chain.get("code")})
            continue
        functions = {str(r.get("function", "")).lower() for r in residue.get("roles") or []}
        role = "unspecified"
        for candidate, markers in ROLE_PRIORITY:  # decisive mechanism wins
            if functions & markers:
                role = candidate
                break
        sites.append({
            "position": position,
            "pdb_auth_resid": chain.get("auth_resid"),
            "role": role,
            "type": "metal_ligand" if role == "metal_ligand" else "active_site",
            "expected_residues": [one],
            "detail": residue.get("roles_summary") or "",
        })
    sites.sort(key=lambda s: (s["position"] is None, s["position"]))
    return {
        "declared_cofactors": [],
        "sites": sites,
        "provenance": {
            "source": "M-CSA (Mechanism and Catalytic Site Atlas)",
            "entry": f"M-CSA:{entry.get('mcsa_id')}",
            "url": f"{API}{entry.get('mcsa_id')}/",
            "reference_pdb": pdb_id, "reference_chain": chain_name,
            "reference_uniprot": entry.get("reference_uniprot_id"),
            "enzyme": entry.get("enzyme_name"), "ec": entry.get("all_ecs"),
            "stratum": classify(entry),
            "skipped_nonstandard_residues": skipped,
            "role_mapping_note": ("M-CSA functions are mapped onto the five roles site-context "
                                  "accepts; unrecognised functions become 'unspecified' rather "
                                  "than being guessed into a role."),
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--entry", type=int, help="M-CSA entry id to convert.")
    ap.add_argument("--out", type=Path, help="Where to write the annotations document.")
    ap.add_argument("--survey", type=int, metavar="N",
                    help="Classify the first N entries and print candidates by stratum.")
    args = ap.parse_args(argv)

    if args.survey:
        page = fetch(f"{API}?format=json&page_size={args.survey}")
        for e in page.get("results", []):
            ref = reference_structure(e)
            if not ref:
                continue
            n = len([r for r in e.get("residues") or [] if r.get("residue_chains")])
            print(f"{classify(e):24s} M-CSA:{e['mcsa_id']:<5} {ref[0]}:{ref[1]} "
                  f"unp={e.get('reference_uniprot_id')} residues={n} {str(e.get('enzyme_name'))[:36]}")
        return 0

    if args.entry is None or args.out is None:
        ap.error("--entry and --out are required unless --survey is used")
    doc = build_annotations(fetch(f"{API}{args.entry}/?format=json"))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    p = doc["provenance"]
    print(f"M-CSA:{args.entry} -> {args.out}  stratum={p['stratum']} "
          f"pdb={p['reference_pdb']}:{p['reference_chain']} sites={len(doc['sites'])}"
          + (f" skipped={len(p['skipped_nonstandard_residues'])}" if p["skipped_nonstandard_residues"] else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
