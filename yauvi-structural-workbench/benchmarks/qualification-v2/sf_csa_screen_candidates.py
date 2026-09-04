#!/usr/bin/env python3
"""Find PDB entries that can actually serve as sf-csa queries in a given SCOP group.

Selecting entries by recall does not work: four of the first twelve chosen this
way turned out to have unresolved residues, and sf-csa refuses any query whose
PDB sequence does not exactly match its FASTA. This enumerates instead.

For each candidate it reports the three things that decide usability:

  * the SCOP fold and superfamily, from the PDBe SIFTS mapping
  * whether the entry is a single unambiguous SCOP domain
  * whether the observed (ATOM) sequence equals the deposited FASTA exactly

Only entries passing all three can be frozen without either substituting a
truncated FASTA or relaxing the module's own input contract.

Curation only. Panel execution forbids network access; the runner never calls this.

    python3 sf_csa_screen_candidates.py "RmlC-like cupins" --limit 12
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
PDBE = "https://www.ebi.ac.uk/pdbe/api/mappings/scop/"
FASTA = "https://www.rcsb.org/fasta/entry/"
PDBFILE = "https://files.rcsb.org/download/"


def _get(url: str, timeout: int = 30) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        return fh.read()


def entries_in(lineage_name: str, rows: int) -> list[str]:
    """PDB entries RCSB annotates with this SCOP lineage name."""
    query = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_instance_annotation.annotation_lineage.name",
            "operator": "exact_match", "value": lineage_name}},
        "return_type": "polymer_instance",
        "request_options": {"paginate": {"start": 0, "rows": rows},
                            "results_content_type": ["experimental"]},
    }
    url = SEARCH + "?json=" + urllib.parse.quote(json.dumps(query))
    # The search endpoint drops large result pages intermittently; retry rather
    # than let a transient disconnect read as "no candidates in this group".
    import time
    for attempt in range(3):
        try:
            data = json.loads(_get(url))
            break
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))
    seen: list[str] = []
    for hit in data.get("result_set", []):
        entry = hit["identifier"].split(".")[0].lower()
        if entry not in seen:
            seen.append(entry)
    return seen


def scop_of(entry: str) -> list[dict]:
    try:
        data = json.loads(_get(PDBE + entry))
    except Exception:
        return []
    return [{"sccs": i.get("sccs"),
             "fold": (i.get("fold") or {}).get("description"),
             "superfamily": (i.get("superfamily") or {}).get("description")}
            for _, b in data.items() for _, i in (b.get("SCOP") or {}).items()]


def _chain_sequence(fasta_text: str, chain: str) -> str | None:
    """The deposited sequence for one chain.

    An RCSB entry FASTA holds one record per distinct entity, headed
    `>1ABC_1|Chains A, B|...`. Concatenating every record -- the obvious way to
    read the file, and how this was first written -- silently sums unrelated
    chains and makes a single-chain comparison look truncated. Select the record
    naming the chain instead.
    """
    records: list[tuple[str, list[str]]] = []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            records.append((line, []))
        elif line.strip() and records:
            records[-1][1].append(line.strip())
    for header, body in records:
        fields = header.split("|")
        if len(fields) < 2:
            continue
        label = fields[1]
        if not label.lower().startswith("chain"):
            continue
        # Labels look like "Chain A", "Chains A, B" or "Chains D[auth A], E[auth B]".
        # The auth identifier is the one the coordinate file uses, so it is the one
        # that must match the chain read out of the ATOM records.
        names: list[str] = []
        for raw in label.split(None, 1)[-1].split(","):
            raw = raw.strip()
            if "[auth" in raw:
                names.append(raw.split("[auth", 1)[1].rstrip("] ").strip())
            elif raw.endswith("]"):
                names.append(raw.rstrip("] ").strip())
            else:
                names.append(raw)
        if chain in names:
            return "".join(body)
    return "".join(records[0][1]) if len(records) == 1 else None


def sequence_exact(entry: str, chain: str = "A", retries: int = 2) -> tuple[str, int, int]:
    """Compare the observed structure sequence with the deposited one for `chain`.

    Returns a status, never a bare boolean. An earlier version collapsed every
    failure into False, so a dropped connection was indistinguishable from a real
    length mismatch -- and four entries were briefly recorded as truncated when
    the fetch had simply failed. A result that cannot be obtained is reported as
    such, not as a negative finding.

    status is one of:
      "exact"        observed sequence equals the deposited chain sequence
      "truncated"    they differ (the usual cause is unresolved residues)
      "no_chain"     the FASTA has no record naming this chain
      "unreadable"   the coordinates could not be parsed for this chain
      "fetch_failed" the files could not be retrieved after retries
    """
    from pathlib import Path
    import tempfile
    import time
    from sf_csa.core import pdb_sequence
    up = entry.upper()
    pdb_bytes = fasta_text = None
    for attempt in range(retries + 1):
        try:
            pdb_bytes = _get(PDBFILE + up + ".pdb")
            fasta_text = _get(FASTA + up + "/download").decode()
            break
        except Exception:
            if attempt == retries:
                return ("fetch_failed", 0, 0)
            time.sleep(1.5 * (attempt + 1))
    seq = _chain_sequence(fasta_text, chain)
    if seq is None:
        return ("no_chain", 0, 0)
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        tmp.write(pdb_bytes)
        path = Path(tmp.name)
    try:
        pseq, _, _ = pdb_sequence(path, chain)
    except Exception:
        return ("unreadable", 0, len(seq))
    finally:
        path.unlink(missing_ok=True)
    return ("exact" if pseq == seq else "truncated", len(pseq), len(seq))


ALPHAFOLD = "https://alphafold.ebi.ac.uk/api/prediction/"
UNIPROT_MAP = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/"


def alphafold_concordance(entry: str, deposited_length: int) -> tuple[str, str, int]:
    """Can this protein be a query as an AlphaFold model, and is it the same entity?

    sf-csa's contract asks for exact predicted monomers, and a predicted model
    trivially satisfies the sequence-exact requirement because every residue has
    coordinates. But the model covers the whole UniProt protein, while SCOP
    classified the crystallised construct. Where those differ the model is a
    different molecule from the one whose fold was classified -- 1SNC is a
    149-residue construct against a 231-residue model, and 1URN is one RRM domain
    against a 282-residue multi-domain protein.

    Concordance is therefore required, not assumed: the model may differ from the
    deposited chain by at most one residue, which absorbs an initiator methionine
    and nothing larger.

    Returns (status, uniprot_accession, model_length).
    """
    try:
        mapping = json.loads(_get(UNIPROT_MAP + entry))
        accessions = list((mapping.get(entry) or {}).get("UniProt") or {})
        if not accessions:
            return ("no_uniprot", "", 0)
        acc = accessions[0]
        meta = json.loads(_get(ALPHAFOLD + acc))
        if not meta:
            return ("no_model", acc, 0)
        model = _get(meta[0]["pdbUrl"])
    except Exception:
        return ("fetch_failed", "", 0)

    from pathlib import Path as _P
    import tempfile
    from sf_csa.core import pdb_sequence
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        tmp.write(model)
        path = _P(tmp.name)
    try:
        pseq, _, _ = pdb_sequence(path, "A")
    except Exception:
        return ("unreadable_model", acc, 0)
    finally:
        path.unlink(missing_ok=True)
    delta = len(pseq) - deposited_length
    return ("concordant" if abs(delta) <= 1 else f"differs{delta:+d}", acc, len(pseq))


UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/"
UNIPROT_PROTEOMES = "https://rest.uniprot.org/proteomes/search"


def reference_proteome(accession: str) -> tuple[str, int]:
    """The reference proteome this protein is actually **in**, if any.

    Corrected 2026-09-01. This function used to search proteomes by the
    organism's taxonomy id and return the largest reference proteome for that
    taxon, which answers a different question: whether the organism has been
    sequenced, not whether this protein is in the result. Five of the twelve
    queries frozen on the old answer -- P00193, P23370, P45850, P00138, P00147 --
    are in no reference proteome at all, and their declared proteome files load
    fine and contain zero occurrences of them. `run_pipeline` raises only when the
    file is missing, so the sequence leg silently vanished for thirteen of sixteen
    records, including two exact self-matches.

    Membership is read from the entry's own Proteomes cross-reference, which is
    the authoritative statement that UniProt places this accession in that
    proteome. An entry listing none returns ("", 0) and is unusable, however well
    sequenced its organism happens to be.

    Returns (upid or "", protein_count).
    """
    try:
        entry = json.loads(_get(UNIPROT_ENTRY + accession + ".json?fields=xref_proteomes"))
    except Exception:
        return ("", 0)
    listed = [x["id"] for x in entry.get("proteomes", [])] or [
        r["id"] for r in entry.get("uniProtKBCrossReferences", [])
        if r.get("database") == "Proteomes"
    ]
    if not listed:
        return ("", 0)
    # Only reference proteomes can actually be downloaded; non-reference ones are
    # unfiltered and their sequences cannot be retrieved. Prefer a reference
    # proteome the entry belongs to, and report its size so the caller can see
    # what the sequence leg will be searched against.
    for upid in listed:
        try:
            record = json.loads(_get(f"https://rest.uniprot.org/proteomes/{upid}.json"))
        except Exception:
            continue
        if record.get("proteomeType") == "Reference proteome":
            return (upid, record.get("proteinCount") or 0)
    return ("", 0)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("lineage", help='SCOP lineage name, e.g. "OB-fold"')
    ap.add_argument("--limit", type=int, default=10, help="candidates to screen")
    ap.add_argument("--rows", type=int, default=60, help="entries to pull before screening")
    ap.add_argument("--exclude-superfamily", default=None,
                    help="skip candidates in this superfamily (to find an analogy partner)")
    args = ap.parse_args()

    try:
        candidates = entries_in(args.lineage, args.rows)
    except Exception as exc:
        print(f"search failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    print(f"{len(candidates)} entries annotated {args.lineage!r}; screening up to {args.limit}\n")
    print(f"{'entry':7} {'sccs':12} {'superfamily':30} {'uniprot':9} {'alphafold':>13} proteome")
    print("-" * 104)

    usable, screened = [], 0
    for entry in candidates:
        if screened >= args.limit:
            break
        domains = scop_of(entry)
        if len(domains) != 1:
            continue                      # ambiguous or unclassified; not worth the fetch
        d = domains[0]
        if args.exclude_superfamily and d["superfamily"] == args.exclude_superfamily:
            continue
        screened += 1
        status, obs, full = sequence_exact(entry)
        af, acc, af_len = alphafold_concordance(entry, full)
        upid, n_prot = reference_proteome(acc) if acc else ("", 0)
        ok = af == "concordant" and bool(upid)
        prot = f"{upid}({n_prot})" if upid else "NO PROTEOME"
        print(f"{entry:7} {str(d['sccs']):12} {str(d['superfamily'])[:30]:30} "
              f"{acc:9} {af:>13} {prot}")
        if ok:
            usable.append((entry, d["sccs"], d["superfamily"], acc, upid))

    print(f"\nusable ({len(usable)}):")
    for e, sccs, sf, acc, upid in usable:
        print(f"  {e}  {sccs:12} {acc:9} {upid:12} {sf}")
    return 0 if usable else 1


if __name__ == "__main__":
    raise SystemExit(main())
