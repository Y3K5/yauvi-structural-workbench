#!/usr/bin/env python3
"""Auditable structure/function comparison across species.

The module deliberately keeps structure search, sequence homology, curated function and
human target decisions in separate columns. A Foldseek title can nominate an architecture;
it can never, by itself, become a function claim.

Every campaign-specific judgement is a manifest entry, not a literal here: the
mechanism-family patterns, the groups whose structural context is contested, the
sets of groups that share a framework but not a substrate, and the title traps
are all supplied by the database manifest. The defaults below are the
periodontal-pathogen values, retained so an existing manifest that omits them
behaves exactly as before.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Mapping


AA3 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M",
}
TSV_FIELDS = [
    "query_id","target_id","database","target_description","evalue","fident",
    "aligned_length","query_length","target_length","query_coverage","target_coverage",
    "query_tm_score","target_tm_score","aligned_tm_score","rmsd","structural_category",
    "function_classification","classification_basis","evidence_class","limitation",
]


class SFCSError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def database_bundle_files(path: Path) -> list[Path]:
    """Return every file that constitutes a frozen structure-search database.

    Foldseek databases are prefix bundles, not single files.  Hashing only the
    main prefix leaves indexes and encoded coordinate/sequence sidecars free to
    drift.  A plain stand-in file remains a one-file bundle for fixtures and
    backward-compatible small adapters.
    """
    if path.is_dir():
        return sorted(item for item in path.rglob("*") if item.is_file())
    matches = sorted(
        item for item in path.parent.glob(path.name + "*")
        if item.is_file()
    )
    return matches or ([path] if path.is_file() else [])


def database_bundle_file_checksums(path: Path) -> dict[str, str]:
    files = database_bundle_files(path)
    if not files:
        raise SFCSError(f"structure database bundle absent: {path}")
    base = path if path.is_dir() else path.parent
    return {
        item.relative_to(base).as_posix(): sha256(item)
        for item in files
    }


def database_bundle_checksum(path: Path) -> str:
    checksums = database_bundle_file_checksums(path)
    if len(checksums) == 1 and path.is_file() and path.name in checksums:
        return checksums[path.name]
    payload = "".join(f"{name}\t{digest}\n" for name, digest in sorted(checksums.items()))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def sequence_sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text())


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=False) + "\n")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_fasta(path: Path) -> list[dict]:
    records, header, seq = [], None, []
    with path.open(errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append({"header": header, "id": fasta_id(header), "sequence": "".join(seq).upper()})
                header, seq = line[1:], []
            elif line:
                seq.append(re.sub(r"[^A-Za-z]", "", line))
    if header is not None:
        records.append({"header": header, "id": fasta_id(header), "sequence": "".join(seq).upper()})
    return records


def fasta_id(header: str) -> str:
    token = header.split()[0]
    bits = token.split("|")
    if len(bits) >= 3 and bits[0] in {"sp", "tr", "ref", "gb"}:
        return bits[1]
    return token


def norm_id(value: str) -> str:
    value = value.rsplit("|", 1)[-1] if value.startswith("p") and "|" in value else value
    bits = value.split("|")
    if len(bits) >= 3 and bits[0] in {"sp", "tr", "ref", "gb"}:
        value = bits[1]
    return value.replace(".pdb", "").replace("_oriented", "")


def select_fasta(path: Path, accession: str) -> dict:
    hits = [r for r in parse_fasta(path) if norm_id(r["id"]) == accession]
    if len(hits) != 1:
        raise SFCSError(f"expected exactly one {accession} record in {path}; found {len(hits)}")
    return hits[0]


def pdb_sequence(path: Path, chain: str = "A") -> tuple[str, list[tuple[int, str]], list[tuple[float,float,float]]]:
    residues, seen, coords = [], set(), []
    with path.open(errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
                continue
            if line[21].strip() not in {chain, ""}:
                continue
            resn = line[17:20].strip().upper()
            if resn not in AA3:
                continue
            key = (line[21].strip(), line[22:27])
            if key not in seen:
                seen.add(key)
                try:
                    resi = int(line[22:26])
                except ValueError:
                    resi = len(residues) + 1
                residues.append((resi, AA3[resn]))
            if line[12:16].strip() == "CA":
                try:
                    coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    pass
    return "".join(x[1] for x in residues), residues, coords


def geometry(coords: list[tuple[float,float,float]]) -> dict:
    if not coords:
        return {"status": "unavailable_no_ca_coordinates"}
    n = len(coords)
    center = [sum(p[i] for p in coords) / n for i in range(3)]
    rg = math.sqrt(sum(sum((p[i]-center[i])**2 for i in range(3)) for p in coords) / n)
    extents = [max(p[i] for p in coords) - min(p[i] for p in coords) for i in range(3)]
    return {
        "status": "computed_coordinate_summary",
        "ca_atoms": n,
        "radius_of_gyration_angstrom": round(rg, 3),
        "axis_aligned_extents_angstrom": [round(x, 3) for x in extents],
        "limit": "coordinate summary is not a lumen, pocket, electrostatic or active-state calculation",
    }


def transform_pdb(source: Path, destination: Path, rotation: str, translation: str) -> None:
    """Apply Foldseek's target-to-query rigid transform to all PDB atom coordinates."""
    u=[float(x) for x in rotation.split(",") if x]
    t=[float(x) for x in translation.split(",") if x]
    if len(u)!=9 or len(t)!=3: raise SFCSError(f"invalid Foldseek transform for {source.name}")
    lines=[]
    for line in source.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM  ","HETATM")) and len(line)>=54:
            try:
                p=[float(line[30:38]),float(line[38:46]),float(line[46:54])]
                v=[sum(u[i*3+j]*p[j] for j in range(3))+t[i] for i in range(3)]
                line=f"{line[:30]}{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}{line[54:]}"
            except ValueError:
                pass
        lines.append(line)
    destination.write_text("\n".join(lines)+"\n")


def write_superposition_html(path: Path, query_pdb: Path, target_pdb: Path, query_id: str,
                             target_id: str, tm_score: str, rmsd: str) -> None:
    """Write the superposition viewer, and say so when it cannot draw.

    The page loads 3Dmol from `assets/vendor/3Dmol-min.js`, relative to the
    release. **No release has ever contained that file** -- nothing in this
    module or in the tooling writes it -- so until one does, every viewer this
    function has produced shows its caption over an empty panel and fails
    silently in the console. That reads as a viewer that broke, which is worse
    than one that was never possible: a reader cannot tell whether the
    superposition failed or the page did.

    So the page checks for the library and, when it is absent, replaces the
    viewer with a plain statement of what is missing and where the coordinates
    are. The structures stay embedded either way: they are the substance, and
    they are readable in a text editor without any viewer at all.

    Keeping the script tag is deliberate. Drop 3Dmol-min.js into
    `assets/vendor/` beside a release and every page in it starts rendering,
    with no regeneration.
    """
    qp=query_pdb.read_text().replace("`","\\`").replace("${","\\${")
    tp=target_pdb.read_text().replace("`","\\`").replace("${","\\${")
    html=f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<style>html,body,#v{{margin:0;width:100%;height:100%;background:#07111f}}.note{{position:absolute;z-index:2;left:12px;top:10px;background:#07111fdd;color:#dbeafe;padding:8px 10px;border-radius:6px;font:12px system-ui}}b{{color:#67e8f9}}i{{color:#f0abfc}}.absent{{position:absolute;z-index:1;inset:0;display:flex;align-items:center;justify-content:center;padding:24px}}.absent p{{max-width:44em;color:#94a3b8;font:13px/1.6 system-ui;text-align:left}}code{{color:#cbd5e1}}</style>
</head><body><div class="note"><b>{query_id}</b> query · <i>{target_id}</i> aligned target · TM {tm_score} · RMSD {rmsd} Å<br>Foldseek rigid superposition; predicted monomers, not an active assembly.</div><div id="v"></div><script src="../../assets/vendor/3Dmol-min.js"></script><script>
if (typeof $3Dmol === "undefined") {{
  document.getElementById("v").innerHTML =
    '<div class="absent"><p><b>No 3D view here.</b> This page expects 3Dmol.js at ' +
    '<code>assets/vendor/3Dmol-min.js</code>, relative to the release, and this release ' +
    'does not ship it. Nothing is wrong with the superposition: the TM-score and RMSD ' +
    'above are the measured result, and both structures are embedded in this file as PDB ' +
    'text, readable in any editor. To render it, put a copy of 3Dmol.js at that path and ' +
    'reload.</p></div>';
}} else {{
  const v=$3Dmol.createViewer(document.getElementById("v"),{{backgroundColor:"#07111f"}});
  v.addModel(`{qp}`,"pdb");v.setStyle({{model:0}},{{cartoon:{{color:"#22d3ee",opacity:.78}}}});
  v.addModel(`{tp}`,"pdb");v.setStyle({{model:1}},{{cartoon:{{color:"#e879f9",opacity:.62}}}});
  v.zoomTo();v.render();
}}
</script></body></html>'''
    path.write_text(html)


def resolve(manifest_path: Path, rel: str, path_base: str) -> Path:
    return (manifest_path.parent / path_base / rel).resolve()


def run_cmd(args: list[str], cwd: Path | None = None) -> None:
    proc = subprocess.run(args, cwd=cwd, text=True, capture_output=True)
    if proc.returncode:
        detail = (proc.stderr or proc.stdout or "no diagnostic output").strip()
        raise SFCSError(f"command failed ({proc.returncode}): {' '.join(args[:3])}\n{detail[-4000:]}")


def tool_version(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise SFCSError(f"required executable not found: {name}")
    if name == "foldseek":
        out = subprocess.run([exe, "version"], text=True, capture_output=True, check=True).stdout.strip()
    elif name == "diamond":
        out = subprocess.run([exe, "version"], text=True, capture_output=True, check=True).stdout.strip()
    else:
        out = exe
    return out


# The mechanism families a PDB title can be assigned to, in order. These are the
# periodontal-pathogen defaults; a campaign targeting other organisms overrides them
# with `mechanism_families` in its database manifest rather than editing this file.
# `refine` re-tests a matched family for a narrower one, which is how a TonB-dependent
# transporter becomes a heme receptor.
# A PDB title that must never carry a hit up to a function claim. "Toluene
# transporter" is the recorded case: same barrel fold, unrelated substrate.
DEFAULT_TITLE_TRAPS = [
    {"substring": "toluene",
     "must_not_promote_to": ["exact_function_supported", "probable_same_function"]},
]

# The closed interpretation vocabulary. Six values, ordered from strongest to
# weakest claim. `classify_hit` may return nothing outside this set, and a release
# whose manifest names a different vocabulary is not comparable with one that does
# not — which is why it is recorded in the manifest rather than assumed.
CLASSIFICATION_VOCABULARY = [
    "exact_function_supported",
    "probable_same_function",
    "same_mechanism_class",
    "structural_analogy_only",
    "candidate_functional_divergence",
    "unresolved_or_conflicted",
]

# Groups whose structural evidence is known to be contested, so a same-group match
# may not be promoted to a shared-function claim. Campaign biology, not a property
# of the algorithm — which is why it is overridable.
DEFAULT_CONTESTED_GROUPS = [
    {"group": "msp_contested",
     "reason": "Msp structural contexts remain contested",
     "evidence_class": "E2_direct_biology_plus_E4_structure",
     "boundary": "no single oligomeric or active pose is asserted"},
]

# Sets of mechanism groups that share an architectural framework but not a
# substrate. A cross-group match inside one of these sets is functional
# divergence, not analogy.
DEFAULT_DIVERGENCE_SETS = [
    {"groups": ["susc_raga_importer", "tonb_heme_receptor", "generic_om_barrel"],
     "reason": "shared outer-membrane transporter framework with differing functional specificity",
     "evidence_class": "E4_computational",
     "boundary": "PDB title or fold cannot establish substrate"},
]

DEFAULT_MECHANISM_FAMILIES = [
    {"group": "omp85_bama", "pattern": r"\bbama\b|omp85|outer membrane protein assembly"},
    {"group": "susc_raga_importer", "pattern": r"\braga\b|\bsusc\b|tonb-dependent|tonb linked",
     "refine": [{"group": "tonb_heme_receptor", "pattern": r"heme|haem|hemoglobin"}]},
    {"group": "t9ss_porg", "pattern": r"\bporg\b|type ix secretion"},
    {"group": "msp_contested", "pattern": r"major surface protein|\bmsp\b"},
    {"group": "generic_om_barrel", "pattern": r"porin|outer[- ]membrane|barrel"},
]


def classify_title(title: str, families: list | None = None) -> str:
    low = title.lower()
    for family in (families if families is not None else DEFAULT_MECHANISM_FAMILIES):
        if not re.search(family["pattern"], low):
            continue
        for narrower in family.get("refine", ()):
            if re.search(narrower["pattern"], low):
                return narrower["group"]
        return family["group"]
    return "unknown"


def structural_category(hit: dict, whole_cov: float, same_tm: float) -> str:
    tm = float(hit.get("alntmscore") or 0)
    qcov = float(hit.get("qcov") or 0)
    tcov = float(hit.get("tcov") or 0)
    if tm >= same_tm and qcov >= whole_cov and tcov >= whole_cov:
        return "whole_architecture_match"
    if tm >= same_tm:
        return "domain_or_partial_match"
    return "below_structural_similarity_threshold"


RESERVED_COMPUTED_FIELDS = (
    "rbh",
    "orthology_status",
    "identity_status",
    "sequence_length",
    "structure_residue_count",
    "geometry",
    "uniprot_annotation",
)


def reject_reserved_fields(record: Mapping, label: str) -> None:
    """Refuse a curator record that carries a field the pipeline computes.

    Curator keys survive verbatim into `target_meta` (`row = {**q, ...}`), so a
    manifest field named like a computed one is indistinguishable from evidence
    downstream. `rbh` is the load-bearing case: asserting it would satisfy or
    violate the panel's false-positive gate without any reciprocal best hit
    having been computed. Rejecting at read time keeps assertion and measurement
    separable, which is the property the gate depends on.
    """
    present = sorted(f for f in RESERVED_COMPUTED_FIELDS if f in record)
    if present:
        raise SFCSError(
            f"{label} declares reserved computed field(s): {', '.join(present)}. "
            "These are produced by the pipeline and may not be supplied by a manifest."
        )


def classify_hit(query: dict, hit: dict, category: str, target_meta: dict | None = None,
                 families: list | None = None, contested: list | None = None,
                 divergence_sets: list | None = None, *, rbh: bool = False) -> tuple[str,str,str,str]:
    qgroup = query["mechanism_group"]
    tgroup = (target_meta or {}).get("mechanism_group") or classify_title(hit.get("theader", ""), families)
    target_id = norm_id(hit.get("target", ""))
    contested = DEFAULT_CONTESTED_GROUPS if contested is None else contested
    divergence_sets = DEFAULT_DIVERGENCE_SETS if divergence_sets is None else divergence_sets

    if target_id == query["accession"]:
        return ("exact_function_supported", "exact accession self-match", "E1_exact_identity", "self-match is a control, not an independent function experiment")
    if category == "below_structural_similarity_threshold":
        return ("unresolved_or_conflicted", "structural threshold not met", "E4_computational", "no same-fold interpretation allowed")
    if qgroup == tgroup and qgroup != "unknown":
        # A contested group may not be promoted, however good the structural match.
        entry = next((c for c in contested if c["group"] == qgroup), None)
        if entry is not None:
            return ("unresolved_or_conflicted", entry["reason"],
                    entry.get("evidence_class", "E2_direct_biology_plus_E4_structure"),
                    entry.get("boundary", "no single active pose is asserted"))
        if rbh and category == "whole_architecture_match":
            return ("probable_same_function", "RBH plus compatible whole architecture and mechanism annotation", "E3_curated_plus_E4_computational", "substrate and native activity still require direct validation")
        return ("same_mechanism_class", "compatible structural architecture and independently named mechanism class", "E3_or_E4_bounded_inference", "substrate, activity, virulence and vaccine suitability are not transferred")
    for entry in divergence_sets:
        if {qgroup, tgroup} <= set(entry["groups"]):
            return ("candidate_functional_divergence", entry["reason"],
                    entry.get("evidence_class", "E4_computational"),
                    entry.get("boundary", "fold cannot establish substrate"))
    if category == "domain_or_partial_match":
        return ("structural_analogy_only", "partial structural similarity without concordant independent function evidence", "E4_computational", "partial fold similarity is not functional equivalence")
    return ("structural_analogy_only", "whole-architecture similarity without concordant independent function evidence", "E4_computational", "fold similarity alone is non-functional evidence")


def foldseek_search(query_pdb: Path, target: Path, out: Path, tmp: Path, max_hits: int, evalue: str,
                    exhaustive: bool = False) -> list[dict]:
    """Search one query against a structure database.

    `exhaustive` exists because two filters act in series, and only the first is
    invisible. Foldseek prefilters before scoring, and the prefilter ignores the
    e-value entirely; the e-value then applies to whatever survives. Measured on
    this campaign, for P00198 against the twelve-entry campaign database:

        -e 0.01                        2 rows
        -e 10000                       2 rows      prefilter, e-value irrelevant
        -e 0.01   --exhaustive-search  2 rows      no prefilter, e-value binds
        -e 10000  --exhaustive-search  12 rows

    So neither setting alone reports a distant pair. That matters to a benchmark,
    because "no row" and "a row saying below threshold" are different facts about
    a comparison, and a panel that cannot tell them apart cannot say whether a
    pair was rejected or never examined. A campaign that needs every declared
    pair reported must set both, and let its own declared thresholds --
    same_fold_tm and whole_architecture_coverage -- do the scientific filtering
    rather than leaving it to a search heuristic.

    It is off by default: exhaustive search is quadratic, and a release scanning
    a large database should keep the prefilter. A qualification campaign over a
    twelve-entry database declares it on.
    """
    fields = "query,target,evalue,fident,alnlen,qlen,tlen,qcov,tcov,qtmscore,ttmscore,alntmscore,rmsd,theader,u,t"
    cmd = ["foldseek","easy-search",str(query_pdb),str(target),str(out),str(tmp),
           "--format-output",fields,"-e",evalue,"--max-seqs",str(max_hits)]
    if exhaustive:
        cmd += ["--exhaustive-search","1"]
    run_cmd(cmd)
    names = fields.split(",")
    rows = []
    if out.exists():
        with out.open(errors="replace") as fh:
            for line in fh:
                vals = line.rstrip("\n").split("\t")
                if len(vals) >= len(names):
                    rows.append(dict(zip(names, vals)))
    return rows


def discover_campaign_structures(db: dict, db_path: Path, root: Path, stage: Path) -> tuple[dict, dict]:
    aliases: dict[str,str] = {}
    for rel in db.get("seqmatch_tables", []):
        p = resolve(db_path, rel, db["path_base"])
        if not p.exists():
            continue
        for r in csv.DictReader(p.open(), delimiter="\t"):
            if r.get("resolved_accession"):
                aliases[r["refseq_id"]] = r["resolved_accession"]
                aliases[r["resolved_accession"]] = r["refseq_id"]
    candidates = []
    for rel in db["campaign_structure_roots"]:
        p = resolve(db_path, rel, db["path_base"])
        if p.exists():
            candidates.extend(p.rglob("*.pdb") if p.is_dir() else [p])
    stage.mkdir(parents=True, exist_ok=True)
    meta, seen = {}, set()
    for src in sorted(candidates):
        if "_afdb_cache" in src.parts:
            continue
        acc = norm_id(src.stem)
        key = aliases.get(acc, acc)
        if key in seen:
            continue
        seen.add(key)
        dest = stage / f"{key}.pdb"
        shutil.copy2(src, dest)
        meta[key] = {"source": str(src.relative_to(root)) if root in src.parents else src.name,
                     "sha256": sha256(src), "accession_alias": acc}
    return meta, aliases


def load_annotations(db: dict, db_path: Path) -> dict:
    annotations = {}
    for rel in db.get("annotation_tables", []):
        p = resolve(db_path, rel, db["path_base"])
        if not p.exists():
            continue
        for r in csv.DictReader(p.open(), delimiter="\t"):
            annotations[r.get("Entry", "")] = r
    return annotations


def build_proteome_universe(db: dict, db_path: Path, out: Path) -> tuple[list[dict], dict[str,dict]]:
    # Paths are recorded relative to the database root, matching stage_campaign
    # above. An absolute path here put the generating machine's home directory
    # into `proteome_denominator.json`, which is release evidence and is
    # published; it is also useless to a reader who cannot resolve it. The
    # sha256 beside it is what identifies the file, so nothing is lost.
    db_root = (db_path.parent / db["path_base"]).resolve()
    files = []
    for pattern in db["proteome_globs"]:
        base = resolve(db_path, pattern.split("*")[0], db["path_base"])
        root = base if base.is_dir() else base.parent
        files.extend(root.glob(Path(pattern).name) if "*" in pattern else [resolve(db_path, pattern, db["path_base"])])
    files = sorted({p.resolve() for p in files if p.exists()})
    records, index = [], {}
    with out.open("w") as fh:
        for pidx, path in enumerate(files):
            proteome_id = path.stem
            count = 0
            for rec in parse_fasta(path):
                uid = f"p{pidx:03d}|{rec['id']}"
                fh.write(f">{uid} {rec['header']}\n{rec['sequence']}\n")
                index[uid] = {**rec, "proteome_id": proteome_id, "source": str(path), "uid": uid}
                count += 1
            records.append({"proteome_id": proteome_id,
                            "path": str(path.relative_to(db_root)) if db_root in path.parents else path.name,
                            "protein_count": count, "sha256": sha256(path)})
    return records, index


def diamond_search(queries_fasta: Path, universe_fasta: Path, out: Path, work: Path, cfg: dict) -> list[dict]:
    dbprefix = work / "campaign_proteomes"
    run_cmd(["diamond","makedb","--in",str(universe_fasta),"--db",str(dbprefix),"--quiet"])
    fields = ["qseqid","sseqid","pident","length","qlen","slen","qcovhsp","evalue","bitscore","stitle"]
    run_cmd(["diamond","blastp","--query",str(queries_fasta),"--db",str(dbprefix),"--out",str(out),
             "--outfmt","6",*fields,"--evalue",str(cfg["sequence_evalue"]),
             "--id",str(cfg["sequence_min_identity"]),"--query-cover",str(cfg["sequence_min_query_coverage"]),
             "--max-target-seqs",str(cfg["sequence_max_hits"]),"--quiet"])
    rows = []
    with out.open() as fh:
        for vals in csv.reader(fh, delimiter="\t"):
            if len(vals) == len(fields): rows.append(dict(zip(fields, vals)))
    return rows


def reciprocal_best_hits(query: dict, grouped: dict[str,list[dict]], seq_index: dict[str,dict],
                         source_proteome: Path, work: Path) -> dict[str,str]:
    """Return target UID -> RBH status for the best hit from each comparison proteome."""
    candidates=[]
    for hits in grouped.values():
        if hits:
            hits.sort(key=lambda x:-float(x["bitscore"]))
            candidates.append(hits[0])
    if not candidates:
        return {}
    qf=work/f"rbh_{query['accession']}.faa"
    with qf.open("w") as fh:
        for h in candidates:
            rec=seq_index[h["sseqid"]]
            fh.write(f">{h['sseqid']}\n{rec['sequence']}\n")
    dbprefix=work/("rbh_source_"+hashlib.sha256(str(source_proteome).encode()).hexdigest()[:12])
    if not (Path(str(dbprefix)+".dmnd")).exists():
        run_cmd(["diamond","makedb","--in",str(source_proteome),"--db",str(dbprefix),"--quiet"])
    out=work/f"rbh_{query['accession']}.tsv"
    run_cmd(["diamond","blastp","--query",str(qf),"--db",str(dbprefix),"--out",str(out),
             "--outfmt","6","qseqid","sseqid","bitscore","--evalue","1e-5","--max-target-seqs","1","--quiet"])
    reverse={}
    if out.exists():
        with out.open() as fh:
            for row in csv.reader(fh,delimiter="\t"):
                if len(row)>=2: reverse[row[0]]=norm_id(row[1])
    return {h["sseqid"]:("reciprocal_best_hit" if reverse.get(h["sseqid"])==query["accession"] else "best_hit_nonreciprocal") for h in candidates}


def run_pipeline(query_manifest_path: Path, db_manifest_path: Path, output: Path) -> dict:
    qman, db = read_json(query_manifest_path), read_json(db_manifest_path)
    qroot = (query_manifest_path.parent / qman["path_base"]).resolve()
    dbroot = (db_manifest_path.parent / db["path_base"]).resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    work = output / "work"
    work.mkdir(exist_ok=True)
    versions = {"sf_csa":"1.0.0","foldseek":tool_version("foldseek"),"diamond":tool_version("diamond")}
    if db["required_foldseek_version"] not in versions["foldseek"]:
        raise SFCSError(f"Foldseek version drift: expected {db['required_foldseek_version']}, got {versions['foldseek']}")

    annotations = load_annotations(db, db_manifest_path)
    validated, query_by_id = [], {}
    all_queries = work / "queries.fasta"
    with all_queries.open("w") as qfh:
        for q in qman["queries"]:
            fasta = (qroot / q["fasta_path"]).resolve()
            pdb = (qroot / q["structure_path"]).resolve()
            if not fasta.exists() or not pdb.exists():
                raise SFCSError(f"missing query input for {q['accession']}")
            reject_reserved_fields(q, f"query {q.get('accession', '?')}")
            rec = select_fasta(fasta, q["accession"])
            if sequence_sha(rec["sequence"]) != q["sequence_sha256"]:
                raise SFCSError(f"FASTA checksum mismatch for {q['accession']}")
            if sha256(pdb) != q["structure_sha256"]:
                raise SFCSError(f"structure checksum mismatch for {q['accession']}")
            pseq, residues, coords = pdb_sequence(pdb, q.get("chain", "A"))
            exact = pseq == rec["sequence"]
            if not exact:
                raise SFCSError(f"structure sequence does not exactly match FASTA for {q['accession']}: {len(pseq)} vs {len(rec['sequence'])}")
            row = {**q,"sequence_length":len(rec["sequence"]),"structure_residue_count":len(residues),
                   "identity_status":"exact_full_sequence","geometry":geometry(coords),
                   "uniprot_annotation":annotations.get(q.get("uniprot_accession", ""), {})}
            validated.append(row); query_by_id[q["accession"]] = row
            qfh.write(f">{q['accession']}\n{rec['sequence']}\n")
    write_json(output / "validated_queries.json", validated)

    campaign_dir = work / "campaign_structures"
    campaign_meta, aliases = discover_campaign_structures(db, db_manifest_path, dbroot, campaign_dir)
    pdb_db = resolve(db_manifest_path, db["pdb_database"], db["path_base"])
    if not pdb_db.exists(): raise SFCSError(f"PDB Foldseek database absent: {pdb_db}")
    actual_database_checksum = database_bundle_checksum(pdb_db)
    if actual_database_checksum != db.get("pdb_database_checksum"):
        raise SFCSError(
            "PDB Foldseek database bundle checksum drift: "
            f"expected {db.get('pdb_database_checksum')}, got {actual_database_checksum}"
        )

    # The sequence leg runs first. Reciprocal best hits are evidence the
    # structural classification is entitled to use, so they have to exist before
    # it runs; computing them afterwards is what left `probable_same_function`
    # unreachable by measurement. See tests/test_rbh_provenance.py.
    universe_fasta = work / "campaign_proteomes.faa"
    proteomes, seq_index = build_proteome_universe(db, db_manifest_path, universe_fasta)
    seq_hits = diamond_search(all_queries, universe_fasta, work/"diamond_hits.tsv", work, db["thresholds"])
    grouped: dict[str,dict[str,list[dict]]] = defaultdict(lambda: defaultdict(list))
    available = set(campaign_meta) | set(aliases)
    for h in seq_hits:
        meta = seq_index.get(h["sseqid"])
        if not meta: continue
        qid = norm_id(h["qseqid"]); acc = norm_id(meta["id"])
        grouped[qid][meta["proteome_id"]].append({**h,"target_accession":acc,"proteome_id":meta["proteome_id"],
            "protein_header":meta["header"],"structure_status":"available_local_structure" if acc in available or aliases.get(acc) in available else "candidate_missing_structure"})

    # RBH is a relation between one query and one target, never a property of the
    # target on its own: a target may be reciprocal for one query and not for
    # another. Keyed by query accession for exactly that reason — a flag stored
    # on a shared target row would be query-independent by construction.
    rbh_by_query: dict[str,dict[str,str]] = {}
    rbh_targets: dict[str,set[str]] = {}
    for q in validated:
        source_proteome = (qroot / q["source_proteome_path"]).resolve()
        if not source_proteome.exists(): raise SFCSError(f"source proteome missing for RBH: {q['accession']}")
        hits = reciprocal_best_hits(q, grouped[q["accession"]], seq_index, source_proteome, work)
        rbh_by_query[q["accession"]] = hits
        rbh_targets[q["accession"]] = {
            norm_id(seq_index[s]["id"]) for s, status in hits.items()
            if status == "reciprocal_best_hit" and s in seq_index
        }

    structure_rows_by_query: dict[str,list[dict]] = defaultdict(list)
    target_meta = {q["accession"]: q for q in validated}
    for q in validated:
        qid, qpdb = q["accession"], (qroot / q["structure_path"]).resolve()
        qdir = output / "targets" / qid
        rawdir = qdir / "raw"; rawdir.mkdir(parents=True, exist_ok=True)
        searches = [("experimental_pdb", pdb_db), ("campaign_models", campaign_dir)]
        for dbname, target in searches:
            raw = rawdir / f"{dbname}.tsv"
            rows = foldseek_search(qpdb, target, raw, work/f"tmp_{qid}_{dbname}",
                                   db["thresholds"]["max_structure_hits"], db["thresholds"]["structure_evalue"],
                                   exhaustive=bool(db["thresholds"].get("exhaustive_structure_search", False)))
            for h in rows:
                category = structural_category(h, db["thresholds"]["whole_architecture_coverage"], db["thresholds"]["same_fold_tm"])
                tid = norm_id(h["target"])
                meta = target_meta.get(tid) or target_meta.get(aliases.get(tid, ""))
                fc,basis,ev,limit = classify_hit(
                    q, h, category, meta,
                    db.get("mechanism_families"),
                    db.get("contested_groups"),
                    db.get("divergence_sets"),
                    rbh=tid in rbh_targets.get(qid, frozenset()),
                )
                structure_rows_by_query[qid].append({
                    "query_id":qid,"target_id":tid,"database":dbname,"target_description":h.get("theader", ""),
                    "evalue":h.get("evalue", ""),"fident":h.get("fident", ""),"aligned_length":h.get("alnlen", ""),
                    "query_length":h.get("qlen", ""),"target_length":h.get("tlen", ""),
                    "query_coverage":h.get("qcov", ""),"target_coverage":h.get("tcov", ""),
                    "query_tm_score":h.get("qtmscore", ""),"target_tm_score":h.get("ttmscore", ""),
                    "aligned_tm_score":h.get("alntmscore", ""),"rmsd":h.get("rmsd", ""),
                    "structural_category":category,"function_classification":fc,"classification_basis":basis,
                    "evidence_class":ev,"limitation":limit,"transform_u":h.get("u", ""),"transform_t":h.get("t", ""),
                })
        rows = sorted(structure_rows_by_query[qid], key=lambda r: (r["database"], -float(r["aligned_tm_score"] or 0)))
        write_tsv(qdir/"structure_hits.tsv", rows, TSV_FIELDS + ["transform_u","transform_t"])
        function_rows = [{k:r[k] for k in ["query_id","target_id","database","function_classification","classification_basis","evidence_class","limitation"]} for r in rows]
        write_tsv(qdir/"functional_evidence.tsv", function_rows, list(function_rows[0]) if function_rows else ["query_id"])
        overlays=[r for r in rows if r["database"]=="campaign_models" and r["target_id"] in target_meta and r["target_id"]!=qid and r["structural_category"]!="below_structural_similarity_threshold"]
        if overlays:
            best=max(overlays,key=lambda r:float(r["aligned_tm_score"] or 0))
            target_source=campaign_dir/f"{best['target_id']}.pdb"
            aligned=qdir/f"{best['target_id']}_aligned_to_{qid}.pdb"
            transform_pdb(target_source,aligned,best["transform_u"],best["transform_t"])
            shutil.copy2(qpdb,qdir/f"{qid}_query.pdb")
            write_superposition_html(qdir/"superposition.html",qdir/f"{qid}_query.pdb",aligned,qid,best["target_id"],best["aligned_tm_score"],best["rmsd"])

    all_species_rows = []
    for q in validated:
        rbh = rbh_by_query[q["accession"]]
        rows = []
        for proteome_id, hits in grouped[q["accession"]].items():
            hits.sort(key=lambda x: -float(x["bitscore"]))
            for rank,h in enumerate(hits[:db["thresholds"]["sequence_hits_per_proteome"]],1):
                h["within_proteome_rank"] = rank
                if rank == 1:
                    h["orthology_status"] = rbh.get(h["sseqid"],"best_hit_unresolved")
                else:
                    h["orthology_status"] = "paralog_or_secondary_candidate"
                h["functional_interpretation"] = ("ortholog_candidate_requires_independent_function_and_structure_evidence"
                    if h["orthology_status"]=="reciprocal_best_hit" else "unresolved_paralog_or_nonreciprocal_candidate")
                rows.append(h)
        fields = ["qseqid","target_accession","proteome_id","within_proteome_rank","orthology_status","pident","length","qlen","slen","qcovhsp","evalue","bitscore","structure_status","functional_interpretation","protein_header"]
        write_tsv(output/"targets"/q["accession"]/"species_comparison.tsv", rows, fields)
        all_species_rows.extend(rows)

        q_struct = structure_rows_by_query[q["accession"]]
        assembly = {
            "query_id":q["accession"],"structure_class":q["structure_class"],"query_state":"predicted_monomer",
            "experimental_exact_assembly":"not_identified_in_current_local_evidence",
            "ligand_state":"not_resolved_from_predicted_monomer",
            "partners_and_stoichiometry":"not_resolved_from_predicted_monomer",
            "closest_experimental_hits":[r for r in q_struct if r["database"]=="experimental_pdb"][:5],
            "boundary":"Foldseek headers do not establish biological assembly, ligand state, partner identity or active pose.",
        }
        write_json(output/"targets"/q["accession"]/"assembly_and_state.json", assembly)
        write_dossier(output/"targets"/q["accession"]/"DOSSIER.md", q, q_struct, rows)

    matrix = []
    for q in validated:
        for r in structure_rows_by_query[q["accession"]]:
            if r["database"] == "campaign_models" and r["target_id"] in target_meta:
                matrix.append(r)
    write_tsv(output/"RELEASE_COMPARISON_MATRIX.tsv", matrix, TSV_FIELDS)
    write_json(output/"proteome_denominator.json", {"proteome_count":len(proteomes),"protein_count":sum(x["protein_count"] for x in proteomes),"proteomes":proteomes})

    manifest = {
        "schema_version":1,"release_id":output.name,"status":"COMPUTATIONAL_STRUCTURE_FUNCTION_COMPARISON_NOT_EXPERIMENTAL_VALIDATION",
        "queries":[q["accession"] for q in validated],"query_count":len(validated),"target_statuses":{q["accession"]:q["decision_status"] for q in validated},
        "tools":versions,"database_manifest_sha256":sha256(db_manifest_path),"query_manifest_sha256":sha256(query_manifest_path),
        "pdb_database":{"path_record":db["pdb_database"],"version":db["pdb_database_version"],"checksum":db["pdb_database_checksum"]},
        "campaign_structure_count":len(campaign_meta),"proteome_count":len(proteomes),"protein_count":sum(x["protein_count"] for x in proteomes),
        "thresholds":db["thresholds"],"classification_vocabulary":db["classification_vocabulary"],
        "scientific_boundaries":["fold similarity does not establish direct functional equivalence","predicted orientation is not native surface-exposure proof","proteins without structures are missing evidence, not structural negatives","target decisions are not changed automatically"],
    }
    write_json(output/"SF_CSA_RELEASE_MANIFEST.json", manifest)
    files = {str(p.relative_to(output)):sha256(p) for p in sorted(output.rglob("*")) if p.is_file() and "work" not in p.parts and p.name != "CHECKSUMS.json"}
    write_json(output/"CHECKSUMS.json", files)
    return manifest


def write_dossier(path: Path, q: dict, structural: list[dict], species: list[dict]) -> None:
    top_pdb = [r for r in structural if r["database"]=="experimental_pdb"][:5]
    top_campaign = [r for r in structural if r["database"]=="campaign_models" and r["target_id"] != q["accession"]][:5]
    lines = [f"# {q['common_name']} — SF-CSA dossier", "",
        f"- Exact target: `{q['accession']}`; {q['organism']} {q['strain']}",
        f"- Decision: **{q['decision_status']}** (unchanged by this analysis)",
        f"- Mechanism group: `{q['mechanism_group']}`", f"- Sequence SHA-256: `{q['sequence_sha256']}`",
        f"- Structure SHA-256: `{q['structure_sha256']}`", f"- Structure class: `{q['structure_class']}`", "",
        "## Closest experimental PDB structures", "",
        "| Target | aln TM | qcov | tcov | classification | limitation |", "|---|---:|---:|---:|---|---|"]
    for r in top_pdb:
        lines.append(f"| `{r['target_id']}` | {r['aligned_tm_score']} | {r['query_coverage']} | {r['target_coverage']} | {r['function_classification']} | {r['limitation']} |")
    lines += ["", "## Closest campaign structures", "", "| Target | aln TM | identity | classification |", "|---|---:|---:|---|"]
    for r in top_campaign:
        lines.append(f"| `{r['target_id']}` | {r['aligned_tm_score']} | {r['fident']} | {r['function_classification']} |")
    missing = sum(r["structure_status"]=="candidate_missing_structure" for r in species)
    lines += ["", "## Species comparison", "", f"Sequence candidates retained: **{len(species)}**; candidates without a local structure: **{missing}**.", "",
        "## Interpretation boundary", "", q["protein_specific_boundary"], "",
        "A shared fold is an architectural hypothesis. It does not transfer substrate, activity, virulence, native exposure, immunogenicity or vaccine suitability.", ""]
    path.write_text("\n".join(lines))


def verify_release(output: Path, databases: Path | None = None) -> list[str]:
    """Audit a finished release against expectations held OUTSIDE it.

    `databases` is the database manifest the run was configured with. Its
    `release_expectations` block is what the release is checked against -- the
    counts deliberately do not come from the release manifest, because a release
    that was allowed to state its own expected shape would agree with itself no
    matter what it produced.
    """
    errors = []
    manifest_p = output/"SF_CSA_RELEASE_MANIFEST.json"
    if not manifest_p.exists(): return ["release manifest missing"]
    m = read_json(manifest_p)

    expected = read_json(databases).get("release_expectations", {}) if databases else None
    if expected is None:
        errors.append("no database manifest supplied: release shape was not verified, only its contents")
    else:
        if "query_count" in expected and m.get("query_count") != expected["query_count"]:
            errors.append(f"expected {expected['query_count']} targets, found {m.get('query_count')}")
        if "proteome_count" in expected and m.get("proteome_count") != expected["proteome_count"]:
            errors.append(f"expected {expected['proteome_count']} proteomes, found {m.get('proteome_count')}")
        for accession, status in (expected.get("target_statuses") or {}).items():
            if m.get("target_statuses",{}).get(accession) != status:
                errors.append(f"{accession} is not {status}")

    traps = (expected or {}).get("title_traps") or DEFAULT_TITLE_TRAPS
    for q in m.get("queries", []):
        qdir=output/"targets"/q
        for name in ("structure_hits.tsv","functional_evidence.tsv","species_comparison.tsv","assembly_and_state.json","DOSSIER.md"):
            if not (qdir/name).exists(): errors.append(f"{q}/{name} missing")
        p=qdir/"structure_hits.tsv"
        if p.exists():
            for r in csv.DictReader(p.open(),delimiter="\t"):
                description = r.get("target_description","").lower()
                for trap in traps:
                    if trap["substring"] in description and r.get("function_classification") in set(trap["must_not_promote_to"]):
                        errors.append(f"PDB title trap promoted for {q}: {trap['substring']}")
    checks=read_json(output/"CHECKSUMS.json") if (output/"CHECKSUMS.json").exists() else {}
    for rel,digest in checks.items():
        p=output/rel
        if not p.exists() or sha256(p)!=digest: errors.append(f"checksum mismatch: {rel}")
    return errors
