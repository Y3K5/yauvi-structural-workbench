#!/usr/bin/env python3
"""Assemble the sf-csa qualification campaign from the frozen selection.

Everything here is derived from `SF_CSA_RECORD_SELECTION.md`: twelve AlphaFold
models, their UniProt sequences, their organisms' reference proteomes, and a
structure database of the twelve crystal entries whose SCOP classification
assigned the strata.

The one decision that is not mechanical is `mechanism_group`, and it decides
whether the panel measures anything. `classify_hit` promotes only when query and
target share a mechanism group, so:

    group at SCOP *fold* granularity  -> a fold analogy shares a group with its
                                         reference and is reported same_mechanism_class,
                                         which is exactly wrong
    group at SCOP *superfamily*       -> homologs share a group and analogues do
                                         not, which is the distinction under test

So the group is the SCOP superfamily, named after it. This also discharges the
Finding 2 precondition: the periodontal default tables are replaced outright, and
`build-manifests` writes the replacement into the database manifest by value.
"""
from __future__ import annotations

import hashlib, json, re, shutil, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# accession -> (pdb entry, SCOP superfamily group, reference proteome)
#
# Reselected 2026-09-02. The previous twelve were chosen against a gate that
# checked whether an entry's *organism* had a reference proteome, not whether the
# entry is *in* one. Five were in none at all -- P00193, P23370, P45850, P00138,
# P00147 -- and a sixth, P26394, was pointed at UP000002695, which loads fine and
# contains zero occurrences of it. The sequence leg silently vanished for thirteen
# of sixteen records, including two exact self-matches, which is what exposed it.
#
# Every entry below passes five gates: one unambiguous SCOP domain; stratum
# derived from SCOP; an AlphaFold model sequence-exact against its UniProt FASTA;
# a reference proteome; and -- new, and the one that was missing -- the accession
# appears in that proteome. Verified 12/12 before this table was written.
QUERIES = {
 "P00198": ("1FDN", "scop_d58_1_ferredoxin_4fe4s",            "UP000006094"),
 "P00208": ("1BLU", "scop_d58_1_ferredoxin_4fe4s",            "UP000001441"),
 "P00818": ("1APS", "scop_d58_10_acylphosphatase",            "UP000002281"),
 "P32081": ("1CSP", "scop_b40_4_nucleic_acid_binding",        "UP000001570"),
 "P0A9X9": ("1MJC", "scop_b40_4_nucleic_acid_binding",        "UP000000625"),
 "P00817": ("1E6A", "scop_b40_5_inorganic_pyrophosphatase",   "UP000002311"),
 "P26394": ("1DZR", "scop_b82_1_rmlc_cupins",                 "UP000001014"),
 "O00625": ("1J1L", "scop_b82_1_rmlc_cupins",                 "UP000005640"),
 "P37610": ("1OS7", "scop_b82_2_clavaminate_synthase_like",   "UP000000625"),
 "Q07688": ("1C02", "scop_a24_10_hpt_domain",                 "UP000002311"),
 "Q9A980": ("2OOC", "scop_a24_10_hpt_domain",                 "UP000001816"),
 "P43934": ("1JOG", "scop_a24_16_nucleotidyltransferase_substrate_binding", "UP000000579"),
}

# One regex per group, matched against target titles. Only the experimental leg
# needs these; on the campaign axis both sides carry the curated group. They are
# stated anyway, because leaving them out silently reinstates the periodontal
# defaults for every experimental hit.
FAMILIES = [
 {"group": "scop_d58_1_ferredoxin_4fe4s",          "pattern": r"ferredoxin|4fe-4s|4fe4s"},
 {"group": "scop_d58_10_acylphosphatase",          "pattern": r"acylphosphatase|acyl[- ]phosphatase|\bacyp\b"},
 {"group": "scop_b40_4_nucleic_acid_binding",      "pattern": r"cold[- ]shock|nucleic acid[- ]binding"},
 {"group": "scop_b40_5_inorganic_pyrophosphatase", "pattern": r"inorganic pyrophosphatase|ppase"},
 {"group": "scop_b82_1_rmlc_cupins",               "pattern": r"cupin|germin|oxalate oxidase|epimerase|rmlc"},
 {"group": "scop_b82_2_clavaminate_synthase_like", "pattern": r"taurine dioxygenase|clavaminate|tfda|\btaud\b"},
 {"group": "scop_a24_16_nucleotidyltransferase_substrate_binding", "pattern": r"nucleotidyltransferase|kanamycin nucleotidyl|\bknt\b"},
 {"group": "scop_a24_10_hpt_domain",               "pattern": r"histidine[- ]containing phosphotransfer|\bhpt\b|ypd1"},
]

BOUNDARY = ("Shared SCOP superfamily is evidence of common ancestry, not of shared "
            "substrate or activity. Fold similarity alone transfers nothing.")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def header_of(fasta: Path) -> str:
    return fasta.read_text(encoding="utf-8").splitlines()[0]


def organism_of(header: str) -> tuple[str, str]:
    """(organism, description) from a UniProt FASTA header."""
    m = re.search(r"OS=(.+?)\s+(?:OX=|GN=|PE=|SV=)", header)
    org = m.group(1) if m else "unknown"
    d = re.match(r">\S+\s+(.*?)\s+OS=", header)
    return org, (d.group(1) if d else "unknown")


def main() -> int:
    models = HERE / "sources/alphafold"
    fastas = HERE / "sources/uniprot"
    campaign = HERE / "sources/campaign_models"
    results = HERE / "results/sf_csa"
    campaign.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)

    seq_rows, ledger_rows, targets, missing = [], [], [], []
    for acc, (entry, group, upid) in QUERIES.items():
        model = next(iter(models.glob(f"AF-{acc}-F1-model_v*.pdb")), None)
        fasta = fastas / f"{acc}.fasta"
        proteome = HERE / f"sources/proteomes/{upid}.fasta"
        for label, path in (("model", model), ("fasta", fasta), ("proteome", proteome)):
            if path is None or not path.exists():
                missing.append(f"{acc} {label}")
        if model is None or not fasta.exists():
            continue

        # Campaign structures are discovered by filename: norm_id() reads the stem,
        # so the model must be named for its accession, not for its AlphaFold file.
        shutil.copy2(model, campaign / f"{acc}.pdb")

        header = header_of(fasta)
        organism, description = organism_of(header)
        seq = "".join(l.strip() for l in fasta.read_text().splitlines()[1:])
        seq_rows.append((acc, f"sources/uniprot/{acc}.fasta",
                         hashlib.sha256(seq.encode("ascii")).hexdigest()))
        ledger_rows.append((acc, "SELECTED"))
        targets.append({
            "accession": acc, "common_name": description, "organism": organism,
            "strain": "reference_proteome:" + upid, "uniprot_accession": acc,
            "mechanism_group": group, "protein_specific_boundary": BOUNDARY,
            "structure_path": f"sources/campaign_models/{acc}.pdb",
            "source_proteome_path": f"sources/proteomes/{upid}.fasta",
            "pdb_entry": entry,
        })

    (results / "SEQUENCE_MANIFEST.tsv").write_text(
        "accession\tsource_file\tsequence_sha256\n" +
        "".join("\t".join(r) + "\n" for r in seq_rows), encoding="utf-8")
    (results / "SELECTION_LEDGER.tsv").write_text(
        "primary_accession\tdecision_status\n" +
        "".join("\t".join(r) + "\n" for r in ledger_rows), encoding="utf-8")
    (results / "orient_manifest.json").write_text("{}\n", encoding="utf-8")
    (HERE / "sources/structure_db.version").write_text(
        "qualification-v2-sf-csa-structdb-1\n"
        "twelve wwPDB entries whose SCOP classification assigned the panel strata\n",
        encoding="utf-8")

    spec = {
        "schema_version": 1,
        "release_scope": "sf-csa qualification panel: four SCOP folds, four relationship strata each",
        "root": ".",
        "sequence_manifest": "results/sf_csa/SEQUENCE_MANIFEST.tsv",
        "decision_ledger": "results/sf_csa/SELECTION_LEDGER.tsv",
        "default_orientation_artifact": "results/sf_csa/orient_manifest.json",
        "targets": targets,
        "database": {
            # run_pipeline enforces the Foldseek pin as a substring check and aborts on
            # drift. It does not check DIAMOND at all, so the DIAMOND pin below is a
            # record, not a guard -- stated because an unenforced pin that looks
            # enforced is worse than an absent one.
            "required_foldseek_version": "10.941cd33",
            "required_diamond_version_unenforced": "2.1.11",
            "pdb_database": "sources/structure_db",
            "pdb_database_version": "qualification-v2-sf-csa-structdb-1",
            "campaign_structure_roots": ["sources/campaign_models"],
            "proteome_globs": ["sources/proteomes/*.fasta"],
            "thresholds": {
                # Collection 2.6. Two search filters act in series and only the
                # first is visible: Foldseek prefilters before scoring and the
                # prefilter ignores the e-value, so a distant pair produces no
                # row at any e-value. Measured for P00198 against this database:
                # e=0.01 -> 2 rows, e=10000 -> 2 rows, e=0.01 exhaustive -> 2
                # rows, e=10000 exhaustive -> 12. The panel gates on the
                # difference between a pair rejected below threshold and a pair
                # never compared, so both are set and the e-value stops acting
                # as a hidden scientific gate. Filtering is left to same_fold_tm
                # and whole_architecture_coverage, which the panel declares.
                "structure_evalue": "10000", "exhaustive_structure_search": True,
                "same_fold_tm": 0.5,
                "whole_architecture_coverage": 0.7, "max_structure_hits": 50,
                "sequence_evalue": "1e-5", "sequence_min_identity": 20,
                "sequence_min_query_coverage": 40, "sequence_max_hits": 5000,
                "sequence_hits_per_proteome": 10,
            },
            "mechanism_families": FAMILIES,
            "contested_groups": [],
            "divergence_sets": [],
            "coverage_boundary":
                "Search is exhaustive only within the twelve-entry structure database and the "
                "ten reference proteomes named here. A protein without a model is reported as "
                "missing evidence, never as a structural negative.",
        },
    }
    spec_path = HERE / "sf_csa_campaign.json"
    spec_path.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")

    print(f"targets: {len(targets)}   groups: {len({t['mechanism_group'] for t in targets})}")
    if missing:
        print(f"MISSING ({len(missing)}): " + ", ".join(missing[:12]))
    print(f"spec: {spec_path.name}")
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
