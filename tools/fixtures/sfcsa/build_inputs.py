#!/usr/bin/env python3
"""Generate the synthetic inputs for the sf-csa offline fixture.

Everything here is synthetic: invented accessions, invented organisms, sequences
built from a repeating pattern, coordinates on a helix. Nothing resembles a real
protein and nothing is downloaded.

The reason this is a generator rather than committed files: the query manifest
carries SHA-256 checksums of the FASTA sequence and the PDB file, and the
pipeline refuses to run if they do not match. Hand-writing those is a guaranteed
source of "fixture is broken and nobody knows why". They are computed here.

    python build_inputs.py            # writes inputs/ next to this file
    python build_inputs.py --check    # verify committed inputs are reproducible

Stdlib only, offline.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"

AA = "ACDEFGHIKLMNPQRSTVWY"

# Two synthetic queries. QRY_A and QRY_B are given the same mechanism_group so
# the fixture can exercise a same-group comparison; their sequences differ.
QUERIES = [
    {
        "accession": "QRY_A",
        "common_name": "synthetic barrel A",
        "organism": "Synthetica exempli",
        "strain": "STUB-1",
        "mechanism_group": "synthetic_barrel",
        "decision_status": "retained_for_review",
        "structure_class": "predicted_monomer",
        "uniprot_accession": "SYNA00001",
        "protein_specific_boundary": (
            "Synthetic fixture protein. No substrate, activity or localisation "
            "claim is possible and none is made."
        ),
        "length": 60,
        "seed": 11,
    },
    {
        "accession": "QRY_B",
        "common_name": "synthetic barrel B",
        "organism": "Synthetica exempli",
        "strain": "STUB-2",
        "mechanism_group": "synthetic_barrel",
        "decision_status": "retained_for_review",
        "structure_class": "predicted_monomer",
        "uniprot_accession": "SYNB00001",
        "protein_specific_boundary": (
            "Synthetic fixture protein. No substrate, activity or localisation "
            "claim is possible and none is made."
        ),
        "length": 55,
        "seed": 29,
    },
]


def synthetic_sequence(length: int, seed: int) -> str:
    """Deterministic pseudo-sequence. No PRNG, so it is stable across versions."""
    return "".join(AA[(seed * (index + 1) + index * index) % len(AA)] for index in range(length))


THREE = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR",
}


def synthetic_pdb(sequence: str, name: str) -> str:
    """CA-only helix. The pipeline reads the CA trace and the residue names.

    `pdb_sequence` requires the structure's residue sequence to match the FASTA
    exactly, so this is generated from the same string.
    """
    lines = [f"HEADER    SYNTHETIC FIXTURE STRUCTURE {name}"]
    for index, letter in enumerate(sequence, start=1):
        angle = index * 100.0 * math.pi / 180.0
        x, y, z = 2.3 * math.cos(angle), 2.3 * math.sin(angle), 1.5 * index
        lines.append(
            f"ATOM  {index:5d}  CA  {THREE[letter]} A{index:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 25.00           C"
        )
    lines += ["TER", "END"]
    return "\n".join(lines) + "\n"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sequence_sha(sequence: str) -> str:
    """Must match `sf_csa.core.sequence_sha` exactly: ASCII-encoded, then hashed.

    Read from the source rather than assumed — an earlier draft here upper-cased
    first, which would have produced a manifest the pipeline rejects.
    """
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def database_bundle_checksum(path: Path) -> str:
    """Match sf_csa.core's deterministic checksum for a directory bundle."""
    files = sorted(item for item in path.rglob("*") if item.is_file())
    if not files:
        raise ValueError(f"synthetic structure database is empty: {path}")
    payload = "".join(
        f"{item.relative_to(path).as_posix()}\t{sha256_file(item)}\n"
        for item in files
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def build(root: Path) -> None:
    if root.exists():
        shutil.rmtree(root)
    queries_dir = root / "queries"
    campaign_dir = root / "campaign_structures"
    proteomes_dir = root / "proteomes"
    pdbdb_dir = root / "pdb_database"
    for directory in (queries_dir, campaign_dir, proteomes_dir, pdbdb_dir):
        directory.mkdir(parents=True)

    # The PDB "database" is a path the pipeline only checks for existence; the
    # stub aligner never reads it.
    (pdbdb_dir / "README.txt").write_text(
        "Synthetic stand-in for a Foldseek PDB database. The stub aligner in\n"
        "stub_bin/ never reads this; the pipeline only checks that it exists.\n",
        encoding="utf-8",
    )

    manifest_queries = []
    for spec in QUERIES:
        sequence = synthetic_sequence(spec["length"], spec["seed"])
        accession = spec["accession"]

        fasta = queries_dir / f"{accession}.faa"
        fasta.write_text(f">{accession} {spec['common_name']}\n{sequence}\n", encoding="utf-8")

        pdb = queries_dir / f"{accession}.pdb"
        pdb.write_text(synthetic_pdb(sequence, accession), encoding="utf-8")

        # The same structure is staged as a campaign model, which is what lets a
        # query appear as another query's target in the comparison matrix.
        shutil.copy2(pdb, campaign_dir / f"{accession}.pdb")

        # Each query's source proteome, used for the reverse RBH search.
        source = queries_dir / f"{accession}_source_proteome.faa"
        source.write_text(f">{accession} {spec['common_name']}\n{sequence}\n", encoding="utf-8")

        manifest_queries.append(
            {
                "accession": accession,
                "common_name": spec["common_name"],
                "organism": spec["organism"],
                "strain": spec["strain"],
                "mechanism_group": spec["mechanism_group"],
                "decision_status": spec["decision_status"],
                "structure_class": spec["structure_class"],
                "uniprot_accession": spec["uniprot_accession"],
                "protein_specific_boundary": spec["protein_specific_boundary"],
                "chain": "A",
                "fasta_path": f"queries/{accession}.faa",
                "structure_path": f"queries/{accession}.pdb",
                "source_proteome_path": f"queries/{accession}_source_proteome.faa",
                "sequence_sha256": sequence_sha(sequence),
                "structure_sha256": sha256_file(pdb),
            }
        )

    # Two comparison proteomes. Sequences are decoys the canned DIAMOND hits can
    # point at, plus -- in proteome 1 -- QRY_B itself.
    #
    # QRY_B is there so a reciprocal best hit can land on a protein that is also
    # a structural target. `probable_same_function` needs both legs on the same
    # target: an RBH from the sequence search and a whole-architecture match from
    # the structure search. Before 2026-09-01 nothing in the pipeline wrote the
    # RBH flag at all, so the trap scenario declared it by hand in the manifest.
    # That route is now refused at manifest read time, and this record is what
    # replaces it: the label is reached by computation, which is what the panel
    # is supposed to be testing.
    for proteome_id, seeds in (("SYN_PROTEOME_1", (11, 29)), ("SYN_PROTEOME_2", (41,))):
        records = []
        for seed in seeds:
            length = 60 if seed == 11 else 55 if seed == 29 else 48
            records.append(
                f">SYN{seed:03d} synthetic protein seed {seed} [{proteome_id}]\n"
                f"{synthetic_sequence(length, seed)}\n"
            )
        if proteome_id == "SYN_PROTEOME_1":
            spec = next(q for q in QUERIES if q["accession"] == "QRY_B")
            records.append(
                f">QRY_B {spec['common_name']} [{proteome_id}]\n"
                f"{synthetic_sequence(spec['length'], spec['seed'])}\n"
            )
        (proteomes_dir / f"{proteome_id}.faa").write_text("".join(records), encoding="utf-8")

    # A UniProt-shaped annotation table. Only `Entry` is keyed on by the pipeline.
    annotations = proteomes_dir / "annotations.tsv"
    annotations.write_text(
        "Entry\tProtein names\tCofactor\tEC number\n"
        "SYNA00001\tSynthetic barrel A\t\t\n"
        "SYNB00001\tSynthetic barrel B\t\t\n",
        encoding="utf-8",
    )

    (root / "query_manifest.json").write_text(
        json.dumps(
            {"path_base": ".", "campaign": "synthetic-fixture", "queries": manifest_queries},
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    database_manifest = {
        "path_base": ".",
        "required_foldseek_version": "0.0.0-stub",
        "pdb_database": "pdb_database",
        "pdb_database_version": "synthetic-fixture-0",
        "pdb_database_checksum": database_bundle_checksum(pdbdb_dir),
        "campaign_structure_roots": ["campaign_structures"],
        "proteome_globs": ["proteomes/*.faa"],
        "annotation_tables": ["proteomes/annotations.tsv"],
        "seqmatch_tables": [],
        "thresholds": {
            "max_structure_hits": 50,
            "structure_evalue": "1e-3",
            "whole_architecture_coverage": 0.8,
            "same_fold_tm": 0.5,
            "sequence_evalue": "1e-5",
            "sequence_min_identity": 30,
            "sequence_min_query_coverage": 60,
            "sequence_max_hits": 25,
            "sequence_hits_per_proteome": 5,
        },
        "classification_vocabulary": [
            "exact_function_supported",
            "probable_same_function",
            "same_mechanism_class",
            "structural_analogy_only",
            "candidate_functional_divergence",
            "unresolved_or_conflicted",
        ],
        # Shape read from `sf_csa.core.DEFAULT_MECHANISM_FAMILIES`: each entry is
        # {"group", "pattern"} with an optional "refine" list of narrower
        # {"group", "pattern"} entries. `classify_title` regex-searches the
        # lower-cased PDB title in list order and returns the first match, so
        # ordering is part of the contract and the fixture exercises it: a title
        # matching both entries below must come back as the narrower group.
        "mechanism_families": [
            {
                "group": "synthetic_barrel",
                "pattern": r"synthetic barrel",
                "refine": [{"group": "synthetic_partial_domain", "pattern": r"partial domain"}],
            },
            {"group": "synthetic_transporter", "pattern": r"transporter"},
        ],
        "contested_groups": [],
        "divergence_sets": [],
        # Expectations live here, outside the release, because verify_release
        # refuses to let a release state its own expected shape.
        "release_expectations": {
            "query_count": 2,
            "proteome_count": 2,
            "target_statuses": {"QRY_A": "retained_for_review", "QRY_B": "retained_for_review"},
        },
    }
    (root / "database_manifest.json").write_text(
        json.dumps(database_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    # ---- Second scenario: a release that must FAIL its own audit. -------------
    #
    # A fixture whose audit always returns clean does not demonstrate that the
    # audit works. This pair of manifests drives one hit up to
    # `probable_same_function` while its PDB title carries a trap substring, so
    # `verify_release` has something to catch.
    #
    # Rewritten 2026-09-01. The label used to be reached by declaring `rbh: true`
    # on QRY_B by hand, because nothing in the pipeline ever wrote the key. Both
    # halves of that changed: the RBH computation now runs before classification
    # and reaches the label as a pairwise fact, and `reject_reserved_fields`
    # refuses a curator-supplied `rbh` at manifest read time. The hand-declared
    # route is not merely unnecessary now, it is rejected -- so the trap manifest
    # is an ordinary manifest, and the promotion it traps is one the module made
    # on its own evidence.
    #
    # That is a stronger test than the one it replaces. A hand-declared promotion
    # is caught upstream by the manifest reader; a computed one is not, and the
    # title trap is genuinely the last line of defence against it.
    trap_queries = [dict(entry) for entry in manifest_queries]
    (root / "query_manifest_trap.json").write_text(
        json.dumps(
            {
                "path_base": ".",
                "campaign": "synthetic-fixture-trap-scenario",
                "queries": trap_queries,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    trap_database = json.loads(json.dumps(database_manifest))
    trap_database["release_expectations"] = {
        "query_count": 2,
        "proteome_count": 2,
        "target_statuses": {
            "QRY_A": "retained_for_review",
            "QRY_B": "retained_for_review",
        },
        # Declared explicitly rather than relying on DEFAULT_TITLE_TRAPS, so the
        # scenario states what it is testing.
        "title_traps": [
            {
                "substring": "toluene",
                "must_not_promote_to": ["exact_function_supported", "probable_same_function"],
            }
        ],
    }
    (root / "database_manifest_trap.json").write_text(
        json.dumps(trap_database, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def ensure_inputs() -> Path:
    """Generate `inputs/` if it is absent, and return its path.

    The inputs are generated rather than committed, so anything that reads them
    -- the fixture tests, the composition example, `run_fixture.sh` -- has to be
    able to bring them into existence. Doing that here rather than in each
    caller keeps one definition of "the inputs are built by this script". This
    deliberately does *not* rebuild an existing tree: `--check` is the way to
    detect drift, and silently regenerating over a modified tree would hide it.
    """
    if not (INPUTS / "query_manifest.json").exists():
        build(INPUTS)
    return INPUTS


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="rebuild into a scratch directory and compare against the committed inputs",
    )
    args = parser.parse_args(argv)

    if not args.check:
        build(INPUTS)
        count = sum(1 for p in INPUTS.rglob("*") if p.is_file())
        print(f"wrote {count} file(s) to {INPUTS}")
        return 0

    import tempfile

    with tempfile.TemporaryDirectory() as scratch:
        fresh = Path(scratch) / "inputs"
        build(fresh)
        committed = {
            str(p.relative_to(INPUTS)): sha256_file(p)
            for p in sorted(INPUTS.rglob("*"))
            if p.is_file()
        }
        rebuilt = {
            str(p.relative_to(fresh)): sha256_file(p)
            for p in sorted(fresh.rglob("*"))
            if p.is_file()
        }

    if committed == rebuilt:
        print(f"inputs are reproducible: {len(committed)} file(s) identical")
        return 0
    for name in sorted(set(committed) | set(rebuilt)):
        if committed.get(name) != rebuilt.get(name):
            print(f"DIFFERS  {name}")
    print("\ninputs are NOT reproducible from build_inputs.py")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
