#!/usr/bin/env python3
"""Show that every sf-csa gate can fail, before any record is curated.

Adoption protocol rule 2: a gate that cannot fail is not a gate. The membrane
panel satisfied every other procedural rule while measuring the wrong quantity,
so the gates are exercised here against a synthetic release whose rows are
written by hand -- one variant that must pass, and for each gate a variant that
must fail. Nothing here touches Foldseek, DIAMOND or a real structure: the
subject under test is the gate logic, not the aligners.

    python3 sf_csa_gate_falsification.py

Exit 0 when every gate both passes on the good release and fails on its
corruption. Exit 1 if any gate is unfalsifiable.
"""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
from pathlib import Path

import run_execution as R

DB_SHA = "f" * 64
STRUCT_FIELDS = ["query_id", "target_id", "database", "target_description", "evalue", "fident",
                 "aligned_length", "query_length", "target_length", "query_coverage",
                 "target_coverage", "query_tm_score", "target_tm_score", "aligned_tm_score",
                 "rmsd", "structural_category", "function_classification",
                 "classification_basis", "evidence_class", "limitation"]
SEQ_FIELDS = ["qseqid", "target_accession", "proteome_id", "within_proteome_rank",
              "orthology_status", "pident", "length", "qlen", "slen", "qcovhsp", "evalue",
              "bitscore", "structure_status", "functional_interpretation", "protein_header"]

GATE_SEMANTICS = {"sf_csa": {"database_manifest_sha256": DB_SHA}}
VOCAB = ["exact_function_supported", "probable_same_function", "same_mechanism_class",
         "structural_analogy_only", "candidate_functional_divergence", "unresolved_or_conflicted"]


def write_release(root: Path, *, query="Q1", target="T1", category="whole_architecture_match",
                  label="same_mechanism_class", orthology="reciprocal_best_hit",
                  db_sha=DB_SHA, structural=True, sequence=True) -> Path:
    out = root / "release"
    qdir = out / "targets" / query
    qdir.mkdir(parents=True, exist_ok=True)

    def tsv(path, fields, row):
        lines = ["\t".join(fields)]
        if row is not None:
            lines.append("\t".join(str(row.get(f, "")) for f in fields))
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    tsv(qdir / "structure_hits.tsv", STRUCT_FIELDS, {
        "query_id": query, "target_id": target, "database": "campaign_models",
        "structural_category": category, "function_classification": label,
        "aligned_tm_score": "0.82", "query_coverage": "0.91", "target_coverage": "0.90",
    } if structural else None)
    tsv(qdir / "species_comparison.tsv", SEQ_FIELDS, {
        "qseqid": query, "target_accession": target, "orthology_status": orthology,
        "pident": "61.4", "evalue": "1e-80",
    } if sequence else None)
    (out / "SF_CSA_RELEASE_MANIFEST.json").write_text(json.dumps({
        "schema_version": 1, "classification_vocabulary": VOCAB,
        "database_manifest_sha256": db_sha, "query_manifest_sha256": "a" * 64,
        "tools": {"foldseek": "10.941cd33", "diamond": "2.1.11"},
    }), encoding="utf-8")
    return out


def judge(out: Path, *, stratum, label, category="whole_architecture_match",
          kind="relationship_judgment", purpose=None, query="Q1", target="T1"):
    """Run the gates over one synthetic judgment. Returns (passed, checks)."""
    record = {"record_id": f"probe-{stratum}", "record_kind": kind, "stratum": stratum,
              "split": "curator_frozen", "_workflow": "sf_csa",
              "_gate_semantics": GATE_SEMANTICS}
    if purpose:
        record["control_purpose"] = purpose
    inv = {"query": query, "target": target}
    expected = {"judgment": {"function_classification": label, "structural_category": category},
                "invocation": inv}
    manifest = R.read_run_manifest(out, inv, "sf_csa")
    ev = json.loads((out / "SF_CSA_RELEASE_MANIFEST.json").read_text(encoding="utf-8"))
    checks: list[dict] = []
    R._gates_sf_csa(record, expected, ev, checks, manifest)
    failed = [c["check"] for c in checks if c.get("required") and not c["passed"]]
    return (not failed), failed


def main() -> int:
    root = Path(tempfile.mkdtemp(prefix="sf-csa-falsify-"))
    results: list[tuple[str, bool, str]] = []

    def record(name, ok, detail=""):
        results.append((name, ok, detail))

    try:
        # --- baseline: a well-formed homolog judgment must pass ------------
        out = write_release(root)
        ok, failed = judge(out, stratum="homologous_superfamily", label="same_mechanism_class")
        record("baseline homolog judgment passes", ok, f"unexpected failures: {failed}")

        # --- gate: the absolute bound on functional promotion --------------
        # An analogy pair promoted to a functional claim must fail.
        shutil.rmtree(root / "release")
        out = write_release(root, label="probable_same_function")
        ok, failed = judge(out, stratum="fold_analogy", label="probable_same_function")
        record("promoted fold analogy is rejected",
               (not ok) and "not_promoted_to_functional_claim" in failed, str(failed))

        shutil.rmtree(root / "release")
        out = write_release(root, label="structural_analogy_only")
        ok, failed = judge(out, stratum="fold_analogy", label="structural_analogy_only")
        record("unpromoted fold analogy passes", ok, str(failed))

        # --- gate: homolog recall ------------------------------------------
        # This is the failure Finding 2 predicts when the tables are not
        # overridden: a true homolog classified as mere analogy.
        shutil.rmtree(root / "release")
        out = write_release(root, label="structural_analogy_only")
        ok, failed = judge(out, stratum="homologous_superfamily", label="structural_analogy_only")
        record("unrecovered homolog is rejected",
               (not ok) and "homologous_superfamily_control_recovered" in failed, str(failed))

        # --- gate: exact controls ------------------------------------------
        shutil.rmtree(root / "release")
        out = write_release(root, label="same_mechanism_class")
        ok, failed = judge(out, stratum="exact", label="same_mechanism_class")
        record("exact control not recovered is rejected",
               (not ok) and "exact_control_recovered" in failed, str(failed))

        # --- gate: closed vocabulary ---------------------------------------
        shutil.rmtree(root / "release")
        out = write_release(root, label="probably_fine_honestly")
        ok, failed = judge(out, stratum="fold_analogy", label="probably_fine_honestly")
        record("label outside the closed vocabulary is rejected",
               (not ok) and "label_within_closed_vocabulary" in failed, str(failed))

        # --- gate: the judgment must not drift ------------------------------
        shutil.rmtree(root / "release")
        out = write_release(root, label="unresolved_or_conflicted")
        ok, failed = judge(out, stratum="homologous_superfamily", label="same_mechanism_class")
        record("changed classification is rejected",
               (not ok) and "function_classification_unchanged" in failed, str(failed))

        # --- gate: table provenance ----------------------------------------
        # Finding 2: which mechanism table was used decides whether the panel
        # measures anything. A release built against a different manifest fails.
        shutil.rmtree(root / "release")
        out = write_release(root, db_sha="0" * 64)
        ok, failed = judge(out, stratum="homologous_superfamily", label="same_mechanism_class")
        record("wrong database manifest is rejected",
               (not ok) and "database_manifest_matches_frozen_tables" in failed, str(failed))

        # An unfrozen expectation must fail closed, not skip the check.
        shutil.rmtree(root / "release")
        out = write_release(root)
        saved = GATE_SEMANTICS["sf_csa"].pop("database_manifest_sha256")
        ok, failed = judge(out, stratum="homologous_superfamily", label="same_mechanism_class")
        GATE_SEMANTICS["sf_csa"]["database_manifest_sha256"] = saved
        record("unfrozen database manifest fails closed",
               (not ok) and "database_manifest_matches_frozen_tables" in failed, str(failed))

        # --- gate: the pair must actually be in the release ------------------
        shutil.rmtree(root / "release")
        out = write_release(root, structural=False)
        ok, failed = judge(out, stratum="homologous_superfamily", label="same_mechanism_class")
        record("absent structural row is rejected",
               (not ok) and "judged_pair_present_in_release" in failed, str(failed))

        # --- derived missing-evidence really is derived ---------------------
        shutil.rmtree(root / "release")
        out = write_release(root, sequence=False)
        manifest = R.read_run_manifest(out, {"query": "Q1", "target": "T1"}, "sf_csa")
        record("absent sequence leg is reported missing",
               manifest["missing_evidence"] == ["sequence_comparison"],
               str(manifest["missing_evidence"]))

        shutil.rmtree(root / "release")
        out = write_release(root, orthology="")
        manifest = R.read_run_manifest(out, {"query": "Q1", "target": "T1"}, "sf_csa")
        record("absent orthology call is reported missing",
               manifest["missing_evidence"] == ["reciprocal_best_hit"],
               str(manifest["missing_evidence"]))

        # --- controls -------------------------------------------------------
        # The computed-path control passes only while the defect is present.
        shutil.rmtree(root / "release")
        out = write_release(root, label="same_mechanism_class", orthology="reciprocal_best_hit")
        ok, failed = judge(out, stratum="homologous_superfamily", label="same_mechanism_class",
                           kind="control_case", purpose="rbh_computed_path")
        record("rbh_computed_path control passes while defect exists", ok, str(failed))

        # ...and fails the day the module is fixed.
        shutil.rmtree(root / "release")
        out = write_release(root, label="probable_same_function", orthology="reciprocal_best_hit")
        ok, failed = judge(out, stratum="homologous_superfamily", label="probable_same_function",
                           kind="control_case", purpose="rbh_computed_path")
        record("rbh_computed_path control fails once the label is reachable",
               (not ok) and "computed_rbh_does_not_reach_function_label" in failed, str(failed))

        # The asserted-rbh control fires only on promotion without computed support.
        shutil.rmtree(root / "release")
        out = write_release(root, label="probable_same_function",
                            orthology="best_hit_nonreciprocal")
        ok, failed = judge(out, stratum="homologous_superfamily", label="probable_same_function",
                           kind="control_case", purpose="rbh_asserted_rejected")
        record("rbh_asserted_rejected control detects unsupported promotion", ok, str(failed))

        shutil.rmtree(root / "release")
        out = write_release(root, label="same_mechanism_class",
                            orthology="best_hit_nonreciprocal")
        ok, failed = judge(out, stratum="homologous_superfamily", label="same_mechanism_class",
                           kind="control_case", purpose="rbh_asserted_rejected")
        record("rbh_asserted_rejected control fails when nothing was promoted",
               (not ok) and "asserted_rbh_detected_as_unsupported" in failed, str(failed))

        # An unknown control purpose must not pass silently.
        shutil.rmtree(root / "release")
        out = write_release(root)
        ok, failed = judge(out, stratum="homologous_superfamily", label="same_mechanism_class",
                           kind="control_case", purpose="something_invented")
        record("unknown control purpose is rejected",
               (not ok) and "control_purpose_declared" in failed, str(failed))
    finally:
        shutil.rmtree(root, ignore_errors=True)

    width = max(len(n) for n, _, _ in results)
    for name, ok, detail in results:
        print(f"{'ok  ' if ok else 'FAIL'}  {name.ljust(width)}  {'' if ok else detail}")
    bad = [n for n, ok, _ in results if not ok]
    print(f"\n{len(results) - len(bad)}/{len(results)} gate probes behaved as declared")
    if bad:
        print("unfalsifiable or misbehaving gates: " + ", ".join(bad), file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
