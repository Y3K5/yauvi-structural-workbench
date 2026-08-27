#!/bin/sh
set -u

EXAMPLE_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PORTFOLIO_ROOT=$(CDPATH= cd -- "$EXAMPLE_DIR/../.." && pwd)
OUTPUT_DIR=${1:-"$EXAMPLE_DIR/output"}
PORTFOLIO_PYTHONPATH="$PORTFOLIO_ROOT/structqc/src:$PORTFOLIO_ROOT/assembly-context/src:$PORTFOLIO_ROOT/site-context/src:$PORTFOLIO_ROOT/state-atlas/src:$PORTFOLIO_ROOT/structcons/src:$PORTFOLIO_ROOT/structrel/src:$PORTFOLIO_ROOT/dockprep/src:$PORTFOLIO_ROOT/pose-evidence/src:$PORTFOLIO_ROOT/structdesign/src:$PORTFOLIO_ROOT/oral-context/src"
export PYTHONPATH="$PORTFOLIO_PYTHONPATH${PYTHONPATH:+:$PYTHONPATH}"

mkdir -p "$OUTPUT_DIR/structqc" "$OUTPUT_DIR/assembly" "$OUTPUT_DIR/site" "$OUTPUT_DIR/state" "$OUTPUT_DIR/conservation" "$OUTPUT_DIR/relationships" "$OUTPUT_DIR/docking" "$OUTPUT_DIR/poses" "$OUTPUT_DIR/design" "$OUTPUT_DIR/oral"

run_analysis() {
  "$@"
  analysis_status=$?
  if [ "$analysis_status" -eq 0 ] || [ "$analysis_status" -eq 1 ]; then
    return 0
  fi
  return "$analysis_status"
}

run_analysis python3 -m structqc.cli run \
  --structure "$EXAMPLE_DIR/query.pdb" \
  --subject-id SYNTHETIC_QUERY \
  --reference-fasta "$EXAMPLE_DIR/reference.fasta" \
  --provenance "$EXAMPLE_DIR/provenance.json" \
  --chain A --out "$OUTPUT_DIR/structqc" || exit $?

run_analysis python3 -m assembly_context.cli run \
  --manifest "$OUTPUT_DIR/structqc/STRUCTURE_EVIDENCE.json" \
  --isolated "$EXAMPLE_DIR/query.pdb" \
  --assembly "$EXAMPLE_DIR/assembly.pdb" \
  --subject-chain A --relationship exact_protein --assembly-id SYNTHETIC_AB \
  --expected-chains A,B --out "$OUTPUT_DIR/assembly" || exit $?

run_analysis python3 -m site_context.cli run \
  --manifest "$OUTPUT_DIR/structqc/STRUCTURE_EVIDENCE.json" \
  --structure "$EXAMPLE_DIR/query.pdb" \
  --annotations "$EXAMPLE_DIR/annotations.json" \
  --pocket-result "$EXAMPLE_DIR/pockets.json" --out "$OUTPUT_DIR/site" || exit $?

run_analysis python3 -m state_atlas.cli run \
  --manifest "$OUTPUT_DIR/structqc/STRUCTURE_EVIDENCE.json" \
  --reference-set "$EXAMPLE_DIR/reference_set.json" \
  --structure "$EXAMPLE_DIR/query.pdb" --chain A \
  --cluster-cutoff-A 0.5 --out "$OUTPUT_DIR/state" || exit $?

run_analysis python3 -m structcons.cli run \
  --manifest "$OUTPUT_DIR/structqc/STRUCTURE_EVIDENCE.json" \
  --msa "$EXAMPLE_DIR/alignment.fasta" --query-id SYNTHETIC_QUERY \
  --variants "$EXAMPLE_DIR/variants.tsv" \
  --layer "$OUTPUT_DIR/assembly/ASSEMBLY_LAYER.json" \
  --layer "$OUTPUT_DIR/site/SITE_LAYER.json" \
  --spatial-threshold-A 8 --out "$OUTPUT_DIR/conservation" || exit $?

run_analysis python3 -m structrel.cli run \
  --manifest "$OUTPUT_DIR/structqc/STRUCTURE_EVIDENCE.json" \
  --relationships "$EXAMPLE_DIR/relationships.json" \
  --equivalences "$EXAMPLE_DIR/equivalences.json" \
  --out "$OUTPUT_DIR/relationships" || exit $?

run_analysis python3 -m dockprep.cli run \
  --config "$EXAMPLE_DIR/docking_config.json" \
  --out "$OUTPUT_DIR/docking" || exit $?

run_analysis python3 -m pose_evidence.cli run \
  --config "$EXAMPLE_DIR/pose_config.json" \
  --out "$OUTPUT_DIR/poses" || exit $?

run_analysis python3 -m structdesign.cli run \
  --config "$EXAMPLE_DIR/design_config.json" \
  --out "$OUTPUT_DIR/design" || exit $?

run_analysis python3 -m oral_context.cli run \
  --config "$EXAMPLE_DIR/oral_context_config.json" \
  --out "$OUTPUT_DIR/oral" || exit $?

python3 "$PORTFOLIO_ROOT/tools/build_structural_report.py" \
  --structure "$EXAMPLE_DIR/query.pdb" \
  --layer "$OUTPUT_DIR/structqc/STRUCTURE_LAYER.json" \
  --layer "$OUTPUT_DIR/assembly/ASSEMBLY_LAYER.json" \
  --layer "$OUTPUT_DIR/site/SITE_LAYER.json" \
  --layer "$OUTPUT_DIR/state/STATE_LAYER.json" \
  --layer "$OUTPUT_DIR/conservation/CONSERVATION_LAYER.json" \
  --layer "$OUTPUT_DIR/relationships/RELATIONSHIP_LAYER.json" \
  --layer "$OUTPUT_DIR/poses/POSE_LAYER.json" \
  --layer "$OUTPUT_DIR/design/DESIGN_LAYER.json" \
  --layer "$OUTPUT_DIR/oral/ORAL_CONTEXT_LAYER.json" \
  --summary "$OUTPUT_DIR/structqc/STRUCTURE_EVIDENCE.json" \
  --summary "$OUTPUT_DIR/assembly/ASSEMBLY_CONTEXT.json" \
  --summary "$OUTPUT_DIR/site/SITE_CONTEXT.json" \
  --summary "$OUTPUT_DIR/state/STATE_ENSEMBLE.json" \
  --summary "$OUTPUT_DIR/conservation/STRUCTURAL_CONSERVATION.json" \
  --summary "$OUTPUT_DIR/relationships/STRUCTURAL_RELATIONSHIPS.json" \
  --summary "$OUTPUT_DIR/docking/DOCKING_JOB.json" \
  --summary "$OUTPUT_DIR/poses/POSE_EVIDENCE.json" \
  --summary "$OUTPUT_DIR/design/DESIGN_CANDIDATES.json" \
  --summary "$OUTPUT_DIR/oral/ORAL_CONTEXT.json" \
  --out "$OUTPUT_DIR/STRUCTURAL_REPORT.html" || exit $?

echo "Structural portfolio example written to $OUTPUT_DIR"
