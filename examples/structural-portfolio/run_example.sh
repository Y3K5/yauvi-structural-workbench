#!/bin/sh
set -eu

EXAMPLE_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [ "$#" -ne 1 ]; then
  echo "usage: $0 NEW_OUTPUT_DIRECTORY" >&2
  echo "Public inputs are written beside it as <output>.inputs." >&2
  exit 2
fi

OUTPUT_DIR=$1
INPUT_DIR="${OUTPUT_DIR}.inputs"
PYTHON=${PYTHON:-python3}

"$PYTHON" "$EXAMPLE_DIR/acquire_public_inputs.py" --out "$INPUT_DIR"
"$PYTHON" "$EXAMPLE_DIR/public_research_case.py" \
  --inputs "$INPUT_DIR" \
  --out "$OUTPUT_DIR"
echo "Public CA II case written to $OUTPUT_DIR"
