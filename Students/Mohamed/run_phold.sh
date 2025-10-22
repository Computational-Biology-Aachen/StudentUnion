#!/usr/bin/env bash
# --------------------------------------------------------------
# PHOLD pipeline wrapper
# Usage: ./run_phold.sh -i <input_gbk> [-o outdir] [-t threads]
# Example: ./run_phold.sh -i output_pharokka/pharokka.gbk -o output_phold -t 8
# --------------------------------------------------------------

set -euo pipefail

# Default values
THREADS=8
OUTDIR="output_phold"

# Parse command-line options
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      INPUT_GBK="$2"
      shift 2
      ;;
    -o|--outdir)
      OUTDIR="$2"
      shift 2
      ;;
    -t|--threads)
      THREADS="$2"
      shift 2
      ;;
    *)
      echo "❌ Unknown option: $1"
      echo "Usage: $0 -i <input_gbk> [-o outdir] [-t threads]"
      exit 1
      ;;
  esac
done

# Validate required argument
if [[ -z "${INPUT_GBK:-}" ]]; then
  echo "❌ Error: Input GenBank file (-i) is required."
  exit 1
fi

# Print summary
echo "--------------------------------------------------------------"
echo "📄 Input GenBank:    $INPUT_GBK"
echo "📦 Output directory: $OUTDIR"
echo "🧠 Threads:          $THREADS"
echo "--------------------------------------------------------------"

# Run PHOLD
phold run \
  -i "$INPUT_GBK" \
  -o "$OUTDIR" \
  -t "$THREADS"

echo "✅ PHOLD completed successfully."
echo "🗃 Results saved in: $OUTDIR"
