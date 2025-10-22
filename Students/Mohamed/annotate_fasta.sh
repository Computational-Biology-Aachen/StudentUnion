#!/usr/bin/env bash
# --------------------------------------------------------------
# Pharokka annotation wrapper
# Usage: ./annotate_fasta.sh <input_fasta> [threads] [database_dir]
# Example: ./annotate_fasta.sh test2.fasta 8 /home/mohamed/pharokka_db
# --------------------------------------------------------------

# Exit on errors
set -euo pipefail

# Input FASTA file (required)
INPUT_FASTA="${1:-}"
if [[ -z "$INPUT_FASTA" ]]; then
    echo "❌ Error: No input FASTA file specified."
    echo "Usage: $0 <input_fasta> [threads] [database_dir]"
    exit 1
fi

# Number of threads (default 8)
THREADS="${2:-8}"

# Database directory (default ./pharokka)
DB_DIR="${3:-./pharokka}"

# Output directory (based on input file name)
BASENAME=$(basename "$INPUT_FASTA" .fasta)
OUTPUT_DIR="output_pharokka"

echo "--------------------------------------------------------------"
echo "📁 Input FASTA:      $INPUT_FASTA"
echo "🧠 Threads:          $THREADS"
echo "🗂️ Database dir:     $DB_DIR"
echo "📦 Output directory: $OUTPUT_DIR"
echo "--------------------------------------------------------------"

# Run Pharokka
pharokka.py \
    -i "$INPUT_FASTA" \
    -o "$OUTPUT_DIR" \
    -d "$DB_DIR" \
    -t "$THREADS" \
    --fast \
    -f

echo "✅ Pharokka annotation completed successfully."
echo "🗃 Results saved in: $OUTPUT_DIR"
