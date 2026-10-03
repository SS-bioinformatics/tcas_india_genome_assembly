#!/usr/bin/env bash
# ==============================================================================
# STEP 00 — Raw data inventory
#
# Thesis (Ch.2 p.18, Appendix I p.67):
#   • Illumina 2×150 bp PE WGS for 8 individuals (Table S1)
#   • Oxford Nanopore MinION long reads from 1 hyper-inbred female (Table S2)
#   • 2 RNA-seq datasets (Imroz_mRNA_27122017) used as MAKER EST evidence:
#       dataset 1 — 12-day-old outbred females (1.7 M reads)
#       dataset 2 — pooled males/females/eggs from wild population (4.1 M reads)
#
# This step performs no tool runs — it only verifies that raw inputs exist and
# prints a catalogue of what is present before anything downstream reads them.
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

banner "STEP 00: Raw data inventory"
echo
echo "Checking canonical raw data directories..."
require_paths "$RAW_ILLUMINA_WGS" "$RAW_NANOPORE" "$RAW_RNASEQ_IMROZ"

echo
echo "── Illumina WGS samples (D:, Table S1) ──────────────────────────"
find "$RAW_ILLUMINA_WGS" -maxdepth 2 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) 2>/dev/null \
  | sort | sed 's/^/  /'

echo
echo "── Nanopore raw reads (F:, Table S2) ────────────────────────────"
find "$RAW_NANOPORE" -maxdepth 2 \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fast5" \) 2>/dev/null \
  | sort | head -30 | sed 's/^/  /'

echo
echo "── RNA-seq (Imroz_mRNA_27122017, E:) ────────────────────────────"
find "$RAW_RNASEQ_IMROZ" -maxdepth 2 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) 2>/dev/null \
  | sort | sed 's/^/  /'

echo
echo "[Step 00 complete — data-only check, no tool run needed]"
