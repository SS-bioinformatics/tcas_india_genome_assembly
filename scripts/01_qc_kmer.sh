#!/usr/bin/env bash
# ==============================================================================
# STEP 01 — Read QC + k-mer genome-size estimation
#
# Thesis (Appendix I p.67):
#   FastQC v0.11.8  →  Trimmomatic v0.38 (LEADING:20 TRAILING:20
#   SLIDINGWINDOW:4:20 MINLEN:50)  →  Jellyfish v2 (k=21, canonical)  →
#   GenomeScope 2.0 (ploidy=2)  →  validated with KMC.
#   Result: estimated genome size ~165 Mbp, heterozygosity ~1.2%.
#
# Conda env: env_qc  (fastqc=0.11.8, trimmomatic=0.38, kmer-jellyfish, kmc,
#                     genomescope2 — installed via work/00_environment/)
#
# Usage:
#   bash 01_qc_kmer.sh           # check: list existing outputs
#   bash 01_qc_kmer.sh run       # run per-sample QC (skips if already done)
#   bash 01_qc_kmer.sh run --force  # force re-run even if done
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/01_qc_kmer"

banner "STEP 01: Read QC + k-mer genome-size estimation [$RUN_MODE]"
require_paths "$RAW_ILLUMINA_WGS" "$QC_KMER_GENOMESCOPE"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing GenomeScope / Jellyfish outputs (D:):"
  find "$QC_KMER_GENOMESCOPE" -maxdepth 2 2>/dev/null | sort | sed 's/^/  /'
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

activate_conda env_qc

# Discover samples: each subdirectory under RAW_ILLUMINA_WGS is one individual
mapfile -t SAMPLES < <(find "$RAW_ILLUMINA_WGS" -mindepth 1 -maxdepth 1 -type d | sort)
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "[ERROR] No sample directories found under $RAW_ILLUMINA_WGS" >&2
  exit 1
fi
echo "Found ${#SAMPLES[@]} sample directories."

for SDIR in "${SAMPLES[@]}"; do
  SNAME="$(basename "$SDIR")"
  SOUT="$WORKDIR/$SNAME"
  mkdir -p "$SOUT"

  # Locate R1/R2 (accept .fastq.gz or .fq.gz)
  R1=$(find "$SDIR" -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" 2>/dev/null | sort | head -1)
  R2=$(find "$SDIR" -name "*_R2*.fastq.gz" -o -name "*_R2*.fq.gz" 2>/dev/null | sort | head -1)
  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "  [SKIP] $SNAME — R1/R2 not found, skipping"
    continue
  fi

  echo
  echo "=== Processing sample: $SNAME ==="

  # --- FastQC (raw) ----------------------------------------------------------
  echo "  FastQC (raw)..."
  fastqc -t 8 -o "$SOUT" "$R1" "$R2"

  # --- Trimmomatic -----------------------------------------------------------
  echo "  Trimmomatic..."
  R1T="$SOUT/${SNAME}_R1.trimmed.fastq.gz"
  R2T="$SOUT/${SNAME}_R2.trimmed.fastq.gz"
  trimmomatic PE -threads 8 -phred33 \
    "$R1" "$R2" \
    "$R1T" "$SOUT/${SNAME}_R1.unpaired.fastq.gz" \
    "$R2T" "$SOUT/${SNAME}_R2.unpaired.fastq.gz" \
    ILLUMINACLIP:/root/miniforge3/envs/env_qc/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

  # --- FastQC (trimmed) -------------------------------------------------------
  echo "  FastQC (trimmed)..."
  fastqc -t 8 -o "$SOUT" "$R1T" "$R2T"

  # --- Jellyfish k=21 ---------------------------------------------------------
  echo "  Jellyfish count (k=21)..."
  jellyfish count -C -m 21 -s 1G -t 8 \
    -o "$SOUT/${SNAME}_k21.jf" \
    <(zcat "$R1T") <(zcat "$R2T")
  jellyfish histo -t 8 "$SOUT/${SNAME}_k21.jf" > "$SOUT/${SNAME}_k21.hist"

  # --- GenomeScope 2.0 --------------------------------------------------------
  echo "  GenomeScope 2.0..."
  mkdir -p "$SOUT/genomescope"
  genomescope2 \
    -i "$SOUT/${SNAME}_k21.hist" \
    -o "$SOUT/genomescope" \
    -k 21 -p 2 \
    --name_prefix "$SNAME"

  echo "  [OK] $SNAME done"
done

mark_done "$WORKDIR"
