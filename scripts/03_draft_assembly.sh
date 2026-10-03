#!/usr/bin/env bash
# ==============================================================================
# STEP 03 — Draft assembly per individual (SPAdes v3.15.2) +
#            haplotig removal (Purge Haplotigs v1.0.0, round 1)
#
# Thesis (Appendix I p.67):
#   SPAdes v3.15.2 per individual (8 samples):  -t 16 -m 128
#   Purge Haplotigs round 1 — per-sample coverage cutoffs:
#     Individual 1A: -l 5 -m 30 -h 190 -a 40
#     Others: cutoffs derived from GenomeScope histogram (see $QC_KMER_GENOMESCOPE)
#   Coverage histogram (hist step) uses minimap2 alignments of reads to draft.
#
# Conda env: env_assembly  (spades=3.15.2, purge_haplotigs=1.0.1, minimap2,
#                           samtools, bedtools — work/00_environment/install_env_assembly.sh)
# Existing outputs on F:  $HYBRIDSPADES_OUTPUT, $PURGE_HAPLOTIGS_OUTPUT
#
# !! EXPENSIVE: SPAdes per sample = hours; 8 samples = ~1–2 days.
#    Only re-run if you need to reproduce from trimmed reads.
#
# Usage:
#   bash 03_draft_assembly.sh           # check: list existing outputs
#   bash 03_draft_assembly.sh run       # skip if already done
#   bash 03_draft_assembly.sh run --force  # force re-run ALL samples
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/03_draft_assembly"

# Per-sample Purge Haplotigs round-1 parameters (from GenomeScope peak analysis)
# Format: SAMPLE_NAME -> "low mid high haplotig_threshold"
declare -A PH_PARAMS=(
  [da01a]="5 30 190 40"   # individual 1A — from thesis text
  [da01b]="5 25 160 40"
  [da01n]="5 25 160 40"
  [12B]="5 25 160 40"
  [12D]="5 25 160 40"
  [12E]="5 25 160 40"
  [13C]="5 25 160 40"
  [13E]="5 25 160 40"
  [18A]="5 25 160 40"
  [18B]="5 25 160 40"
)

banner "STEP 03: SPAdes draft assembly + Purge Haplotigs round 1 [$RUN_MODE]"
require_paths "$RAW_ILLUMINA_WGS" "$HYBRIDSPADES_OUTPUT" "$PURGE_HAPLOTIGS_OUTPUT"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing SPAdes assemblies (F:):"
  find "$HYBRIDSPADES_OUTPUT" -maxdepth 1 -type d 2>/dev/null | sort | sed 's/^/  /'
  echo
  echo "Existing Purge Haplotigs outputs (F:):"
  find "$PURGE_HAPLOTIGS_OUTPUT" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

activate_conda env_assembly

mapfile -t SAMPLES < <(find "$RAW_ILLUMINA_WGS" -mindepth 1 -maxdepth 1 -type d | sort)

for SDIR in "${SAMPLES[@]}"; do
  SNAME="$(basename "$SDIR")"
  SOUT="$WORKDIR/$SNAME"
  mkdir -p "$SOUT"

  R1=$(find "$SDIR" \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" \) 2>/dev/null | sort | head -1)
  R2=$(find "$SDIR" \( -name "*_R2*.fastq.gz" -o -name "*_R2*.fq.gz" \) 2>/dev/null | sort | head -1)
  # Fall back to trimmed reads in work/01 if available
  [[ -z "$R1" ]] && R1=$(find "$PIPELINE_ROOT/work/01_qc_kmer/$SNAME" -name "*_R1.trimmed.fastq.gz" 2>/dev/null | head -1)
  [[ -z "$R2" ]] && R2=$(find "$PIPELINE_ROOT/work/01_qc_kmer/$SNAME" -name "*_R2.trimmed.fastq.gz" 2>/dev/null | head -1)

  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "  [SKIP] $SNAME — R1/R2 not found"; continue
  fi

  echo; echo "=== SPAdes: $SNAME ==="
  SPADES_OUT="$SOUT/spades"
  spades.py \
    --pe1-1 "$R1" --pe1-2 "$R2" \
    -o "$SPADES_OUT" \
    -t 16 -m 128

  DRAFT="$SPADES_OUT/scaffolds.fasta"
  [[ -f "$DRAFT" ]] || { echo "[ERROR] SPAdes did not produce scaffolds.fasta" >&2; exit 1; }

  echo; echo "=== Purge Haplotigs round 1: $SNAME ==="
  PH_OUT="$SOUT/purge_haplotigs"
  mkdir -p "$PH_OUT"

  # Align reads back to draft for coverage histogram
  minimap2 -ax sr -t 16 "$DRAFT" "$R1" "$R2" \
    | samtools sort -@ 8 -o "$PH_OUT/${SNAME}.sorted.bam"
  samtools index "$PH_OUT/${SNAME}.sorted.bam"

  # Coverage histogram
  purge_haplotigs hist \
    -b "$PH_OUT/${SNAME}.sorted.bam" \
    -g "$DRAFT" -t 16 \
    -d "$PH_OUT"

  # Get per-sample params (default fallback if not defined)
  PARAMS="${PH_PARAMS[$SNAME]:-5 30 190 40}"
  read -r L M H A <<< "$PARAMS"

  purge_haplotigs contigcov \
    -i "$PH_OUT/${SNAME}.sorted.bam.gencov" \
    -l "$L" -m "$M" -h "$H" \
    -o "$PH_OUT/coverage_stats.csv"

  purge_haplotigs purge \
    -g "$DRAFT" \
    -c "$PH_OUT/coverage_stats.csv" \
    -a "$A" \
    -o "$PH_OUT/${SNAME}_purged" \
    -t 16

  echo "  [OK] $SNAME purged assembly: $PH_OUT/${SNAME}_purged.fasta"
done

mark_done "$WORKDIR"
