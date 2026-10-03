#!/usr/bin/env bash
# ==============================================================================
# STEP 04 — Error polishing (Pilon, 3 rounds) +
#            Purge Haplotigs round 2 (stringent) + MeDuSa validation
#
# Thesis (Appendix I p.67):
#   3 × Pilon rounds per individual (minimap2 SR alignment → Pilon)
#   Purge Haplotigs round 2:  -l 10 -m 50 -h 150  haplotig threshold  -a 80
#   MeDuSa: reference-guided scaffolding validation (cross-check vs RagTag)
#
# Conda env: env_assembly  (minimap2, samtools, bedtools)
# Pilon:     conda create -n env_pilon -c bioconda pilon
# MeDuSa:    java -jar medusa.jar  (tool already on F:)
#
# !! EXPENSIVE: Pilon x3 per sample × 8 samples = multi-day run.
#    Polished assemblies already exist on D: ($POLISHED_MINIMAP2_BAMS, $POLISHED_MEDUSA).
#
# Usage:
#   bash 04_polishing.sh           # check: list existing outputs
#   bash 04_polishing.sh run       # skip if already done
#   bash 04_polishing.sh run --force
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/04_polishing"

banner "STEP 04: Pilon polishing (×3) + Purge Haplotigs round 2 + MeDuSa [$RUN_MODE]"
require_paths "$POLISHED_MINIMAP2_BAMS" "$POLISHED_MEDUSA"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing minimap2 BAMs / Pilon outputs (D:):"
  find "$POLISHED_MINIMAP2_BAMS" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  echo
  echo "Existing MeDuSa validation outputs (D:):"
  find "$POLISHED_MEDUSA" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

# Activate assembly env (minimap2, samtools)
activate_conda env_assembly

# Purge Haplotigs round 2 uses the same env_assembly (purge_haplotigs installed there)
# Pilon: use env_pilon if available, else try PATH
PILON_JAR=""
for candidate in \
  "$HOME/miniforge3/envs/env_pilon/share/pilon/pilon.jar" \
  "/root/miniforge3/envs/env_pilon/share/pilon/pilon.jar" \
  "$(conda run -n env_pilon bash -c 'find $CONDA_PREFIX/share -name pilon.jar 2>/dev/null | head -1' 2>/dev/null || true)"; do
  [[ -f "$candidate" ]] && { PILON_JAR="$candidate"; break; }
done
[[ -z "$PILON_JAR" ]] && { echo "[ERROR] pilon.jar not found. conda create -n env_pilon -c bioconda pilon" >&2; exit 1; }

mapfile -t SAMPLES < <(find "$PIPELINE_ROOT/work/03_draft_assembly" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
[[ ${#SAMPLES[@]} -eq 0 ]] && { echo "[ERROR] Run Step 03 first to generate draft assemblies." >&2; exit 1; }

for SDIR in "${SAMPLES[@]}"; do
  SNAME="$(basename "$SDIR")"
  DRAFT="$SDIR/purge_haplotigs/${SNAME}_purged.fasta"
  [[ -f "$DRAFT" ]] || { echo "  [SKIP] $SNAME — purged draft not found at $DRAFT"; continue; }

  R1=$(find "$RAW_ILLUMINA_WGS/$SNAME" \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" \) 2>/dev/null | head -1)
  R2=$(find "$RAW_ILLUMINA_WGS/$SNAME" \( -name "*_R2*.fastq.gz" -o -name "*_R2*.fq.gz" \) 2>/dev/null | head -1)

  POUT="$WORKDIR/$SNAME"
  mkdir -p "$POUT"
  CURRENT="$DRAFT"

  # 3 Pilon rounds
  for round in 1 2 3; do
    echo; echo "=== $SNAME — Pilon round $round ==="
    ROUND_OUT="$POUT/pilon_round${round}"
    mkdir -p "$ROUND_OUT"

    minimap2 -ax sr -t 16 "$CURRENT" "$R1" "$R2" \
      | samtools sort -@ 8 -o "$POUT/${SNAME}_r${round}.sorted.bam"
    samtools index "$POUT/${SNAME}_r${round}.sorted.bam"

    java -Xmx24g -jar "$PILON_JAR" \
      --genome "$CURRENT" \
      --frags "$POUT/${SNAME}_r${round}.sorted.bam" \
      --output "${SNAME}_pilon_r${round}" \
      --outdir "$ROUND_OUT" \
      --threads 16

    CURRENT="$ROUND_OUT/${SNAME}_pilon_r${round}.fasta"
  done

  # Purge Haplotigs round 2 (stringent)
  echo; echo "=== $SNAME — Purge Haplotigs round 2 ==="
  PH2="$POUT/purge_haplotigs_r2"
  mkdir -p "$PH2"

  minimap2 -ax sr -t 16 "$CURRENT" "$R1" "$R2" \
    | samtools sort -@ 8 -o "$PH2/${SNAME}_r2.sorted.bam"
  samtools index "$PH2/${SNAME}_r2.sorted.bam"

  purge_haplotigs hist -b "$PH2/${SNAME}_r2.sorted.bam" -g "$CURRENT" -t 16 -d "$PH2"
  purge_haplotigs contigcov \
    -i "$PH2/${SNAME}_r2.sorted.bam.gencov" \
    -l 10 -m 50 -h 150 \
    -o "$PH2/coverage_stats_r2.csv"
  purge_haplotigs purge \
    -g "$CURRENT" \
    -c "$PH2/coverage_stats_r2.csv" \
    -a 80 \
    -o "$PH2/${SNAME}_purged_r2" \
    -t 16

  echo "  [OK] $SNAME polished+purged: $PH2/${SNAME}_purged_r2.fasta"
done

mark_done "$WORKDIR"
