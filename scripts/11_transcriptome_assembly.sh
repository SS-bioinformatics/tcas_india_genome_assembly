#!/usr/bin/env bash
# ==============================================================================
# STEP 11 — Transcriptome assembly (SOAPdenovo-Trans, MAKER EST evidence)
#
# Thesis (Appendix I p.68, Table S7):
#   2 RNA-seq datasets from $RAW_RNASEQ_IMROZ:
#     dataset 1: 12-day-old outbred females (1.7 M reads)
#     dataset 2: pooled males/females/eggs, wild population (4.1 M reads)
#   QC: FastQC → Trimmomatic (same params as Step 01: LEADING:20 TRAILING:20 MINLEN:50)
#   SOAPdenovo-Trans k-mer sweep:
#     k = 21, 23, 25, 27, 29, 31 (step 2) and k = 41, 51, 61, 71, 81, 91, 101, 111, 121 (step 10)
#   Best k: k=21 for dataset 1, k=31 for dataset 2 (chosen by total length + BUSCO score)
#   Output transcripts used as EST evidence (-est_gff) in MAKER round 1 (Step 13).
#
# Conda env: env_transcriptome
#   conda create -n env_transcriptome -c bioconda soapdenovo-trans busco
#   (FastQC/Trimmomatic already in env_qc)
#
# Existing outputs: $TRANSCRIPTOME_ASSEMBLY_OUTPUT, $TRANSCRIPTOME_JNS_RNASEQ (F:)
#
# Usage:
#   bash 11_transcriptome_assembly.sh           # check: list existing outputs
#   bash 11_transcriptome_assembly.sh run       # run k-mer sweep (both datasets)
#   bash 11_transcriptome_assembly.sh run --force
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/11_transcriptome_assembly"

# k-mer sweep ranges (step 2 and step 10 per thesis)
K_SMALL=(21 23 25 27 29 31)
K_LARGE=(41 51 61 71 81 91 101 111 121)
ALL_K=("${K_SMALL[@]}" "${K_LARGE[@]}")

banner "STEP 11: Transcriptome assembly (SOAPdenovo-Trans) [$RUN_MODE]"
require_paths "$RAW_RNASEQ_IMROZ" "$TRANSCRIPTOME_ASSEMBLY_OUTPUT"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing transcriptome assembly outputs (F:):"
  find "$TRANSCRIPTOME_ASSEMBLY_OUTPUT" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  echo
  echo "JNS RNA-seq assembly outputs (F:):"
  find "$TRANSCRIPTOME_JNS_RNASEQ" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

# QC: use env_qc for FastQC + Trimmomatic
activate_conda env_qc

# Discover RNA-seq files (dataset 1 and 2 under RAW_RNASEQ_IMROZ)
mapfile -t ALL_R1 < <(find "$RAW_RNASEQ_IMROZ" -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" | sort)
mapfile -t ALL_R2 < <(find "$RAW_RNASEQ_IMROZ" -name "*_R2*.fastq.gz" -o -name "*_R2*.fq.gz" | sort)

[[ ${#ALL_R1[@]} -eq 0 ]] && { echo "[ERROR] No RNA-seq R1 files found under $RAW_RNASEQ_IMROZ" >&2; exit 1; }

TRIM_DIR="$WORKDIR/trimmed"
mkdir -p "$TRIM_DIR"

for i in "${!ALL_R1[@]}"; do
  R1="${ALL_R1[$i]}"
  R2="${ALL_R2[$i]:-}"
  DSNAME="dataset$((i+1))"

  echo; echo "=== Trimmomatic: $DSNAME ==="
  if [[ -n "$R2" && -f "$R2" ]]; then
    trimmomatic PE -threads 8 -phred33 \
      "$R1" "$R2" \
      "$TRIM_DIR/${DSNAME}_R1.trimmed.fastq.gz" "$TRIM_DIR/${DSNAME}_R1.unpaired.fastq.gz" \
      "$TRIM_DIR/${DSNAME}_R2.trimmed.fastq.gz" "$TRIM_DIR/${DSNAME}_R2.unpaired.fastq.gz" \
      ILLUMINACLIP:/root/miniforge3/envs/env_qc/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
      LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
  else
    trimmomatic SE -threads 8 -phred33 \
      "$R1" "$TRIM_DIR/${DSNAME}_R1.trimmed.fastq.gz" \
      LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
  fi
done

# Now switch to env_transcriptome for SOAPdenovo-Trans
conda deactivate 2>/dev/null || true
activate_conda env_transcriptome

for i in "${!ALL_R1[@]}"; do
  DSNAME="dataset$((i+1))"
  DS_DIR="$WORKDIR/$DSNAME"
  mkdir -p "$DS_DIR"

  # Write SOAPdenovo-Trans config
  CONF="$DS_DIR/soap_config.txt"
  {
    echo "max_rd_len=150"
    echo "[LIB]"
    echo "avg_ins=300"
    echo "reverse_seq=0"
    echo "asm_flags=3"
    echo "rank=1"
    echo "q1=$TRIM_DIR/${DSNAME}_R1.trimmed.fastq.gz"
    [[ -f "$TRIM_DIR/${DSNAME}_R2.trimmed.fastq.gz" ]] && echo "q2=$TRIM_DIR/${DSNAME}_R2.trimmed.fastq.gz"
  } > "$CONF"

  for k in "${ALL_K[@]}"; do
    echo "  SOAPdenovo-Trans k=$k ($DSNAME)..."
    mkdir -p "$DS_DIR/k${k}"
    SOAPdenovo-Trans-127mer all \
      -s "$CONF" \
      -K "$k" \
      -o "$DS_DIR/k${k}/trans_k${k}" \
      -p 16
  done

  echo "  k-mer sweep complete for $DSNAME."
  echo "  Best k (thesis): dataset1=k21, dataset2=k31"
  echo "  BUSCO each assembly if needed: busco -i trans_k<k>.scafSeq -l insecta_odb10 -m transcriptome"
done

# Record the best assemblies
echo
echo "Best transcriptome assemblies for MAKER EST evidence (thesis Table S7):"
echo "  Dataset 1 (k=21): $WORKDIR/dataset1/k21/trans_k21.scafSeq"
echo "  Dataset 2 (k=31): $WORKDIR/dataset2/k31/trans_k31.scafSeq"

mark_done "$WORKDIR"
