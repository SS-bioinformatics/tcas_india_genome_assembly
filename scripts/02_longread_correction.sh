#!/usr/bin/env bash
# ==============================================================================
# STEP 02 — Nanopore long-read correction (Canu)
#
# Thesis (Appendix I p.67-68): ONT reads from the hyper-inbred female were
# corrected with Canu before use as gap-filling donor reads in Step 06.
# Parameters: genomeSize=165m, -nanopore, useGrid=false, maxMemory=28.
#
# Conda env: env_canu  (conda create -n env_canu -c bioconda canu)
# Corrected reads already exist at: $NANOPORE_CANU_CORRECTED
#
# !! EXPENSIVE: takes several hours on raw ONT data.
#    Only re-run if you need to reproduce from raw reads.
#
# Usage:
#   bash 02_longread_correction.sh           # check: list existing outputs
#   bash 02_longread_correction.sh run       # skip if already done
#   bash 02_longread_correction.sh run --force  # force re-run
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/02_longread_correction"

banner "STEP 02: Nanopore long-read correction (Canu) [$RUN_MODE]"
require_paths "$RAW_NANOPORE" "$NANOPORE_CANU_CORRECTED"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing Canu-corrected reads (F:):"
  find "$NANOPORE_CANU_CORRECTED" -maxdepth 2 2>/dev/null | sort | sed 's/^/  /'
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

activate_conda env_canu

mapfile -t ONT_READS < <(find "$RAW_NANOPORE" \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)
if [[ ${#ONT_READS[@]} -eq 0 ]]; then
  echo "[ERROR] No .fastq.gz reads found under $RAW_NANOPORE" >&2; exit 1
fi
echo "Found ${#ONT_READS[@]} Nanopore read files."

canu -correct \
  -p inTcas1_ont \
  -d "$WORKDIR/canu_correction" \
  genomeSize=165m \
  maxThreads=16 \
  maxMemory=28 \
  useGrid=false \
  -nanopore "${ONT_READS[@]}"

echo
echo "Corrected reads: $WORKDIR/canu_correction/inTcas1_ont.correctedReads.fasta.gz"
mark_done "$WORKDIR"
