#!/usr/bin/env bash
# ==============================================================================
# STEP 09 — Assembly QC: QUAST + BUSCO + nucmer/MUMmer structural comparison
#
# Thesis (Methods p.18-21, Table 1, Table S5, Figure 1):
#   QUAST  : basic stats vs Tcas5.2 reference
#   BUSCO  : insecta_odb10 (database obtained 2021-12-10)
#              Result (Table S5): C:98.32% [S:97.88%, D:0.44%], F:1.32%, M:0.37%
#              Total genes checked: 1367
#   nucmer : --maxmatch, then MUMmer / delta-filter / show-coords for synteny plots
#
# Conda env: env_qc   (quast, busco — add if not already: conda install -c bioconda quast busco)
#            env_assembly or env_mummer for nucmer
#
# HISTORICAL note:
#   The existing BUSCO summary at $BUSCO_FINAL_SUMMARY is a frozen historical
#   result that matches published Table S5 exactly. It was computed on an
#   intermediate genome (tcas6_finalGapfillScaffCorr). When this step is
#   re-run, it targets the final genome ($FINAL_GENOME_FASTA) instead.
#
# Usage:
#   bash 09_assembly_qc.sh           # check: list existing historical QC outputs
#   bash 09_assembly_qc.sh run       # run QUAST + BUSCO + nucmer on final genome
#   bash 09_assembly_qc.sh run --force
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/09_assembly_qc"

banner "STEP 09: Assembly QC (QUAST / BUSCO / nucmer) [$RUN_MODE]"
require_paths "$FINAL_GENOME_FASTA" "$PUBREF_TCAS52" \
              "$BUSCO_FINAL_SUMMARY" "$GENOME_COMPARISONS"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Historical BUSCO summary (frozen result, exact match to thesis Table S5):"
  cat "$BUSCO_FINAL_SUMMARY" 2>/dev/null | sed 's/^/  /'
  echo
  echo "QUAST historical outputs (F:):"
  find "$QUAST_HISTORICAL" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  echo
  echo "nucmer/MUMmer comparison outputs (D:):"
  find "$GENOME_COMPARISONS" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
if ! step_done "$WORKDIR" || [[ "$FORCE" == "--force" ]]; then
  mkdir -p "$WORKDIR"
  activate_conda env_qc

  REF_FASTA=$(find "$PUBREF_TCAS52" -name "*.fna" 2>/dev/null | head -1)
  [[ -z "$REF_FASTA" ]] && { echo "[ERROR] Tcas5.2 .fna not found" >&2; exit 1; }

  # --- QUAST -----------------------------------------------------------------
  echo; echo "=== QUAST ==="
  quast.py \
    "$FINAL_GENOME_FASTA" \
    --reference "$REF_FASTA" \
    --threads 16 \
    --min-contig 500 \
    -o "$WORKDIR/quast_final"
  echo "QUAST output: $WORKDIR/quast_final/report.txt"
  cat "$WORKDIR/quast_final/report.txt" | sed 's/^/  /'

  # --- BUSCO -----------------------------------------------------------------
  echo; echo "=== BUSCO (insecta_odb10) ==="
  busco \
    -i "$FINAL_GENOME_FASTA" \
    -l insecta_odb10 \
    -o busco_inTcas1_final \
    --out_path "$WORKDIR" \
    -m genome \
    -c 16 \
    --offline 2>/dev/null || \
  busco \
    -i "$FINAL_GENOME_FASTA" \
    -l insecta_odb10 \
    -o busco_inTcas1_final \
    --out_path "$WORKDIR" \
    -m genome \
    -c 16

  echo; echo "BUSCO summary:"
  cat "$WORKDIR/busco_inTcas1_final/short_summary"*.txt 2>/dev/null | sed 's/^/  /'

  echo; echo "Thesis Table S5 expected: C:98.32% [S:97.88%, D:0.44%], F:1.32%, M:0.37%"

  # --- nucmer structural comparison ------------------------------------------
  echo; echo "=== nucmer (inTcas1 vs Tcas5.2) ==="
  mkdir -p "$WORKDIR/nucmer"
  nucmer \
    --maxmatch \
    -p "$WORKDIR/nucmer/inTcas1_vs_tcas52" \
    "$REF_FASTA" \
    "$FINAL_GENOME_FASTA"
  delta-filter -r -q "$WORKDIR/nucmer/inTcas1_vs_tcas52.delta" \
    > "$WORKDIR/nucmer/inTcas1_vs_tcas52.filtered.delta"
  show-coords -rcl "$WORKDIR/nucmer/inTcas1_vs_tcas52.filtered.delta" \
    > "$WORKDIR/nucmer/inTcas1_vs_tcas52.coords"
  mummerplot \
    --filter --png \
    -p "$WORKDIR/nucmer/inTcas1_vs_tcas52" \
    "$WORKDIR/nucmer/inTcas1_vs_tcas52.filtered.delta" || true

  mark_done "$WORKDIR"
else
  echo "[SKIP] Already done: $WORKDIR/.done — pass '--force' to re-run"
fi
