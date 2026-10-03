#!/usr/bin/env bash
# ==============================================================================
# STEP 06 — Iterative gap filling (TGS-GapCloser, 7-step Table S6 matrix)
#
# Thesis (Methods p.18, Appendix I p.67, Table S6):
#   Gaps in the scaffolded 1A genome were filled iteratively using 6 other
#   draft genomes as sequence donors, then finally with Canu-corrected ONT reads.
#   Excluded: 1B and 18B (low size / high fragmentation).
#   Fill order (Table S6): 12B → 12D → 12E → 13C → 13E → 18A → ONT
#   Flag --ne (no extension beyond gap) throughout.
#   ONT fill: --tgstype ont (long-read mode).
#   Validated with GapCloser, LR-Gapcloser, GAPPadder (17 python scripts on D:).
#
# Tool: TGS-GapCloser  (conda create -n env_tgsgapcloser -c bioconda tgsgapcloser)
# Existing outputs: $GAPFILL_OUTPUT (D:)
#
# !! EXPENSIVE: 7 sequential fill rounds, each re-aligns all reads.
#
# Usage:
#   bash 06_gapfilling.sh           # check: list existing outputs
#   bash 06_gapfilling.sh run       # run (skip if done)
#   bash 06_gapfilling.sh run --force
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/06_gapfilling"

# Gap-fill donor order (Table S6 — 18A used, 18B excluded)
FILL_ORDER=(12B 12D 12E 13C 13E 18A)

banner "STEP 06: TGS-GapCloser iterative gap filling (7 steps) [$RUN_MODE]"
require_paths "$GAPFILL_OUTPUT" "$NANOPORE_CANU_CORRECTED" "$FINAL_GENOMES_FASTA"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing gap-fill outputs (D:):"
  find "$GAPFILL_OUTPUT" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

activate_conda env_tgsgapcloser

# Starting scaffold: output of Step 05
CURRENT_SCAFF="$PIPELINE_ROOT/work/05_scaffolding_ragtag/ragtag_output/ragtag.scaffold.fasta"
if [[ ! -f "$CURRENT_SCAFF" ]]; then
  CURRENT_SCAFF=$(find "$RAGTAG_OUTPUT" -name "ragtag.scaffold.fasta" 2>/dev/null | head -1)
fi
[[ -f "$CURRENT_SCAFF" ]] || { echo "[ERROR] Step 05 scaffold not found. Run Step 05 first." >&2; exit 1; }
echo "Starting scaffold: $CURRENT_SCAFF"

# Fills 1–6: draft genome donors (short reads, default mode)
for DONOR_NAME in "${FILL_ORDER[@]}"; do
  DONOR_FASTA=$(find "$FINAL_GENOMES_FASTA" -name "*${DONOR_NAME}*.fasta" 2>/dev/null | head -1)
  if [[ -z "$DONOR_FASTA" ]]; then
    # Try purge_haplotigs outputs if available
    DONOR_FASTA=$(find "$PURGE_HAPLOTIGS_OUTPUT" -name "*${DONOR_NAME}*purged*.fasta" 2>/dev/null | head -1)
  fi
  if [[ -z "$DONOR_FASTA" ]]; then
    echo "  [WARN] Donor $DONOR_NAME not found — skipping this fill step"
    continue
  fi

  echo; echo "=== Gap fill with $DONOR_NAME draft genome ==="
  OUT_PREFIX="$WORKDIR/fill_${DONOR_NAME}"
  mkdir -p "$OUT_PREFIX"

  tgsgapcloser \
    --scaff "$CURRENT_SCAFF" \
    --reads "$DONOR_FASTA" \
    --output "$OUT_PREFIX/filled" \
    --ne \
    --thread 16

  NEXT="$OUT_PREFIX/filled.scaff_seqs"
  [[ -f "$NEXT" ]] && CURRENT_SCAFF="$NEXT" || echo "  [WARN] TGS-GapCloser produced no output for $DONOR_NAME"
done

# Fill 7: Canu-corrected ONT long reads
echo; echo "=== Gap fill step 7: ONT long reads (Canu-corrected) ==="
ONT_READS=$(find "$NANOPORE_CANU_CORRECTED" -name "*.fasta.gz" -o -name "*.fasta" 2>/dev/null | head -1)
[[ -z "$ONT_READS" ]] && { echo "[ERROR] Corrected ONT reads not found in $NANOPORE_CANU_CORRECTED" >&2; exit 1; }

OUT_ONT="$WORKDIR/fill_ONT"
mkdir -p "$OUT_ONT"

tgsgapcloser \
  --scaff "$CURRENT_SCAFF" \
  --reads "$ONT_READS" \
  --output "$OUT_ONT/filled" \
  --ne \
  --tgstype ont \
  --thread 16

echo
echo "Final gap-filled assembly: $OUT_ONT/filled.scaff_seqs"
fasta_summary "$OUT_ONT/filled.scaff_seqs"

mark_done "$WORKDIR"
