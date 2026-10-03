#!/usr/bin/env bash
# ==============================================================================
# STEP 08 — Contamination screening (FCS-GX, NCBI)
#
# Thesis (Methods p.18): confirmed no contamination in inTcas1 using FCS-GX
# (Astashyn et al., 2023) against the full NCBI genome database. T. castaneum
# taxid = 7070.
#
# Actual workflow (6 iterative rounds, more granular than thesis text):
#   inTcas1_draft.fasta
#     → FCS-GX screen → Contamination1.txt
#     → exclude_contam_2.txt → inTcas1_contam_removed_1.fasta
#     → FCS-GX screen → Contamination2.txt
#     → ... (6 rounds total)
#     → inTcas1_final_removedContam.fasta  ← CONFIRMED final genome (Step 07)
#
# All intermediate files exist at: $CONTAM_WORKDIR (E:)
# FCS-GX requires a large local database (multi-GB); run on a machine where it
# is already staged. The existing 6-round trail is the definitive result.
#
# Usage:
#   bash 08_contamination_screening.sh        # check: list contamination trail
#   bash 08_contamination_screening.sh run    # re-run (requires FCS-GX DB)
#   bash 08_contamination_screening.sh run --force
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/08_contamination_screening"

TAXID=7070   # Tribolium castaneum NCBI taxid

banner "STEP 08: Contamination screening (FCS-GX) [$RUN_MODE]"
require_paths "$CONTAM_WORKDIR" "$CONTAM_FCSGX_DIR"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing contamination screening trail (E:):"
  find "$CONTAM_WORKDIR" -maxdepth 1 \
    \( -name "Contamination*.txt" \
    -o -name "exclude*contam*.txt" \
    -o -name "trim*contam*" \
    -o -name "*contam_removed*.fasta" \
    -o -name "inTcas1_final_removedContam*" \
    -o -name "convert_contam_to_N.sh" \
    -o -name "remove_sequence.sh" \) 2>/dev/null \
    | sort | sed 's/^/  /'
  echo
  echo "Final output (confirmed final genome):"
  ls -lh "$FINAL_GENOME_FASTA" 2>/dev/null | sed 's/^/  /'
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

# Locate FCS-GX Python wrapper
FCS_PY=$(find "$CONTAM_FCSGX_DIR" -name "fcs.py" 2>/dev/null | head -1)
[[ -z "$FCS_PY" ]] && FCS_PY=$(command -v fcs.py 2>/dev/null || true)
[[ -z "$FCS_PY" ]] && { echo "[ERROR] fcs.py not found. Check $CONTAM_FCSGX_DIR or PATH." >&2; exit 1; }

# FCS-GX database path (must exist locally)
GXDB="${FCS_GX_DB:-}"
if [[ -z "$GXDB" ]]; then
  echo "[ERROR] Set the FCS_GX_DB environment variable to the path of the GX database." >&2
  echo "        e.g.  export FCS_GX_DB=/path/to/gxdb" >&2
  exit 1
fi

# Starting draft: output of Step 06 (or Step 07 draft before contamination removal)
DRAFT="$PIPELINE_ROOT/work/06_gapfilling/fill_ONT/filled.scaff_seqs"
[[ -f "$DRAFT" ]] || { echo "[ERROR] Step 06 output not found. Run Step 06 first." >&2; exit 1; }

echo "Input draft: $DRAFT"
echo "FCS-GX taxid: $TAXID"
echo "FCS-GX DB: $GXDB"

OUT="$WORKDIR/fcsgx_round1"
mkdir -p "$OUT"

python3 "$FCS_PY" screen genome \
  --fasta "$DRAFT" \
  --out-dir "$OUT" \
  --gx-db "$GXDB" \
  --tax-id "$TAXID"

echo
echo "FCS-GX screen output: $OUT"
echo "Review the report, create exclude_contam_2.txt, then run the removal step:"
echo "  bash $CONTAM_WORKDIR/remove_sequence.sh  <draft>  <exclude_list>  <output.fasta>"
echo
echo "Repeat for each contaminated sequence batch (thesis used 6 rounds)."
echo "Final clean genome → $WORKDIR/inTcas1_final_removedContam.fasta"

mark_done "$WORKDIR"
