#!/usr/bin/env bash
# ==============================================================================
# STEP 15 — Whole-genome synteny (SibeliaZ) + Figure 1 generation (R scripts)
#
# Thesis (Methods p.18, Results p.21, Figure 1A-D):
#   SibeliaZ (Minkin & Medvedev, 2020) with LCB parameters:
#     -k 25 (k-mer size)   -t 16 (threads)
#     Blocks reported at >95% identity and >60% coverage.
#   Comparisons: inTcas1 vs Tcas5.2, vs T.freemani (Tfree1.0), vs T.confusum (Tcon1.0)
#   Output → maf2synteny → LCB table → R scripts → Figure 1A-D, 1E
#
# Figure-generation R scripts (F:, confirmed as source of published Fig 1):
#   $FIGURE_SCRIPTS_DIR/dotplot.R       → Fig 1C-D
#   $FIGURE_SCRIPTS_DIR/synteny_LCB.R  → Fig 1A (LCB proportion heatmap)
#   $FIGURE_SCRIPTS_DIR/ideogram.R     → Fig 1B (sequence addition/deletion schematic)
#   $FIGURE_SCRIPTS_DIR/y_chr_synteny.R → Fig 1E (Y chromosome)
#   $FIGURE_SCRIPTS_DIR/chrY_consensus.R
#
# Conda env: env_sibeliaz (conda create -n env_sibeliaz -c bioconda sibeliaz)
# R: r-base + required packages (ggplot2, ggbio, etc.) in env_sibeliaz or env_assembly
#
# Existing SibeliaZ output: $SIBELIAZ_OUTPUT (E:)
#
# Usage:
#   bash 15_synteny_sibeliaz.sh           # check: list existing outputs + R scripts
#   bash 15_synteny_sibeliaz.sh run       # run SibeliaZ + R figure scripts
#   bash 15_synteny_sibeliaz.sh run --force
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/15_synteny_sibeliaz"

# Comparison genomes (Tcas5.2 is the reference; add Tfree/Tcon if available)
TFREE_FASTA=$(find "$TCON_TFREE_DATASETS" -name "*freemani*.fna" -o -name "*Tfree*.fna" 2>/dev/null | head -1)
TCON_FASTA=$(find "$TCON_TFREE_DATASETS"  -name "*confusum*.fna" -o -name "*Tcon*.fna"  2>/dev/null | head -1)

banner "STEP 15: SibeliaZ synteny + Figure 1 R scripts [$RUN_MODE]"
require_paths "$SIBELIAZ_OUTPUT" "$FIGURE_SCRIPTS_DIR" "$FINAL_GENOME_FASTA" "$PUBREF_TCAS52"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing SibeliaZ output (E:):"
  find "$SIBELIAZ_OUTPUT" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  echo
  echo "Figure generation scripts (F:):"
  find "$FIGURE_SCRIPTS_DIR" -maxdepth 2 -name "*.R" -o -name "*.py" 2>/dev/null | sort | sed 's/^/  /'
  echo
  echo "T.freemani fasta: ${TFREE_FASTA:-(not found in $TCON_TFREE_DATASETS)}"
  echo "T.confusum fasta: ${TCON_FASTA:-(not found in $TCON_TFREE_DATASETS)}"
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

activate_conda env_sibeliaz

REF_FASTA=$(find "$PUBREF_TCAS52" -name "*.fna" 2>/dev/null | head -1)

# Build genome list: inTcas1 + Tcas5.2 (+ Tfree/Tcon if found)
GENOME_LIST=("$FINAL_GENOME_FASTA" "$REF_FASTA")
[[ -n "$TFREE_FASTA" && -f "$TFREE_FASTA" ]] && GENOME_LIST+=("$TFREE_FASTA")
[[ -n "$TCON_FASTA"  && -f "$TCON_FASTA"  ]] && GENOME_LIST+=("$TCON_FASTA")

echo "Genomes for SibeliaZ comparison:"
printf "  %s\n" "${GENOME_LIST[@]}"

# --- SibeliaZ ----------------------------------------------------------------
echo; echo "=== SibeliaZ (k=25) ==="
sibeliaz \
  -k 25 \
  -t 16 \
  -o "$WORKDIR/sibeliaz_out" \
  "${GENOME_LIST[@]}"

# --- maf2synteny (convert MAF blocks to synteny coordinates) ----------------
echo; echo "=== maf2synteny (convert to synteny table) ==="
if command -v maf2synteny &>/dev/null; then
  maf2synteny \
    "$WORKDIR/sibeliaz_out/blocks_coords.gff" \
    -b 100 \
    -o "$WORKDIR/synteny_blocks"
else
  echo "  [WARN] maf2synteny not found — parse $WORKDIR/sibeliaz_out/blocks_coords.gff manually"
fi

# --- Figure generation (R scripts) ------------------------------------------
echo; echo "=== Figure 1 generation (R scripts) ==="
for Rscript_name in dotplot.R synteny_LCB.R ideogram.R y_chr_synteny.R chrY_consensus.R; do
  RSCRIPT=$(find "$FIGURE_SCRIPTS_DIR" -name "$Rscript_name" 2>/dev/null | head -1)
  if [[ -z "$RSCRIPT" ]]; then
    echo "  [SKIP] $Rscript_name not found under $FIGURE_SCRIPTS_DIR"
    continue
  fi
  echo "  Running: $Rscript_name"
  Rscript "$RSCRIPT" \
    --sibeliaz_out "$WORKDIR/sibeliaz_out" \
    --out_dir "$WORKDIR/figures" 2>/dev/null || \
  Rscript "$RSCRIPT" 2>/dev/null || \
  echo "    [WARN] $Rscript_name requires manual argument setup — see script header"
done

echo
echo "SibeliaZ outputs: $WORKDIR/sibeliaz_out/"
echo "Figures:          $WORKDIR/figures/"

mark_done "$WORKDIR"
