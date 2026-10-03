#!/usr/bin/env bash
# ==============================================================================
# STEP 16 — Secondary synteny/orthology analysis: 1SynChro (Opscan-based)
#
# NOT reported in the thesis methods text. Included as a secondary, exploratory
# synteny/orthology analysis run alongside (not instead of) the SibeliaZ-based
# synteny in Step 15.
#
# Tool: SynChro (Opscan-based whole-genome synteny/ortholog block detection).
# Location: $SYNCHRO_TOOL_DIR (F:\P1\tools\1SynChro\1SynChro_dir\)
# Scripts: n01Genomes.py, n02Homologs.py, n03Summary.py, n04SVGPicture.py,
#           n05DotPlot.py, n06G1fG2.py, Opscan/ (must be compiled first)
# Requires: Python 2.7 + Opscan binary compiled from Opscan/Makefile
#
# To build Opscan:
#   cd $SYNCHRO_TOOL_DIR/Opscan && make
#
# Usage:
#   bash 16_synteny_1synchro_secondary.sh           # check: list tool files
#   bash 16_synteny_1synchro_secondary.sh run       # run SynChro pipeline
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
WORKDIR="$PIPELINE_ROOT/work/16_synteny_synchro"

banner "STEP 16 [SECONDARY/EXPLORATORY]: 1SynChro synteny [$RUN_MODE]"
require_paths "$SYNCHRO_TOOL_DIR" "$FINAL_GENOME_FASTA" "$PUBREF_TCAS52"

if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Tool README:"
  cat "$SYNCHRO_TOOL_DIR/README_SynChro.txt" 2>/dev/null | head -30 | sed 's/^/  /' \
    || echo "  (README not readable from here)"
  echo
  echo "SynChro tool files (F:):"
  find "$SYNCHRO_TOOL_DIR" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  echo
  echo "Opscan binary:"
  find "$SYNCHRO_TOOL_DIR/Opscan" -name "opscan" 2>/dev/null | sed 's/^/  /' \
    || echo "  (not compiled — run: cd $SYNCHRO_TOOL_DIR/Opscan && make)"
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
mkdir -p "$WORKDIR"

# Verify Opscan binary exists
OPSCAN_BIN=$(find "$SYNCHRO_TOOL_DIR/Opscan" -name "opscan" 2>/dev/null | head -1)
if [[ -z "$OPSCAN_BIN" ]]; then
  echo "[ERROR] Opscan binary not compiled."
  echo "        Run: cd $SYNCHRO_TOOL_DIR/Opscan && make"
  exit 1
fi

# Verify Python 2 available
PYTHON2=$(command -v python2 || command -v python2.7 || true)
if [[ -z "$PYTHON2" ]]; then
  echo "[ERROR] Python 2.7 not found. SynChro requires Python 2."
  echo "        conda create -n env_python2 python=2.7"
  echo "        then: conda run -n env_python2 python <script>"
  exit 1
fi

REF_FASTA=$(find "$PUBREF_TCAS52" -name "*.fna" 2>/dev/null | head -1)
REF_GFF=$(find "$PUBREF_TCAS52" -name "*.gff" -o -name "*.gff3" 2>/dev/null | head -1)
MAKER_GFF="$FINAL_ANNOTATION_MAKER_GFF"

echo "Running SynChro pipeline (n01 → n06)..."
cd "$SYNCHRO_TOOL_DIR"

"$PYTHON2" n01Genomes.py inTcas1 Tcas5.2 \
  --genome1 "$FINAL_GENOME_FASTA" --gff1 "$MAKER_GFF" \
  --genome2 "$REF_FASTA" --gff2 "$REF_GFF" \
  --outdir "$WORKDIR" 2>/dev/null || \
"$PYTHON2" n01Genomes.py inTcas1 Tcas5.2

"$PYTHON2" n02Homologs.py inTcas1 Tcas5.2 2>/dev/null || "$PYTHON2" n02Homologs.py inTcas1 Tcas5.2
"$PYTHON2" n03Summary.py  inTcas1 Tcas5.2 2>/dev/null || "$PYTHON2" n03Summary.py  inTcas1 Tcas5.2
"$PYTHON2" n04SVGPicture.py inTcas1 Tcas5.2 2>/dev/null || "$PYTHON2" n04SVGPicture.py inTcas1 Tcas5.2
"$PYTHON2" n05DotPlot.py  inTcas1 Tcas5.2 2>/dev/null || "$PYTHON2" n05DotPlot.py  inTcas1 Tcas5.2
"$PYTHON2" n06G1fG2.py    inTcas1 Tcas5.2 2>/dev/null || "$PYTHON2" n06G1fG2.py    inTcas1 Tcas5.2

echo
echo "[Step 16 complete — SynChro outputs in $SYNCHRO_TOOL_DIR (tool's own output location)]"
echo "[SECONDARY/EXPLORATORY — not part of the published thesis pipeline]"
