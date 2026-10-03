#!/usr/bin/env bash
# ==============================================================================
# STEP 07 — Final assembled genome checkpoint (inTcas1)
#
# Thesis (Methods p.18, Table 1):
#   After gap-filling, all short reads from all 8 individuals were merged for
#   one final Pilon polishing round → inTcas1 reference.
#   Final genome: 165.5 Mbp, 8551 scaffolds (11 chromosomes + 8540 unplaced),
#   N50 = 14.47 Mbp, L50 = 5, BUSCO = 98.3% (insecta_odb10).
#
# CANONICAL FILE:
#   $FINAL_GENOME_FASTA = E:\P1\_pipeline_steps\07_final_genome\inTcas1_final_removedContam.fasta
#   This is the genome that produced the published Table 1 — near-exact match
#   to all 6 reported metrics.
#
# This step computes assembly statistics on the confirmed final genome and
# compares them to the thesis Table 1 values.
#
# Usage:
#   bash 07_final_genome.sh        # verify + print stats (no tool run)
#   bash 07_final_genome.sh run    # same (stats-only step, always runs)
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

WORKDIR="$PIPELINE_ROOT/work/07_final_genome"
mkdir -p "$WORKDIR"

banner "STEP 07: Final genome checkpoint (inTcas1)"
require_paths "$FINAL_GENOME_FASTA" "$FINAL_GENOME_TCAS6_FINAL_ALL"

echo
echo "CONFIRMED final genome (E:):"
ls -lh "$FINAL_GENOME_FASTA" 2>/dev/null | sed 's/^/  /'

echo
echo "Assembly statistics (computed directly from FASTA):"
python3 "$PIPELINE_ROOT/work/compute_assembly_stats.py" "$FINAL_GENOME_FASTA" 2>/dev/null \
  || fasta_summary "$FINAL_GENOME_FASTA"

echo
echo "Thesis Table 1 published values for comparison:"
cat <<'TABLE'
  Total length    : 165,527,733 bp
  Scaffolds (total): 8551  (11 chromosomes + 8540 unplaced)
  Largest scaffold : 20,697,200 bp (LG4)
  N50              : 14.47 Mbp
  L50              : 5
  Ns               : 2.86 Mbp
TABLE

echo
echo "Earlier candidate (tcas6_final_all.fasta, F: — paired with annotation, kept for Step 13):"
ls -lh "$FINAL_GENOME_TCAS6_FINAL_ALL" 2>/dev/null | sed 's/^/  /'

echo
echo "Per-individual final assemblies (D:, used for pangenome — see improvement report):"
find "$FINAL_GENOMES_FASTA" -maxdepth 1 -name "*.fasta" 2>/dev/null | sort | sed 's/^/  /'

echo
echo "[Step 07 complete — stats logged to $WORKDIR]"
