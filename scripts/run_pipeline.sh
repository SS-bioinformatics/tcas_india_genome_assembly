#!/usr/bin/env bash
# ==============================================================================
# run_pipeline.sh — master runner for the Chapter 2 (inTcas1) pipeline
#
# Usage:
#   ./run_pipeline.sh                       # check mode: verify all file paths
#   ./run_pipeline.sh check                 # same as above
#   ./run_pipeline.sh run                   # run all steps (expensive steps skip if done)
#   ./run_pipeline.sh check 05              # check mode, step 05 only
#   ./run_pipeline.sh run 05                # run mode, step 05 only
#   ./run_pipeline.sh run 05 --force        # force re-run step 05 even if done
#
# Execution environment:
#   Run from WSL2 Ubuntu (tools live in conda envs inside WSL2).
#   Paths on D/E/F are accessed as /mnt/d, /mnt/e, /mnt/f (auto-detected in paths.config.sh).
#   Scripts on C: are at /mnt/c/Users/csc/Desktop/PhD_work/ch1/scripts/
#
# Pipeline steps:
#   00  Raw data inventory (no tools)
#   01  FastQC + Trimmomatic + Jellyfish + GenomeScope2      [env_qc]
#   02  Canu long-read correction                             [env_canu]       !! expensive
#   03  SPAdes draft assembly + Purge Haplotigs round 1       [env_assembly]   !! expensive
#   04  Pilon polishing (×3) + Purge Haplotigs round 2        [env_assembly]   !! expensive
#   05  RagTag reference-guided scaffolding (1A vs Tcas5.2)   [env_assembly]
#   06  TGS-GapCloser iterative gap-fill (7 steps)            [env_tgsgapcloser]
#   07  Final genome checkpoint + Table 1 stats (no tools)
#   08  FCS-GX contamination screening                        [external/NCBI DB]
#   09  QUAST + BUSCO + nucmer/MUMmer                         [env_qc]
#   10  DiscoverY + RagTag + BLASTn + findZX (Y chromosome)  [env_assembly]
#   11  SOAPdenovo-Trans transcriptome assembly               [env_transcriptome]
#   12  RepeatModeler + Blastx filter + RepeatMasker          [env_repeatmodeler/repeatmasker]
#   13  MAKER2 (×3 rounds) + Liftoff + merge                  [env_maker/annotation] !! expensive
#   14  Final annotation curation                             [python3]
#   15  SibeliaZ synteny + Figure 1 R scripts                 [env_sibeliaz]
#   16  1SynChro synteny (secondary/exploratory, not in thesis)[python2 + Opscan]
# ==============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

MODE="${1:-check}"
ONLY="${2:-}"
FORCE="${3:-}"

STEPS=(
  00_raw_data_check.sh
  01_qc_kmer.sh
  02_longread_correction.sh
  03_draft_assembly.sh
  04_polishing.sh
  05_scaffolding_ragtag.sh
  06_gapfilling.sh
  07_final_genome.sh
  08_contamination_screening.sh
  09_assembly_qc.sh
  10_y_chromosome.sh
  11_transcriptome_assembly.sh
  12_repeat_library.sh
  13_annotation_maker.sh
  14_final_curation.sh
  15_synteny_sibeliaz.sh
  16_synteny_1synchro_secondary.sh
)

source ../config/paths.config.sh 2>/dev/null || true
mkdir -p ../logs
LOGFILE="../logs/pipeline_${MODE}_$(date +%Y%m%d_%H%M%S).log"
echo "Log: $LOGFILE"
echo "Mode: $MODE | Filter: ${ONLY:-(all)} | Force: ${FORCE:-(no)}"
echo "================================================================"

FAILED=0
for s in "${STEPS[@]}"; do
  if [[ -n "$ONLY" ]] && [[ "$s" != "$ONLY"* ]]; then
    continue
  fi
  echo
  echo "### $s [$MODE] ###" | tee -a "$LOGFILE"
  if bash "./$s" "$MODE" "$FORCE" 2>&1 | tee -a "$LOGFILE"; then
    echo "  [OK] $s" | tee -a "$LOGFILE"
  else
    echo "  [FAIL] $s — see log" | tee -a "$LOGFILE"
    FAILED=$((FAILED + 1))
  fi
done

echo
echo "================================================================"
if [[ $FAILED -eq 0 ]]; then
  echo "All steps completed successfully."
else
  echo "$FAILED step(s) failed — review $LOGFILE"
fi
echo "Full log: $LOGFILE"
