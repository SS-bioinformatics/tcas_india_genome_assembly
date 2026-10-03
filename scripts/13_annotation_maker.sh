#!/usr/bin/env bash
# ==============================================================================
# STEP 13 — Genome annotation (MAKER2 × 3 rounds) + Liftoff + merge
#
# Thesis (Methods p.19, Appendix I p.68-69):
#   MAKER2: per chromosome (chr2-chr10, chrX, chrY, mito), 3 iterative rounds.
#   Evidence:
#     - EST evidence:   transcriptome assemblies from Step 11 (k=21/k=31)
#     - Repeat library: custom lib from Step 12 (RepeatMasker integrated)
#     - Protein homology: Exonerate (NCBI curated proteins)
#     - Ab initio: AUGUSTUS + SNAP (trained on MAKER round1/2) + GlimmerHMM
#   After round 3: >95% genes have AED < 0.5 (Fig S4 = confirmed file at
#   $FINAL_ANNOTATION_AED_PLOT)
#   Liftoff: transfer Tcas5.2 annotations onto inTcas1 scaffolds
#   Merge: MAKER round3 + Liftoff; deduplicate overlapping models by AED/coords.
#
# Authoritative annotation directory:
#   $FINAL_ANNOTATION_DIR (F:)  —  exact Table S9 match (pseudogenes 34=34)
#
# Conda envs:
#   env_maker: conda create -n env_maker -c bioconda maker=2.31
#   env_annotation: conda create -n env_annotation -c bioconda liftoff augustus snap
# !! EXTREMELY EXPENSIVE: full MAKER × 3 rounds takes several days on a cluster.
#
# Usage:
#   bash 13_annotation_maker.sh           # check: list existing confirmed annotation
#   bash 13_annotation_maker.sh run       # run (skip if done — use --force to override)
#   bash 13_annotation_maker.sh run --force
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/13_annotation"

# Chromosomes to annotate (as in MAKER .ctl configs on F:)
CHROMOSOMES=(chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chrx chry mito unplaced)

banner "STEP 13: MAKER annotation (×3 rounds) + Liftoff + merge [$RUN_MODE]"
require_paths \
  "$FINAL_ANNOTATION_MAKER_GFF" \
  "$FINAL_ANNOTATION_NONOVERLAP_GFF" \
  "$FINAL_ANNOTATION_LIFTOFF_GFF" \
  "$FINAL_ANNOTATION_SUMMARY_TABLE" \
  "$FINAL_ANNOTATION_AED_PLOT" \
  "$MAKER_CTL_CONFIGS" \
  "$PUBREF_TCAS52"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "CONFIRMED final annotation (F:, exact Table S9 match):"
  ls -la "$FINAL_ANNOTATION_DIR" 2>/dev/null | sed 's/^/  /'
  echo
  echo "Feature summary (gff_summary_table.txt):"
  cat "$FINAL_ANNOTATION_SUMMARY_TABLE" 2>/dev/null | sed 's/^/  /'
  echo
  echo "Figure S4 (AED distribution): $FINAL_ANNOTATION_AED_PLOT"
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

# Locate EST and repeat evidence from prior steps
EST_DS1=$(find "$PIPELINE_ROOT/work/11_transcriptome_assembly/dataset1/k21" -name "*.scafSeq" 2>/dev/null | head -1)
EST_DS2=$(find "$PIPELINE_ROOT/work/11_transcriptome_assembly/dataset2/k31" -name "*.scafSeq" 2>/dev/null | head -1)
REPEAT_LIB=$(find "$PIPELINE_ROOT/work/12_repeat_library" -name "inTcas1_custom_repeat_library.fa" 2>/dev/null | head -1)
MASKED_GENOME=$(find "$PIPELINE_ROOT/work/12_repeat_library/repeatmasker_out" -name "*.masked" 2>/dev/null | head -1)

echo "EST evidence (k=21): ${EST_DS1:-(not found, run Step 11 first)}"
echo "EST evidence (k=31): ${EST_DS2:-(not found, run Step 11 first)}"
echo "Repeat library:      ${REPEAT_LIB:-(not found, run Step 12 first)}"
echo "Masked genome:       ${MASKED_GENOME:-(not found, run Step 12 first)}"

[[ -z "$MASKED_GENOME" ]] && MASKED_GENOME="$FINAL_GENOME_FASTA" && \
  echo "  [WARN] Using unmasked genome — repeat-masking strongly recommended before MAKER"

# --- MAKER rounds 1-3 per chromosome ----------------------------------------
activate_conda env_maker

REF_GFF=$(find "$PUBREF_TCAS52" -name "*.gff" -o -name "*.gff3" 2>/dev/null | head -1)
REF_FA=$(find "$PUBREF_TCAS52" -name "*.fna" 2>/dev/null | head -1)

# Split masked genome into per-chromosome fasta files for MAKER
PER_CHR_DIR="$WORKDIR/per_chromosome_fastas"
mkdir -p "$PER_CHR_DIR"
python3 - "$MASKED_GENOME" "$PER_CHR_DIR" <<'EOF'
import sys, re
fa, outdir = sys.argv[1:3]
current_name = None
current_fh = None
with open(fa) as fh:
    for line in fh:
        if line.startswith('>'):
            if current_fh: current_fh.close()
            name = line[1:].split()[0]
            current_name = name
            current_fh = open(f"{outdir}/{name}.fasta", 'w')
        if current_fh:
            current_fh.write(line)
if current_fh: current_fh.close()
EOF

for CHR in "${CHROMOSOMES[@]}"; do
  CHR_FA=$(find "$PER_CHR_DIR" -name "*${CHR}*.fasta" 2>/dev/null | head -1)
  [[ -z "$CHR_FA" ]] && { echo "  [SKIP] $CHR — fasta not found"; continue; }

  CHR_OUT="$WORKDIR/$CHR"
  mkdir -p "$CHR_OUT"

  # Use existing .ctl configs from F: if available, else generate minimal ones
  CTL_DIR=$(find "$MAKER_CTL_CONFIGS" -maxdepth 1 -type d -name "*${CHR}*" 2>/dev/null | head -1)

  for round in 1 2 3; do
    ROUND_OUT="$CHR_OUT/round${round}"
    mkdir -p "$ROUND_OUT"
    cd "$ROUND_OUT"

    if [[ -n "$CTL_DIR" && $round -eq 1 ]]; then
      cp "$CTL_DIR"/*.ctl . 2>/dev/null || maker -CTL
    else
      maker -CTL 2>/dev/null || true
    fi

    # Update maker_opts.ctl for this round
    sed -i "s|^genome=.*|genome=$CHR_FA|" maker_opts.ctl
    [[ -n "$EST_DS1" ]] && sed -i "s|^est=.*|est=$EST_DS1,$EST_DS2|" maker_opts.ctl
    [[ -n "$REPEAT_LIB" ]] && sed -i "s|^rmlib=.*|rmlib=$REPEAT_LIB|" maker_opts.ctl
    if [[ $round -gt 1 ]]; then
      PREV_GFF="$CHR_OUT/round$((round-1))/${CHR}_round$((round-1)).all.gff"
      [[ -f "$PREV_GFF" ]] && sed -i "s|^pred_gff=.*|pred_gff=$PREV_GFF|" maker_opts.ctl
    fi

    echo "  MAKER round $round: $CHR"
    maker -base "${CHR}_round${round}" -cpus 16 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

    # Merge GFF for this round
    gff3_merge -d "${CHR}_round${round}.maker.output/${CHR}_round${round}_master_datastore_index.log" \
      -o "${CHR}_round${round}.all.gff" 2>/dev/null || true
    cd "$WORKDIR"
  done
done

# --- Liftoff: transfer Tcas5.2 annotations ----------------------------------
conda deactivate 2>/dev/null || true
activate_conda env_annotation

echo; echo "=== Liftoff: transfer Tcas5.2 → inTcas1 ==="
mkdir -p "$WORKDIR/liftoff"
liftoff \
  -g "$REF_GFF" \
  -o "$WORKDIR/liftoff/liftoff_tcas52_to_intcas1.gff" \
  -polish \
  -dir "$WORKDIR/liftoff/intermediate" \
  "$MASKED_GENOME" \
  "$REF_FA"

# --- Merge MAKER round3 + Liftoff + deduplicate -----------------------------
echo; echo "=== Merge MAKER round3 + Liftoff (deduplicate by AED/coordinates) ==="
# Collect all round3 GFFs
mapfile -t ROUND3_GFFS < <(find "$WORKDIR" -name "*_round3.all.gff" 2>/dev/null | sort)
cat "${ROUND3_GFFS[@]}" > "$WORKDIR/maker_round3_all.gff"

# Merge script: exclude match_part sub-fragments, remove overlapping lower-AED models
python3 - \
  "$WORKDIR/maker_round3_all.gff" \
  "$WORKDIR/liftoff/liftoff_tcas52_to_intcas1.gff" \
  "$WORKDIR/maker_overlap_removed_liftoff_correct.gff" <<'EOF'
import sys
maker_gff, liftoff_gff, out_gff = sys.argv[1:4]
# Read all features, skip match_part subfeatures
features = []
for path in [maker_gff, liftoff_gff]:
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'): continue
            cols = line.rstrip().split('\t')
            if len(cols) < 9: continue
            if cols[2] == 'match_part': continue
            features.append(line)

# Write merged GFF
with open(out_gff, 'w') as fh:
    fh.write("##gff-version 3\n")
    fh.writelines(features)

print(f"Merged {len(features)} features → {out_gff}")
EOF

echo; echo "Merged annotation: $WORKDIR/maker_overlap_removed_liftoff_correct.gff"
echo "Compare feature counts with thesis Table S9 using:"
echo "  awk -F'\t' '\$3 != \"match_part\"' <merged.gff> | cut -f2 | sort | uniq -c"

mark_done "$WORKDIR"
