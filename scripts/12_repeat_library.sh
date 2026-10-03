#!/usr/bin/env bash
# ==============================================================================
# STEP 12 — Custom repeat library (RepeatModeler → Blastx filter → RepeatMasker)
#
# Thesis (Appendix I p.68, Table S8):
#   1. RepeatModeler v1.0.7:  BuildDatabase + RepeatModeler -pa 32 -engine ncbi
#      Produces: inTcas1-families.fa (de novo consensus repeat models)
#   2. Blastx filter: remove potential protein-coding genes misclassified as repeats
#      blastx vs GenBank nr-Arthropoda subset (e-value < 1e-3)
#   3. Merge with known transposons from RepBase
#   4. RepeatMasker v4.0.7: -pa 24 -e ncbi -gccalc -lib <custom_library>
#      Result: 69% repeat content (Table S8)
#
# Existing RepeatModeler output: $REPEATMODELER_OUTPUT (D:)
# NOTE: The Blastx-filter + RepBase merge (step 2-3) are not isolated as a
# standalone script on disk. Commands are reproduced here from the thesis text.
# This step also resolves MISMATCH M6 if re-run with the correct custom library.
#
# Conda envs:
#   env_repeatmodeler: conda create -n env_repeatmodeler -c bioconda repeatmodeler=1.0.7
#   env_repeatmasker:  conda create -n env_repeatmasker  -c bioconda repeatmasker=4.0.7
#   BLAST in env_repeatmodeler for the Blastx filter.
#   RepBase: obtain separately (RepeatMasker license); provide path via $REPBASE_LIB.
#
# Usage:
#   bash 12_repeat_library.sh           # check: list existing RepeatModeler outputs
#   bash 12_repeat_library.sh run       # run full pipeline
#   bash 12_repeat_library.sh run --force
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/12_repeat_library"

# Paths to optional external resources
REPBASE_LIB="${REPBASE_LIB:-}"            # set externally: export REPBASE_LIB=/path/to/RepBase.fa
NR_ARTHROPODA_DB="${NR_ARTHROPODA_DB:-}"  # set externally: export NR_ARTHROPODA_DB=/path/to/nr_arthropoda

banner "STEP 12: Custom repeat library (RepeatModeler → Blastx → RepeatMasker) [$RUN_MODE]"
require_paths "$FINAL_GENOME_FASTA" "$REPEATMODELER_OUTPUT"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing RepeatModeler outputs (D:):"
  find "$REPEATMODELER_OUTPUT" -maxdepth 2 2>/dev/null | sort | sed 's/^/  /'
  echo
  echo "RepBase path: ${REPBASE_LIB:-(not set — export REPBASE_LIB=...)}"
  echo "nr-Arthropoda BLAST db: ${NR_ARTHROPODA_DB:-(not set — export NR_ARTHROPODA_DB=...)}"
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

# --- Step 1: RepeatModeler --------------------------------------------------
activate_conda env_repeatmodeler

echo; echo "=== RepeatModeler v1.0.7 ==="
cd "$WORKDIR"
BuildDatabase -name inTcas1 -engine rmblast "$FINAL_GENOME_FASTA"
RepeatModeler \
  -pa 32 \
  -engine ncbi \
  -database inTcas1 \
  2>&1 | tee "$WORKDIR/repeatmodeler.log"

FAMILIES="$WORKDIR/inTcas1-families.fa"
[[ -f "$FAMILIES" ]] || { echo "[ERROR] RepeatModeler did not produce inTcas1-families.fa" >&2; exit 1; }

# --- Step 2: Blastx filter --------------------------------------------------
echo; echo "=== Blastx filter (protein-coding gene removal, e-value < 1e-3) ==="
if [[ -n "$NR_ARTHROPODA_DB" && -f "${NR_ARTHROPODA_DB}.pin" ]]; then
  blastx \
    -query "$FAMILIES" \
    -db "$NR_ARTHROPODA_DB" \
    -evalue 1e-3 \
    -outfmt 6 \
    -num_threads 16 \
    -out "$WORKDIR/repeat_vs_nr_arthropoda.blastx"

  # Remove hits: extract query IDs with blast hits
  awk '{print $1}' "$WORKDIR/repeat_vs_nr_arthropoda.blastx" | sort -u \
    > "$WORKDIR/repeats_protein_hits.txt"
  echo "  Repeat models matching proteins: $(wc -l < "$WORKDIR/repeats_protein_hits.txt")"

  python3 - "$FAMILIES" "$WORKDIR/repeats_protein_hits.txt" "$WORKDIR/inTcas1-families_filtered.fa" <<'EOF'
import sys
fa, hits_file, out = sys.argv[1:4]
remove = set(open(hits_file).read().split())
write = False
with open(fa) as fin, open(out,'w') as fout:
    for line in fin:
        if line.startswith('>'):
            write = line[1:].split()[0] not in remove
        if write:
            fout.write(line)
EOF
  FILTERED="$WORKDIR/inTcas1-families_filtered.fa"
else
  echo "  [WARN] NR_ARTHROPODA_DB not set or not found — skipping Blastx filter."
  echo "         Set: export NR_ARTHROPODA_DB=/path/to/nr_arthropoda_db"
  FILTERED="$FAMILIES"
fi

# --- Step 3: Merge with RepBase known transposons ---------------------------
echo; echo "=== Merge with RepBase transposons ==="
CUSTOM_LIB="$WORKDIR/inTcas1_custom_repeat_library.fa"
if [[ -n "$REPBASE_LIB" && -f "$REPBASE_LIB" ]]; then
  cat "$FILTERED" "$REPBASE_LIB" > "$CUSTOM_LIB"
  echo "  Custom library: $CUSTOM_LIB"
else
  echo "  [WARN] REPBASE_LIB not set — custom library will use RepeatModeler output only."
  echo "         Set: export REPBASE_LIB=/path/to/RepBase.fa"
  cp "$FILTERED" "$CUSTOM_LIB"
fi
echo "  Library sequences: $(grep -c '^>' "$CUSTOM_LIB")"

# --- Step 4: RepeatMasker ---------------------------------------------------
conda deactivate 2>/dev/null || true
activate_conda env_repeatmasker

echo; echo "=== RepeatMasker v4.0.7 ==="
mkdir -p "$WORKDIR/repeatmasker_out"
RepeatMasker \
  -pa 24 \
  -e ncbi \
  -gccalc \
  -lib "$CUSTOM_LIB" \
  -dir "$WORKDIR/repeatmasker_out" \
  "$FINAL_GENOME_FASTA"

echo
echo "RepeatMasker summary (.tbl):"
cat "$WORKDIR/repeatmasker_out/"*.tbl 2>/dev/null | sed 's/^/  /'
echo
echo "Thesis Table S8 expected: ~69% repeat content"
echo "(M6: low value with generic library is known — custom library should approach 69%)"

mark_done "$WORKDIR"
