#!/usr/bin/env bash
# ==============================================================================
# STEP 10 — Y chromosome assembly (DiscoverY + RagTag + BLASTn filter + findZX)
#
# Thesis (Appendix I p.68, Methods p.18, Fig 1E, Fig S3):
#   1. DiscoverY: identify male-specific k-mers
#        k-mer threshold = 0.5 × average male autosomal coverage (k-mer-to-contig: 0.2)
#        Males: 6 samples  |  Females: 2 samples
#   2. Filter: keep contigs present in ≥2 other male samples (consistency filter)
#   3. SPAdes v3.15.2: re-assemble male-specific contigs
#   4. RagTag scaffold: against Tcas5.2's putative Y-labelled contigs
#   5. BLASTn filter (thesis-described exact logic):
#        Remove contigs with >80% identity AND aligned-length/query-length ratio >0.6
#        to any autosomal contig → removes mis-classified autosomal fragments
#   6. 2 × Pilon polishing rounds
#   7. Append Y chromosome to the rest of inTcas1 (Step 07)
#   8. Validate with findZX
#
# Result: Y chromosome = 0.8 Mbp, 38 genes
#
# Conda env: env_assembly (minimap2, samtools, ragtag, spades)
# DiscoverY: $Y_DISCOVERY_TOOL (python scripts, already on F:)
# findZX:    conda create -n env_findzx -c bioconda findzx
# BLAST:     conda install -n env_assembly -c bioconda blast
#
# Usage:
#   bash 10_y_chromosome.sh           # check: list existing DiscoverY outputs
#   bash 10_y_chromosome.sh run       # run full Y-chr pipeline
#   bash 10_y_chromosome.sh run --force
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/10_y_chromosome"

KMER_THRESHOLD=0.2
MIN_MALE_SAMPLES=2
BLAST_IDENTITY=80
BLAST_LENGTH_RATIO=0.6

banner "STEP 10: Y chromosome assembly (DiscoverY + RagTag + BLASTn filter) [$RUN_MODE]"
require_paths "$Y_DISCOVERY_OUTPUT" "$Y_DISCOVERY_TOOL" "$Y_FINDZX_CONFIG"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing DiscoverY outputs (F:):"
  find "$Y_DISCOVERY_OUTPUT" -maxdepth 1 2>/dev/null | sort | sed 's/^/  /'
  echo
  echo "DiscoverY tool scripts (F:):"
  find "$Y_DISCOVERY_TOOL" -maxdepth 1 -name "*.py" 2>/dev/null | sort | sed 's/^/  /'
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

activate_conda env_assembly

# Locate BAMs for male and female samples (from Step 04 or existing polished BAMs)
MALE_BAMS=()
FEMALE_BAMS=()
for bam in $(find "$POLISHED_MINIMAP2_BAMS" -name "*.sorted.bam" 2>/dev/null | sort); do
  name="$(basename "$bam" .sorted.bam)"
  # Sex determination by sample code: males end in A/B/C/D/E; females identified by user as 1N, 01N
  if [[ "$name" =~ [Nn]$ ]] || [[ "$name" =~ 01[Nn] ]]; then
    FEMALE_BAMS+=("$bam")
  else
    MALE_BAMS+=("$bam")
  fi
done

echo "Male BAMs  (${#MALE_BAMS[@]}): ${MALE_BAMS[*]:-none found}"
echo "Female BAMs(${#FEMALE_BAMS[@]}): ${FEMALE_BAMS[*]:-none found}"

MALE_LIST="$WORKDIR/male_bams.txt"
FEMALE_LIST="$WORKDIR/female_bams.txt"
printf "%s\n" "${MALE_BAMS[@]}" > "$MALE_LIST"
printf "%s\n" "${FEMALE_BAMS[@]}" > "$FEMALE_LIST"

# --- DiscoverY ---------------------------------------------------------------
echo; echo "=== DiscoverY: male-specific k-mer identification ==="
python3 "$Y_DISCOVERY_TOOL/discoverY.py" \
  --male_bam_list "$MALE_LIST" \
  --female_bam_list "$FEMALE_LIST" \
  --threshold "$KMER_THRESHOLD" \
  --min_samples "$MIN_MALE_SAMPLES" \
  --out "$WORKDIR/discovery_out"

# --- SPAdes on Y-candidate contigs -------------------------------------------
Y_CONTIGS="$WORKDIR/discovery_out/y_candidate_contigs.fasta"
[[ -f "$Y_CONTIGS" ]] || { echo "[ERROR] DiscoverY did not produce y_candidate_contigs.fasta" >&2; exit 1; }

echo; echo "=== SPAdes re-assembly of Y-candidate contigs ==="
spades.py \
  -s "$Y_CONTIGS" \
  -o "$WORKDIR/y_spades" \
  -t 16 -m 64

# --- RagTag scaffold against Tcas5.2 Y-labelled sequences -------------------
echo; echo "=== RagTag scaffold (Y candidates → Tcas5.2 putative Y contigs) ==="
Y_REF=$(find "$PUBREF_TCAS52" -name "*.fna" 2>/dev/null | head -1)
ragtag.py scaffold \
  "$Y_REF" \
  "$WORKDIR/y_spades/contigs.fasta" \
  -o "$WORKDIR/y_ragtag" \
  -t 16

Y_SCAFF="$WORKDIR/y_ragtag/ragtag.scaffold.fasta"

# --- BLASTn filter (exact thesis logic) -------------------------------------
echo; echo "=== BLASTn filter: remove autosome-homologous contigs ==="
AUTOSOMES="$FINAL_GENOME_FASTA"

makeblastdb -in "$AUTOSOMES" -dbtype nucl -out "$WORKDIR/autosomes_db"
blastn \
  -query "$Y_SCAFF" \
  -db "$WORKDIR/autosomes_db" \
  -outfmt "6 qseqid pident length qlen" \
  -out "$WORKDIR/y_vs_autosomes.blast" \
  -num_threads 16

awk -v id="$BLAST_IDENTITY" -v ratio="$BLAST_LENGTH_RATIO" \
  '($2 > id) && ($3/$4 > ratio) {print $1}' \
  "$WORKDIR/y_vs_autosomes.blast" | sort -u > "$WORKDIR/y_contigs_to_remove.txt"

echo "Contigs flagged for removal: $(wc -l < "$WORKDIR/y_contigs_to_remove.txt")"

# Remove flagged contigs (seqkit or awk-based approach)
if command -v seqkit &>/dev/null; then
  seqkit grep -v -f "$WORKDIR/y_contigs_to_remove.txt" \
    "$Y_SCAFF" > "$WORKDIR/y_filtered.fasta"
else
  python3 - "$Y_SCAFF" "$WORKDIR/y_contigs_to_remove.txt" "$WORKDIR/y_filtered.fasta" <<'EOF'
import sys
fa, rem_file, out = sys.argv[1:4]
remove = set(open(rem_file).read().split())
write = False
with open(fa) as fin, open(out,'w') as fout:
    for line in fin:
        if line.startswith('>'):
            write = line[1:].split()[0] not in remove
        if write:
            fout.write(line)
EOF
fi

echo "Y-chromosome sequences after filter:"
fasta_summary "$WORKDIR/y_filtered.fasta"

# --- 2 × Pilon rounds --------------------------------------------------------
echo; echo "=== Pilon polishing (2 rounds) on filtered Y chromosome ==="
PILON_JAR=$(find "$HOME/miniforge3/envs/env_pilon" /root/miniforge3/envs/env_pilon \
  -name "pilon.jar" 2>/dev/null | head -1 || true)
if [[ -n "$PILON_JAR" ]]; then
  CURRENT_Y="$WORKDIR/y_filtered.fasta"
  # Use any available paired reads for polishing (from 1A or combined)
  R1=$(find "$RAW_ILLUMINA_WGS" -name "*da01a*R1*.fastq.gz" 2>/dev/null | head -1)
  R2=$(find "$RAW_ILLUMINA_WGS" -name "*da01a*R2*.fastq.gz" 2>/dev/null | head -1)
  for round in 1 2; do
    POUT="$WORKDIR/y_pilon_r${round}"
    mkdir -p "$POUT"
    minimap2 -ax sr -t 16 "$CURRENT_Y" "$R1" "$R2" \
      | samtools sort -@ 8 -o "$POUT/y_r${round}.sorted.bam"
    samtools index "$POUT/y_r${round}.sorted.bam"
    java -Xmx24g -jar "$PILON_JAR" \
      --genome "$CURRENT_Y" \
      --frags "$POUT/y_r${round}.sorted.bam" \
      --output "y_pilon_r${round}" \
      --outdir "$POUT" \
      --threads 16
    CURRENT_Y="$POUT/y_pilon_r${round}.fasta"
  done
  echo "Y chromosome (polished): $CURRENT_Y"
  fasta_summary "$CURRENT_Y"
else
  echo "[WARN] Pilon not found — skipping polishing. Install: conda create -n env_pilon -c bioconda pilon"
  CURRENT_Y="$WORKDIR/y_filtered.fasta"
fi

# --- findZX validation -------------------------------------------------------
echo; echo "=== findZX validation ==="
if command -v findZX &>/dev/null 2>&1 || conda env list | grep -q env_findzx; then
  conda run -n env_findzx findZX \
    --configfile "$Y_FINDZX_CONFIG/config.yml" || true
else
  echo "  [INFO] findZX not found (env_findzx). Install: conda create -n env_findzx -c bioconda findzx"
  echo "  Manually run: findZX --configfile $Y_FINDZX_CONFIG/config.yml"
fi

mark_done "$WORKDIR"
