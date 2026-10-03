#!/usr/bin/env bash
# ==============================================================================
# STEP 05 — Reference-guided scaffolding (RagTag v2.1.0)
#
# Thesis (Methods p.18, Appendix I p.67):
#   Individual 1A polished+purged assembly scaffolded against Tcas5.2
#   (GCF_000002335.3) using RagTag scaffold.
#   The -C flag causes unplaced contigs to be concatenated into a single
#   "unplaced_contigs" record in some intermediate outputs — downstream steps
#   use the version that keeps unplaced contigs as separate records, matching
#   the published scaffold count.
#   RagTag v2.1.0 confirmed via Figure S1.
#
# Conda env: env_assembly  (ragtag=2.1.0, minimap2)
#
# Usage:
#   bash 05_scaffolding_ragtag.sh           # check: list existing outputs
#   bash 05_scaffolding_ragtag.sh run       # run (skip if done)
#   bash 05_scaffolding_ragtag.sh run --force
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config/paths.config.sh"
source "$SCRIPT_DIR/_lib.sh"

RUN_MODE=$(mode "${1:-check}")
FORCE="${2:-}"
WORKDIR="$PIPELINE_ROOT/work/05_scaffolding_ragtag"

banner "STEP 05: RagTag reference-guided scaffolding [$RUN_MODE]"
require_paths "$RAGTAG_OUTPUT" "$PUBREF_TCAS52"

# ── CHECK mode ──────────────────────────────────────────────────────────────
if [[ "$RUN_MODE" == "check" ]]; then
  echo
  echo "Existing RagTag output (D:):"
  find "$RAGTAG_OUTPUT" -maxdepth 2 2>/dev/null | sort | sed 's/^/  /'
  exit 0
fi

# ── RUN mode ────────────────────────────────────────────────────────────────
guard_expensive "$WORKDIR" "$FORCE"
mkdir -p "$WORKDIR"

activate_conda env_assembly

# Locate the 1A polished+purged assembly from Step 04
DRAFT_1A=""
for candidate in \
  "$PIPELINE_ROOT/work/04_polishing/da01a/purge_haplotigs_r2/da01a_purged_r2.fasta" \
  "$PIPELINE_ROOT/work/04_polishing/01A/purge_haplotigs_r2/01A_purged_r2.fasta" \
  "$PIPELINE_ROOT/work/04_polishing/1A/purge_haplotigs_r2/1A_purged_r2.fasta"; do
  [[ -f "$candidate" ]] && { DRAFT_1A="$candidate"; break; }
done

if [[ -z "$DRAFT_1A" ]]; then
  # Fall back to the existing purge_haplotigs output on F:
  DRAFT_1A=$(find "$PURGE_HAPLOTIGS_OUTPUT" -name "*da01a*purged*.fasta" -o \
                                             -name "*1A*purged*.fasta" 2>/dev/null | head -1)
fi

[[ -z "$DRAFT_1A" ]] && { echo "[ERROR] 1A purged+polished assembly not found. Run Step 03/04 first." >&2; exit 1; }
echo "Input 1A assembly: $DRAFT_1A"

# Reference fasta
REF_FASTA=$(find "$PUBREF_TCAS52" -name "*.fna" -o -name "*genomic*.fna" 2>/dev/null | head -1)
[[ -z "$REF_FASTA" ]] && { echo "[ERROR] Tcas5.2 reference fna not found in $PUBREF_TCAS52" >&2; exit 1; }
echo "Reference: $REF_FASTA"

ragtag.py scaffold \
  "$REF_FASTA" \
  "$DRAFT_1A" \
  -o "$WORKDIR/ragtag_output" \
  -t 16 \
  -C

echo
echo "Scaffolded assembly: $WORKDIR/ragtag_output/ragtag.scaffold.fasta"
fasta_summary "$WORKDIR/ragtag_output/ragtag.scaffold.fasta"

mark_done "$WORKDIR"
