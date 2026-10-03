#!/usr/bin/env bash
# Shared helpers sourced by every numbered step script. Not run directly.

banner() {
  echo "================================================================"
  echo "$1"
  echo "================================================================"
}

require_paths() {
  local missing=0
  for p in "$@"; do
    if [ ! -e "$p" ]; then
      echo "  [MISSING] $p" >&2
      missing=1
    else
      echo "  [OK]      $p" >&2
    fi
  done
  return $missing
}

mode() {
  if [ "${1:-check}" == "run" ]; then echo "run"; else echo "check"; fi
}

# Activate a conda environment in a non-interactive (script) shell
activate_conda() {
  local env_name="$1"
  local conda_sh=""
  for candidate in \
    "$HOME/miniforge3/etc/profile.d/conda.sh" \
    "/root/miniforge3/etc/profile.d/conda.sh" \
    "/opt/conda/etc/profile.d/conda.sh" \
    "$HOME/miniconda3/etc/profile.d/conda.sh"; do
    if [ -f "$candidate" ]; then conda_sh="$candidate"; break; fi
  done
  if [ -z "$conda_sh" ]; then
    echo "[ERROR] conda not found. Install Miniforge3 first (see work/00_environment/)." >&2
    exit 1
  fi
  # shellcheck disable=SC1090
  source "$conda_sh"
  conda activate "$env_name" 2>/dev/null || {
    echo "[ERROR] conda env '$env_name' not found. Check work/00_environment/ install scripts." >&2
    exit 1
  }
}

# Run a command inside a conda env without polluting the current shell
conda_run() {
  local env_name="$1"; shift
  if command -v conda &>/dev/null; then
    conda run --no-capture-output -n "$env_name" "$@"
  else
    activate_conda "$env_name"
    "$@"
  fi
}

# Return true (0) if this step already completed (done sentinel exists)
step_done() {
  [[ -f "${1}/.done" ]]
}

# Write the sentinel that marks a step complete
mark_done() {
  touch "${1}/.done"
  echo "[DONE] Step complete. Outputs in: $1"
}

# Guard for expensive steps: require --force to re-run completed work
guard_expensive() {
  local workdir="$1"
  local force="${2:-}"
  if step_done "$workdir" && [[ "$force" != "--force" ]]; then
    echo "[SKIP] Already done: ${workdir}/.done"
    echo "       Re-running is expensive. Add '--force' as 3rd arg to override."
    exit 0
  fi
}

# Print a quick summary of a FASTA file using only python3 (no .fai needed)
fasta_summary() {
  local fa="$1"
  [ -f "$fa" ] || { echo "  (not found: $fa)"; return; }
  python3 - "$fa" <<'EOF'
import sys, gzip
f = sys.argv[1]
opener = gzip.open if f.endswith('.gz') else open
nseq = tot = ns = 0
with opener(f,'rt') as fh:
    for line in fh:
        if line.startswith('>'):
            nseq += 1
        else:
            s = line.rstrip()
            tot += len(s)
            ns  += s.count('N') + s.count('n')
print(f"  sequences : {nseq:,}")
print(f"  total bp  : {tot:,}")
print(f"  Ns        : {ns:,}")
print(f"  ungapped  : {tot-ns:,}")
EOF
}
