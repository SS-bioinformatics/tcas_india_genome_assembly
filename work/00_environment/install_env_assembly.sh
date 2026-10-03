#!/bin/bash
# env_assembly: purge_haplotigs, ragtag, samtools, minimap2, bedtools — for the Step 05/06 re-run
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"

ENV_NAME="env_assembly"
if conda env list | grep -q "^$ENV_NAME "; then
  echo "$ENV_NAME already exists"
else
  # purge_haplotigs: thesis pins v1.0.0; closest available on bioconda is 1.0.1 (logged as a
  # version-drift note per Phase 0 discipline). ragtag: thesis doesn't pin a version, using latest stable.
  mamba create -y -n "$ENV_NAME" -c bioconda -c conda-forge \
    purge_haplotigs=1.0.1 \
    ragtag=2.1.0 \
    samtools \
    minimap2 \
    bedtools \
    r-base
fi

conda activate "$ENV_NAME"
echo
echo "=== Installed versions (env_assembly) ==="
purge_haplotigs --version 2>&1 | head -1 || conda list -n "$ENV_NAME" | grep purge_haplotigs
ragtag.py --version 2>&1 | head -1
samtools --version 2>&1 | head -1
minimap2 --version 2>&1 | head -1
