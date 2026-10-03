#!/bin/bash
# env_qc: FastQC, Trimmomatic, Jellyfish, KMC, GenomeScope2 — Steps 00/01
# Thesis versions: FastQC v0.11.8 (well, methods say "v.0.11.8" Appendix I), Trimmomatic v0.38
set -euo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"

ENV_NAME="env_qc"
if conda env list | grep -q "^$ENV_NAME "; then
  echo "$ENV_NAME already exists"
else
  # NOTE: "jellyfish" on conda-forge is an unrelated string-similarity library.
  # The genomic k-mer counter (Marcais & Kingsford 2011) is "kmer-jellyfish" on bioconda.
  mamba create -y -n "$ENV_NAME" -c bioconda -c conda-forge \
    fastqc=0.11.8 \
    trimmomatic=0.38 \
    kmer-jellyfish=2.3.0 \
    kmc=3.2.1 \
    genomescope2=2.0
fi

conda activate "$ENV_NAME"
echo
echo "=== Installed versions (env_qc) ==="
fastqc --version
trimmomatic --version 2>&1 | head -1
jellyfish --version
kmc -h 2>&1 | head -3
genomescope2 --version 2>&1 | head -1 || conda list -n "$ENV_NAME" | grep genomescope2
