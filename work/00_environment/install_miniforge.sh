#!/bin/bash
# Installs Miniforge (conda + mamba, conda-forge default channel) into WSL2 Ubuntu.
# Idempotent: skips download/install if already present.
set -euo pipefail

INSTALL_DIR="$HOME/miniforge3"

if [ -d "$INSTALL_DIR" ]; then
  echo "Miniforge already installed at $INSTALL_DIR"
else
  echo "Downloading Miniforge installer..."
  cd /tmp
  curl -fsSL -o Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  bash Miniforge3.sh -b -p "$INSTALL_DIR"
  rm -f Miniforge3.sh
fi

# Initialize for bash (idempotent — conda init checks for existing block)
"$INSTALL_DIR/bin/conda" init bash

echo
echo "Installed versions:"
"$INSTALL_DIR/bin/conda" --version
"$INSTALL_DIR/bin/mamba" --version
