#!/usr/bin/env bash
# Install skani into the 'syn2bani' conda environment.
# Run on the login or I/O node.
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate syn2bani

if command -v skani >/dev/null 2>&1; then
  echo "skani already installed: $(skani --version)"
  exit 0
fi

echo "Installing skani into the syn2bani conda environment..."
conda install -c bioconda skani -y

echo "Done: $(skani --version)"
