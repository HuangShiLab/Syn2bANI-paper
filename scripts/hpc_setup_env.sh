#!/bin/bash -l
#SBATCH --job-name=s2b_setup_env
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/setup_env_%j.out
#SBATCH --error=logs/setup_env_%j.err

set -euo pipefail

JOB_ID="${SLURM_JOB_ID:-login}"
LOG="logs/setup_env_${JOB_ID}.log"

echo "Setting up Syn2bANI environment on $(hostname) at $(date)" | tee -a "${LOG}"

# Python packages: install in user space using the base Python 3.12
python3 -m pip install --user numpy pandas scikit-learn matplotlib seaborn 2>&1 | tee -a "${LOG}"

# Install skani and fastANI in a minimal conda environment
ENV_NAME=syn2bani
conda remove -y -n "${ENV_NAME}" --all 2>/dev/null || true
conda create -y -n "${ENV_NAME}" python=3.11 2>&1 | tee -a "${LOG}"
conda install -y -n "${ENV_NAME}" -c bioconda -c conda-forge skani fastani 2>&1 | tee -a "${LOG}"

# Verify
echo "=== Verification ===" | tee -a "${LOG}"
python3 - <<'PY' 2>&1 | tee -a "${LOG}"
import numpy, pandas, sklearn, matplotlib, seaborn
print('numpy', numpy.__version__)
print('pandas', pandas.__version__)
print('sklearn', sklearn.__version__)
print('matplotlib', matplotlib.__version__)
print('seaborn', seaborn.__version__)
PY
/home/shihuang/.conda/envs/syn2bani/bin/skani --version 2>&1 | tee -a "${LOG}"
/home/shihuang/.conda/envs/syn2bani/bin/fastANI -v 2>&1 | tee -a "${LOG}"

echo "Setup complete at $(date)" | tee -a "${LOG}"
