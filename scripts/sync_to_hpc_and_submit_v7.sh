#!/bin/bash -l
# Sync Syn2bANI code + paper scripts to HPC2021 and submit v7 mid-ANI benchmark.
# Run this from the Mac Studio (or any host that can ssh to hpc2021.hku.hk).

set -euo pipefail

HPC_USER=shihuang
HPC_HOST=hpc2021.hku.hk
HPC_WORK=/lustre1/g/aos_shihuang
LOCAL_SYN2BANI="${LOCAL_SYN2BANI:-$HOME/Downloads/Syn2bANI}"
LOCAL_PAPER="${LOCAL_PAPER:-$HOME/Downloads/Syn2bANI-paper}"

echo "=== Syncing Syn2bANI code to HPC ==="
rsync -avz --exclude='target' --exclude='.git' \
  "${LOCAL_SYN2BANI}/" "${HPC_USER}@${HPC_HOST}:${HPC_WORK}/Syn2bANI/"

echo ""
echo "=== Syncing Syn2bANI-paper scripts to HPC ==="
rsync -avz --exclude='.git' \
  "${LOCAL_PAPER}/scripts/" "${HPC_USER}@${HPC_HOST}:${HPC_WORK}/Syn2bANI-paper/scripts/"

echo ""
echo "=== Submitting v7 benchmark job ==="
ssh "${HPC_USER}@${HPC_HOST}" "bash -l -c 'cd ${HPC_WORK}/Syn2bANI-paper && mkdir -p logs && sbatch scripts/slurm_v7_mid_ani.sh'"

echo ""
echo "=== Done. Check queue with: squeue -u ${HPC_USER} ==="
