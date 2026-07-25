#!/bin/bash -l
# Submit v7 feature extraction array job + training job on HPC2021.
# Run this on a login node after syncing the latest Syn2bANI code.

set -euo pipefail

cd /lustre1/g/aos_shihuang/Syn2bANI-paper
mkdir -p logs

echo "=== Building Syn2bANI release ==="
cd /lustre1/g/aos_shihuang/Syn2bANI
cargo build --release

echo ""
echo "=== Submitting v7 feature extraction array job ==="
cd /lustre1/g/aos_shihuang/Syn2bANI-paper
FEATURES_JOBID=$(sbatch --parsable scripts/slurm_extract_v7_features.sh)
echo "Features job ID: ${FEATURES_JOBID}"

echo ""
echo "=== Submitting v7 training job (depends on features) ==="
TRAIN_JOBID=$(sbatch --parsable --dependency=afterany:${FEATURES_JOBID} scripts/slurm_train_v7_model.sh)
echo "Training job ID: ${TRAIN_JOBID}"

echo ""
echo "=== Submitted ==="
echo "Monitor with: squeue -u shihuang"
echo "Features logs: logs/s2b_v7_features_${FEATURES_JOBID}_*.out"
echo "Training log:  logs/s2b_v7_train_${TRAIN_JOBID}.out"
