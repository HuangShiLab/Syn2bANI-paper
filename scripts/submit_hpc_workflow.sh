#!/bin/bash -l
# Submit the full Syn2bANI HPC benchmark workflow with SLURM dependencies.
#
# Usage:
#   cd /lustre1/g/aos_shihuang/Syn2bANI-paper
#   bash scripts/submit_hpc_workflow.sh

set -euo pipefail

mkdir -p logs

echo "Submitting Syn2bANI HPC workflow..."

MATRIX=$(sbatch scripts/slurm_matrix_array.sh | awk '{print $4}')
echo "  Matrix array job: $MATRIX"

MERGE=$(sbatch --dependency=afterok:"$MATRIX" scripts/slurm_merge_phase1.sh | awk '{print $4}')
echo "  Merge+sample job: $MERGE"

FASTANI=$(sbatch --dependency=afterok:"$MERGE" scripts/slurm_fastani_array.sh | awk '{print $4}')
echo "  FastANI array job: $FASTANI"

TRAIN=$(sbatch --dependency=afterok:"$FASTANI" scripts/slurm_train_final.sh | awk '{print $4}')
echo "  Final train job: $TRAIN"

echo ""
echo "Workflow submitted. Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/matrix_${MATRIX}_*.err logs/merge_p1_${MERGE}.err logs/fastani_${FASTANI}_*.err logs/train_${TRAIN}.err"
