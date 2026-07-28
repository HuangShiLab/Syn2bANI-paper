#!/usr/bin/env bash
# Submit the full Syn2bANI v8 HPC benchmark workflow with SLURM dependencies.
#
# Usage:
#   cd /lustre1/g/aos_shihuang/Syn2bANI-paper
#   bash scripts/submit_hpc_workflow_v8.sh

set -euo pipefail

mkdir -p logs

echo "Submitting Syn2bANI v8 HPC workflow..."

MATRIX=$(sbatch scripts/slurm_matrix_v8_array.sh | awk '{print $4}')
echo "  Matrix v8 array job: $MATRIX"

MERGE=$(sbatch --dependency=afterok:"$MATRIX" scripts/slurm_merge_v8.sh | awk '{print $4}')
echo "  Merge+sample job: $MERGE"

FASTANI=$(sbatch --dependency=afterok:"$MERGE" scripts/slurm_fastani_v8_array.sh | awk '{print $4}')
echo "  FastANI array job: $FASTANI"

FINAL=$(sbatch --dependency=afterok:"$FASTANI" scripts/slurm_merge_v8_final.sh | awk '{print $4}')
echo "  Final merge job: $FINAL"

echo ""
echo "Workflow submitted. Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/matrix_v8_${MATRIX}_*.err logs/merge_v8_p1_${MERGE}.err logs/fastani_v8_${FASTANI}_*.err logs/merge_v8_final_${FINAL}.err"
