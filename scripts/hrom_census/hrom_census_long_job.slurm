#!/bin/bash
# HROM structural census: three idempotent long workers are normally ample for
# its 62.87M within-species pairs. Resubmitting the same array resumes safely.
#SBATCH --job-name=hrom-census
#SBATCH --array=0-2%3
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/databases/HROM/census/logs/census_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/databases/HROM/census/logs/census_%A_%a.err

set -euo pipefail
REPO=/lustre1/g/aos_shihuang/Syn2b-paper
WORKDIR=/lustre1/g/aos_shihuang/databases/HROM/census
SYN2B=/lustre1/g/aos_shihuang/Syn2b/target/release/syn2b
TMPDIR_LOCAL=/tmp/${SLURM_JOB_ID}
mkdir -p "$WORKDIR/logs" "$TMPDIR_LOCAL"
cd "$REPO"
python3 scripts/gtdb_census/census_worker.py \
    --workdir "$WORKDIR" \
    --tasks "$WORKDIR/tasks.jsonl" \
    --syn2b "$SYN2B" \
    --ensure-filename-genome-id \
    --tmpdir "$TMPDIR_LOCAL" \
    --workers "$SLURM_CPUS_PER_TASK" \
    --max-mega 12 \
    --soft-gb 60 \
    --hard-gb 80 \
    --min-free-gb 50 \
    --drain
