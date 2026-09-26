#!/bin/bash
#SBATCH --job-name=syn2b-census
# FEW LONG JOBS policy: a small array of long-running nodes; tasks are claimed
# atomically from the shared workdir, so jobs may start late, die, or overlap
# without any lost or duplicated work. Resubmit the same array to resume.
#SBATCH --array=0-5%6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4-00:00:00
#SBATCH --mem=180G
#SBATCH --output=logs/census_%A_%a.out
#SBATCH --error=logs/census_%A_%a.err

# --- site-specific paths ----------------------------------------------------
WORKDIR=/lustre1/g/aos_shihuang/data/syn2b_census
SYN2B=/lustre1/g/aos_shihuang/Syn2b/target/release/syn2b
TMPDIR_LOCAL=/tmp/$SLURM_JOB_ID # fast local scratch for task directories
# Space audit thresholds (GB, on the shared workdir filesystem):
SOFT_GB=250
HARD_GB=320
MIN_FREE_GB=50
# ----------------------------------------------------------------------------

mkdir -p "$WORKDIR/logs" "$TMPDIR_LOCAL"
python3 scripts/gtdb_census/census_worker.py \
    --workdir "$WORKDIR" \
    --tasks "$WORKDIR/tasks.jsonl" \
    --syn2b "$SYN2B" \
    --tmpdir "$TMPDIR_LOCAL" \
    --workers "$SLURM_CPUS_PER_TASK" --max-mega 12 \
    --soft-gb "$SOFT_GB" --hard-gb "$HARD_GB" --min-free-gb "$MIN_FREE_GB" \
    --drain
