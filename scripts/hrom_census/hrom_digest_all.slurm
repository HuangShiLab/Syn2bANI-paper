#!/bin/bash
# One-time HROM pre-digest. 20 shards fit Lustre metadata and syn2b I/O overlap
# within the current account submit limit.
#SBATCH --job-name=hrom-digest
#SBATCH --array=0-19%20
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/databases/HROM/census/logs/digest_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/databases/HROM/census/logs/digest_%A_%a.err

set -euo pipefail
REPO=/lustre1/g/aos_shihuang/Syn2b-paper
WORKDIR=/lustre1/g/aos_shihuang/databases/HROM/census
SYN2B=/lustre1/g/aos_shihuang/Syn2b/target/release/syn2b
mkdir -p "$WORKDIR/logs"
cd "$REPO"
python3 scripts/gtdb_census/digest_all.py \
    --workdir "$WORKDIR" \
    --syn2b "$SYN2B" \
    --shard "$SLURM_ARRAY_TASK_ID" \
    --nshards 20 \
    --ensure-filename-genome-id
