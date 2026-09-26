#!/bin/bash
# Idempotent within-cluster skani ANI pass. This can run before, with, or after
# the structural array.
#SBATCH --job-name=hrom-ani
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/databases/HROM/census/logs/ani_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/databases/HROM/census/logs/ani_%j.err

set -euo pipefail
REPO=/lustre1/g/aos_shihuang/Syn2b-paper
WORKDIR=/lustre1/g/aos_shihuang/databases/HROM/census
SKANI=/lustre1/g/aos_shihuang/tools/bin/skani
mkdir -p "$WORKDIR/logs"
cd "$REPO"
python3 scripts/gtdb_census/skani_ani_pass.py \
    --workdir "$WORKDIR" \
    --skani "$SKANI" \
    --workers 16
