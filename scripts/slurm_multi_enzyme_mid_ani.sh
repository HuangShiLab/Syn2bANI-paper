#!/bin/bash -l
#SBATCH --job-name=s2b_multi_mid
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=logs/multi_enzyme_mid_ani_%j.out
#SBATCH --error=logs/multi_enzyme_mid_ani_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

VAL_DIR=/lustre1/g/aos_shihuang/data/validation_mid_ani
OUT_FILE="${VAL_DIR}/multi_enzyme_mid_ani.tsv"

python3 scripts/run_multi_enzyme_mid_ani.py \
  --pairs "${VAL_DIR}/mid_ani_pairs_85_95.tsv" \
  --output "${OUT_FILE}"

echo "=== Multi-enzyme mid-ANI complete ==="
ls -lh "${OUT_FILE}"
