#!/bin/bash -l
#SBATCH --job-name=s2b_enzyme_cmp
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/enzyme_comparison_%j.out
#SBATCH --error=logs/enzyme_comparison_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

VAL_DIR=/lustre1/g/aos_shihuang/data/validation_mid_ani
OUT_FILE="${VAL_DIR}/enzyme_comparison.tsv"

python3 scripts/run_enzyme_comparison.py \
  --pairs "${VAL_DIR}/mid_ani_pairs_85_95.tsv" \
  --output "${OUT_FILE}"

echo "=== Enzyme comparison complete ==="
ls -lh "${OUT_FILE}"
