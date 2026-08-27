#!/bin/bash -l
#SBATCH --job-name=s2b_train_v5
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --output=logs/train_gbrt_v5_bcgI_%j.out
#SBATCH --error=logs/train_gbrt_v5_bcgI_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

python3 scripts/train_gbrt_v5.py \
  --matrix results/matrix_gtdb_r207_100k.tsv \
  --output results/gbrt_model_v5_bcgI.json \
  --report results/gbrt_v5_bcgI_report.txt

echo "=== Training complete ==="
ls -lh results/gbrt_model_v5_bcgI.json
ls -lh results/gbrt_model_v5_bcgI.pkl
