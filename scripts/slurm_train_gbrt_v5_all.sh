#!/bin/bash -l
#SBATCH --job-name=s2b_train_v5_all
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=3:00:00
#SBATCH --output=logs/train_gbrt_v5_all_%j.out
#SBATCH --error=logs/train_gbrt_v5_all_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

TRAIN_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207

echo "=== Train GBRT v5 BcgI ==="
python3 scripts/train_gbrt_v5_combined.py \
  --matrix "${TRAIN_DIR}/train_bcgi.tsv" \
  --mode bcgi \
  --output results/gbrt_model_v5_bcgi.json \
  --report results/gbrt_v5_bcgi_report.txt

echo ""
echo "=== Train GBRT v5 CjePI ==="
python3 scripts/train_gbrt_v5_combined.py \
  --matrix "${TRAIN_DIR}/train_cjepi.tsv" \
  --mode cjepi \
  --output results/gbrt_model_v5_cjepi.json \
  --report results/gbrt_v5_cjepi_report.txt

echo ""
echo "=== Train GBRT v5 Combined ==="
python3 scripts/train_gbrt_v5_combined.py \
  --matrix "${TRAIN_DIR}/train_combined.tsv" \
  --mode combined \
  --output results/gbrt_model_v5_combined.json \
  --report results/gbrt_v5_combined_report.txt

echo ""
echo "=== All trainings complete ==="
ls -lh results/gbrt_model_v5_*.json results/gbrt_model_v5_*.pkl
