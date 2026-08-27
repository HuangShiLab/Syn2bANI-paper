#!/bin/bash -l
#SBATCH --job-name=s2b_validate
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/validate_oral_gut_%j.out
#SBATCH --error=logs/validate_oral_gut_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

VAL_DIR=/lustre1/g/aos_shihuang/data/validation_oral_gut
GENOMES_DIR="${VAL_DIR}/genomes"
PAIRS_FILE="${VAL_DIR}/validation_pairs.tsv"
OUT_FILE="${VAL_DIR}/validation_matrix.tsv"
SYN2BANI=/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani
SKANI=/home/shihuang/.conda/envs/syn2bani/bin/skani
FASTANI=/home/shihuang/.conda/envs/syn2bani/bin/fastANI

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

echo "=== Build validation pairs ==="
python3 scripts/build_validation_pairs.py \
  --manifest "${VAL_DIR}/manifest.tsv" \
  --output "${PAIRS_FILE}"

echo "=== Run benchmark (FastANI + skani + Syn2bANI) ==="
/usr/bin/time -v python3 scripts/run_benchmark_matrix_v2.py \
  --pairs "${PAIRS_FILE}" \
  --genomes "${GENOMES_DIR}" \
  --syn2bani "${SYN2BANI}" \
  --skani "${SKANI}" \
  --fastani "${FASTANI}" \
  --output "${OUT_FILE}" \
  --threads 8 \
  --tools all \
  --chunk-size 500

echo "=== Validation complete ==="
ls -lh "${OUT_FILE}"
