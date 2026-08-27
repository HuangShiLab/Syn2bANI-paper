#!/bin/bash -l
#SBATCH --job-name=s2b_extract_og
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/s2b_extract_og_%j.out
#SBATCH --error=logs/s2b_extract_og_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

SYN2BANI_DIR=/lustre1/g/aos_shihuang/Syn2bANI
PAPER_DIR=/lustre1/g/aos_shihuang/Syn2bANI-paper
VAL_DIR=/lustre1/g/aos_shihuang/data/validation_oral_gut
GENOMES_DIR="${VAL_DIR}/genomes"
PAIRS_FILE="${VAL_DIR}/validation_pairs.tsv"
OUTPUT_DIR="${VAL_DIR}/v6_features"

mkdir -p "${OUTPUT_DIR}" logs

echo "=== Build Syn2bANI release ==="
cd "${SYN2BANI_DIR}"
cargo build --release

echo ""
echo "=== Extract v6 features for oral/gut validation pairs ==="
cd "${PAPER_DIR}"
python3 scripts/extract_syn2bani_features.py \
  --pairs "${PAIRS_FILE}" \
  --genomes "${GENOMES_DIR}" \
  --syn2bani "${SYN2BANI_DIR}/target/release/syn2bani" \
  --output "${OUTPUT_DIR}/oral_gut_v6_features.tsv" \
  --threads 8

echo ""
echo "=== Done ==="
ls -lh "${OUTPUT_DIR}"
