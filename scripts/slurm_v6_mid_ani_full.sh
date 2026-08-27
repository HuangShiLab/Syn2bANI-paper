#!/bin/bash -l
#SBATCH --job-name=s2b_v6_full
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/s2b_v6_full_%j.out
#SBATCH --error=logs/s2b_v6_full_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

SYN2BANI_DIR=/lustre1/g/aos_shihuang/Syn2bANI
PAPER_DIR=/lustre1/g/aos_shihuang/Syn2bANI-paper
VAL_DIR=/lustre1/g/aos_shihuang/data/validation_mid_ani
GENOMES_DIR="${VAL_DIR}/genomes"
PAIRS_FILE="${VAL_DIR}/mid_ani_pairs.tsv"
FASTANI_MATRIX="${VAL_DIR}/mid_ani_matrix_fastani.tsv"
SKANI=/lustre1/g/aos_shihuang/tools/anaconda3/pkgs/skani-0.2.2-ha6fb395_2/bin/skani
OUTPUT_DIR="${PAPER_DIR}/results/validation/v6_mid_ani_full"

mkdir -p "${OUTPUT_DIR}" logs

echo "=== Build Syn2bANI release ==="
cd "${SYN2BANI_DIR}"
cargo build --release

echo ""
echo "=== Benchmark v6 on full within-genus validation set ==="
cd "${PAPER_DIR}"
python3 scripts/benchmark_v6_mid_ani.py \
  --pairs "${PAIRS_FILE}" \
  --genomes "${GENOMES_DIR}" \
  --syn2bani "${SYN2BANI_DIR}/target/release/syn2bani" \
  --skani "${SKANI}" \
  --fastani-matrix "${FASTANI_MATRIX}" \
  --output-dir "${OUTPUT_DIR}" \
  --threads 8

echo ""
echo "=== Validation complete ==="
ls -lh "${OUTPUT_DIR}"
