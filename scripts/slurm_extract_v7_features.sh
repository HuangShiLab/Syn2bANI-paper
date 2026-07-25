#!/bin/bash -l
#SBATCH --job-name=s2b_v7_features
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --array=1-25
#SBATCH --output=logs/s2b_v7_features_%A_%a.out
#SBATCH --error=logs/s2b_v7_features_%A_%a.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

SYN2BANI_DIR=/lustre1/g/aos_shihuang/Syn2bANI
GTDB_DIR=/lustre1/g/aos_shihuang/databases/GTDB/GTDBr207
GENOMES_DIR="${GTDB_DIR}/genomes"
PAIR_CHUNK="${GTDB_DIR}/pair_chunks/pairs_chunk_${SLURM_ARRAY_TASK_ID}.tsv"
OUTPUT_DIR="${GTDB_DIR}/v7_features"
OUTPUT_FILE="${OUTPUT_DIR}/features_chunk_${SLURM_ARRAY_TASK_ID}.tsv"

mkdir -p "${OUTPUT_DIR}" logs

if [[ ! -f "${PAIR_CHUNK}" ]]; then
    echo "Pair chunk not found: ${PAIR_CHUNK}"
    exit 0
fi

echo "=== Extracting v7 features for chunk ${SLURM_ARRAY_TASK_ID} ==="
python3 /lustre1/g/aos_shihuang/Syn2bANI-paper/scripts/extract_syn2bani_features.py \
  --pairs "${PAIR_CHUNK}" \
  --genomes "${GENOMES_DIR}" \
  --syn2bani "${SYN2BANI_DIR}/target/release/syn2bani" \
  --output "${OUTPUT_FILE}" \
  --threads 8

echo "=== Done: ${OUTPUT_FILE} ==="
