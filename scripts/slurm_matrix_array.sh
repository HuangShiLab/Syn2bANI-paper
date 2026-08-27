#!/bin/bash -l
#SBATCH --job-name=s2b_matrix
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --array=1-25
#SBATCH --output=logs/matrix_%A_%a.out
#SBATCH --error=logs/matrix_%A_%a.err

set -euo pipefail

# Activate conda environment with tools
source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

# Configuration
GENOMES_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
PAIRS_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/pair_chunks
OUT_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/matrix_chunks
SYN2BANI=/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani
SKANI=/home/shihuang/.conda/envs/syn2bani/bin/skani

mkdir -p "${OUT_DIR}"

CHUNK_FILE="${PAIRS_DIR}/pairs_chunk_${SLURM_ARRAY_TASK_ID}.tsv"
OUTPUT_FILE="${OUT_DIR}/chunk_${SLURM_ARRAY_TASK_ID}.tsv"

if [ ! -f "${CHUNK_FILE}" ]; then
    echo "Chunk file not found: ${CHUNK_FILE} (skipping)"
    exit 0
fi

# Phase 1: run Syn2bANI + skani on all pairs. FastANI is run separately on a
# stratified subset in phase 2 to save CPU time.
/usr/bin/time -v python3 scripts/run_benchmark_chunk.py \
  --pairs "${CHUNK_FILE}" \
  --genomes "${GENOMES_DIR}" \
  --syn2bani "${SYN2BANI}" \
  --skani "${SKANI}" \
  --output "${OUTPUT_FILE}" \
  --threads 8 \
  --tools syn2bani,skani \
  --chunk-size 1000

echo "Matrix task ${SLURM_ARRAY_TASK_ID} complete"
