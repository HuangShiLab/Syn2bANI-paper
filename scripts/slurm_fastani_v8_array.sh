#!/bin/bash -l
#SBATCH --job-name=s2b_v8_fastani
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=1-5
#SBATCH --output=logs/fastani_v8_%A_%a.out
#SBATCH --error=logs/fastani_v8_%A_%a.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

GENOMES_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
PAIRS_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/fastani_v8_pair_chunks
OUT_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/fastani_v8_chunks
FASTANI=/home/shihuang/.conda/envs/syn2bani/bin/fastANI

mkdir -p "${OUT_DIR}"

CHUNK_FILE="${PAIRS_DIR}/pairs_chunk_${SLURM_ARRAY_TASK_ID}.tsv"
OUTPUT_FILE="${OUT_DIR}/chunk_${SLURM_ARRAY_TASK_ID}.tsv"

if [ ! -f "${CHUNK_FILE}" ]; then
    echo "Chunk file not found: ${CHUNK_FILE} (skipping)"
    exit 0
fi

/usr/bin/time -v python3 scripts/run_benchmark_chunk.py \
  --pairs "${CHUNK_FILE}" \
  --genomes "${GENOMES_DIR}" \
  --fastani "${FASTANI}" \
  --output "${OUTPUT_FILE}" \
  --threads 4 \
  --tools fastani \
  --chunk-size 1000

echo "FastANI v8 task ${SLURM_ARRAY_TASK_ID} complete"
