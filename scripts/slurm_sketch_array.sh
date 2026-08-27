#!/bin/bash -l
#SBATCH --job-name=s2b_sketch
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-1000%50
#SBATCH --output=logs/sketch_%A_%a.out
#SBATCH --error=logs/sketch_%A_%a.err

set -euo pipefail

# Activate conda environment
source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

# Configuration
GENOMES_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
SKETCH_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/sketches
MANIFEST=/lustre1/g/aos_shihuang/data/gtdb-r207/metadata/genome_manifest.txt
SYN2BANI=/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani
CHUNK_SIZE=65

mkdir -p "${SKETCH_DIR}"

START=$(((SLURM_ARRAY_TASK_ID - 1) * CHUNK_SIZE + 1))
END=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE))

if [ ! -f "${MANIFEST}" ]; then
    echo "Manifest not found: ${MANIFEST}"
    exit 1
fi

sed -n "${START},${END}p" "${MANIFEST}" | while read -r ACC PATH; do
    [ -z "${ACC}" ] && continue
    if [ ! -f "${PATH}" ]; then
        echo "Genome not found: ${PATH} (skipping)"
        continue
    fi
    "${SYN2BANI}" sketch "${PATH}" -o "${SKETCH_DIR}/${SLURM_ARRAY_TASK_ID}/" -e BcgI -t 8 -p || true
done

echo "Sketch task ${SLURM_ARRAY_TASK_ID} complete"
