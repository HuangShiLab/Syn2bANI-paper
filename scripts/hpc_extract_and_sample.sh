#!/bin/bash -l
#SBATCH --job-name=s2b_extract
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=logs/extract_sample_%j.out
#SBATCH --error=logs/extract_sample_%j.err

set -euo pipefail

# Activate conda environment
source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

DATA_ROOT=/lustre1/g/aos_shihuang/data/gtdb-r207
TAR_FILE="${DATA_ROOT}/gtdb_genomes_reps_r207.tar.gz"
GENOMES_DIR="${DATA_ROOT}/genomes_all"
METADATA_DIR="${DATA_ROOT}/metadata"
PAIR_CHUNKS_DIR="${DATA_ROOT}/pair_chunks"
N_PER_LABEL=25000
CHUNK_SIZE=1800

export GENOMES_DIR METADATA_DIR PAIR_CHUNKS_DIR N_PER_LABEL CHUNK_SIZE

echo "Checking tar file: ${TAR_FILE}"
if [ ! -f "${TAR_FILE}" ]; then
    echo "ERROR: Tar file not found"
    exit 1
fi

mkdir -p "${GENOMES_DIR}"

# Check if genomes already extracted as flat .fna files
genome_count=$(find "${GENOMES_DIR}" -maxdepth 1 -name "*.fna" | wc -l)
if [ "${genome_count}" -gt 0 ]; then
    echo "Found ${genome_count} existing .fna files in ${GENOMES_DIR}; skipping extraction."
else
    echo "Cleaning any partial extraction..."
    rm -rf "${GENOMES_DIR}"/*
    STAGING_DIR="${DATA_ROOT}/genomes_staging"
    rm -rf "${STAGING_DIR}"
    mkdir -p "${STAGING_DIR}"

    echo "Extracting tar archive..."
    tar -xzf "${TAR_FILE}" -C "${STAGING_DIR}/"

    echo "Decompressing and flattening .fna.gz files with 8 parallel gunzip jobs..."
    count=0
    max_jobs=8
    while IFS= read -r src; do
        base=$(basename "${src}" .gz)
        # Remove _genomic suffix from filename so sample script can match by accession
        dest_base=${base%_genomic.fna}.fna
        dest="${GENOMES_DIR}/${dest_base}"
        gunzip -c "${src}" > "${dest}" &
        count=$((count + 1))
        if [ "${count}" -ge "${max_jobs}" ]; then
            wait
            count=0
        fi
    done < <(find "${STAGING_DIR}" -type f -name "*.fna.gz")
    wait

    echo "Cleaning up staging directory..."
    rm -rf "${STAGING_DIR}"

    genome_count=$(find "${GENOMES_DIR}" -maxdepth 1 -name "*.fna" | wc -l)
    echo "Extracted ${genome_count} genome files"
fi

if [ "${genome_count}" -eq 0 ]; then
    echo "ERROR: No genome files found after extraction"
    exit 1
fi

# Generate manifest
python3 - <<'PY'
from pathlib import Path
import os

genomes_dir = Path(os.environ['GENOMES_DIR'])
manifest_path = Path(os.environ['METADATA_DIR']) / 'genome_manifest.txt'
files = sorted(genomes_dir.glob('*.fna'))
print(f'Found {len(files)} genome files')
with open(manifest_path, 'w') as f:
    for p in files:
        stem = p.stem
        f.write(f'{stem}\t{p}\n')
print(f'Manifest written: {manifest_path}')
PY

# Sample 100k pairs
echo "Sampling pairs..."
python3 scripts/sample_gtdb_r207_pairs_v2.py \
  --genomes "${GENOMES_DIR}" \
  --bac-metadata "${METADATA_DIR}/bac120_taxonomy_r207.tsv" \
  --ar-metadata "${METADATA_DIR}/ar53_taxonomy_r207.tsv" \
  --output results/pairs_gtdb_r207_100k.tsv \
  --n-per-label "${N_PER_LABEL}" \
  --seed 42

# Split into chunks
echo "Splitting pairs into chunks..."
python3 scripts/split_pair_chunks.py \
  --pairs results/pairs_gtdb_r207_100k.tsv \
  --output-dir "${PAIR_CHUNKS_DIR}" \
  --chunk-size "${CHUNK_SIZE}"

echo "Extraction and sampling complete at $(date)"
