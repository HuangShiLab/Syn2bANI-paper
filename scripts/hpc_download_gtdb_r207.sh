#!/bin/bash -l
#SBATCH --job-name=s2b_download_gtdb
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/download_gtdb_%j.out
#SBATCH --error=logs/download_gtdb_%j.err

set -euo pipefail

JOB_ID="${SLURM_JOB_ID:-login}"

DATA_ROOT=/lustre1/g/aos_shihuang/data/gtdb-r207
GENOMES_DIR="${DATA_ROOT}/genomes_all"
METADATA_DIR="${DATA_ROOT}/metadata"
LOG_FILE="${DATA_ROOT}/download.log"

mkdir -p "${GENOMES_DIR}" "${METADATA_DIR}" "$(dirname "$LOG_FILE")"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================"
echo "GTDB-R207 Download Script (HPC)"
echo "Data root: ${DATA_ROOT}"
echo "Start time: $(date)"
echo "========================================"

# Download metadata
echo "[1/3] Downloading GTDB-R207 metadata..."
METADATA_URL="https://data.gtdb.ecogenomic.org/releases/release207/207.0"
wget -q --show-progress -c "${METADATA_URL}/bac120_metadata_r207.tsv" \
    -O "${METADATA_DIR}/bac120_metadata_r207.tsv" || true
wget -q --show-progress -c "${METADATA_URL}/ar53_metadata_r207.tsv" \
    -O "${METADATA_DIR}/ar53_metadata_r207.tsv" || true

echo "[1/3] Metadata download complete."

# Download representative genomes via GTDB tar (avoids NCBI datasets dependency)
echo "[2/3] Downloading GTDB-R207 representative genomes..."
TAR_URL="https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/gtdb_genomes_reps_r207.tar.gz"
TAR_FILE="${DATA_ROOT}/gtdb_genomes_reps_r207.tar.gz"

if [ ! -f "${TAR_FILE}" ]; then
    wget -q --show-progress -c "${TAR_URL}" -O "${TAR_FILE}"
fi

echo "[2/3] Extracting tar archive..."
tar -xzf "${TAR_FILE}" -C "${GENOMES_DIR}/"

echo "[2/3] Genome download and extraction complete."

# Validate
echo "[3/3] Validating downloaded genomes..."
genome_count=$(find "${GENOMES_DIR}" -name "*.fna" -o -name "*.fasta" -o -name "*.fa" | wc -l)
echo "Found ${genome_count} genome files"
find "${GENOMES_DIR}" -name "*.fna" -size 0 -delete
find "${GENOMES_DIR}" -name "*.fasta" -size 0 -delete

echo "========================================"
echo "Download Summary"
echo "End time: $(date)"
echo "Genomes: ${genome_count}"
echo "========================================"
