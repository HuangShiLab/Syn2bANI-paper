#!/bin/bash -l
# Auto-continue the Syn2bANI HPC workflow after the GTDB-R207 tar download finishes.
# Run this on the HPC login node once wget has completed, or schedule it via cron.

set -euo pipefail

TAR_FILE=/lustre1/g/aos_shihuang/data/gtdb-r207/gtdb_genomes_reps_r207.tar.gz
EXPECTED_SIZE=65396676602
LOG=/lustre1/g/aos_shihuang/Syn2bANI-paper/logs/auto_continue.log
PAIR_CHUNKS_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/pair_chunks

echo "$(date): Checking tar file..." | tee -a "${LOG}"

if [ ! -f "${TAR_FILE}" ]; then
    echo "$(date): Tar file not found yet. Exiting." | tee -a "${LOG}"
    exit 0
fi

ACTUAL_SIZE=$(stat -c%s "${TAR_FILE}")
if [ "${ACTUAL_SIZE}" -lt "${EXPECTED_SIZE}" ]; then
    echo "$(date): Tar file incomplete (${ACTUAL_SIZE}/${EXPECTED_SIZE}). Exiting." | tee -a "${LOG}"
    exit 0
fi

echo "$(date): Tar file complete. Submitting extract+sample job." | tee -a "${LOG}"

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

EXTRACT_JOB=$(sbatch scripts/hpc_extract_and_sample.sh | awk '{print $4}')
echo "$(date): Extract job: ${EXTRACT_JOB}" | tee -a "${LOG}"

# Determine array size from existing chunks, or use safe default if not yet created
if [ -d "${PAIR_CHUNKS_DIR}" ]; then
    N_CHUNKS=$(ls "${PAIR_CHUNKS_DIR}"/pairs_chunk_*.tsv 2>/dev/null | wc -l)
else
    N_CHUNKS=0
fi
if [ "${N_CHUNKS}" -eq 0 ]; then
    ARRAY_SPEC="1-100"
    echo "$(date): Chunks not yet created, using default array ${ARRAY_SPEC}." | tee -a "${LOG}"
else
    ARRAY_SPEC="1-${N_CHUNKS}"
    echo "$(date): Detected ${N_CHUNKS} pair chunks, array ${ARRAY_SPEC}." | tee -a "${LOG}"
fi

# Submit the full workflow with dependency on extraction
MATRIX=$(sbatch --dependency=afterok:"${EXTRACT_JOB}" --array="${ARRAY_SPEC}" scripts/slurm_matrix_array.sh | awk '{print $4}')
echo "$(date): Matrix array job: ${MATRIX}" | tee -a "${LOG}"

MERGE=$(sbatch --dependency=afterok:"${MATRIX}" scripts/slurm_merge_phase1.sh | awk '{print $4}')
echo "$(date): Merge+sample job: ${MERGE}" | tee -a "${LOG}"

FASTANI=$(sbatch --dependency=afterok:"${MERGE}" scripts/slurm_fastani_array.sh | awk '{print $4}')
echo "$(date): FastANI array job: ${FASTANI}" | tee -a "${LOG}"

TRAIN=$(sbatch --dependency=afterok:"${FASTANI}" scripts/slurm_train_final.sh | awk '{print $4}')
echo "$(date): Final train job: ${TRAIN}" | tee -a "${LOG}"

echo "$(date): Full workflow submitted." | tee -a "${LOG}"
