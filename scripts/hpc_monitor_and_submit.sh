#!/bin/bash -l
# Monitor the GTDB-R207 tar download and auto-submit the workflow when complete.
# Run on the HPC login node with nohup.

set -euo pipefail

TAR_FILE=/lustre1/g/aos_shihuang/data/gtdb-r207/gtdb_genomes_reps_r207.tar.gz
EXPECTED_SIZE=65396676602
LOG=/lustre1/g/aos_shihuang/Syn2bANI-paper/logs/monitor_submit.log
INTERVAL_SEC=300

echo "$(date): Starting monitor loop." | tee -a "${LOG}"

while true; do
    if [ -f "${TAR_FILE}" ]; then
        ACTUAL_SIZE=$(stat -c%s "${TAR_FILE}")
        echo "$(date): Tar size ${ACTUAL_SIZE}/${EXPECTED_SIZE} ($((ACTUAL_SIZE * 100 / EXPECTED_SIZE))%)." | tee -a "${LOG}"
        if [ "${ACTUAL_SIZE}" -ge "${EXPECTED_SIZE}" ]; then
            echo "$(date): Tar complete. Running auto-continue." | tee -a "${LOG}"
            cd /lustre1/g/aos_shihuang/Syn2bANI-paper
            bash scripts/hpc_auto_continue.sh
            echo "$(date): Auto-continue finished. Exiting monitor." | tee -a "${LOG}"
            exit 0
        fi
    else
        echo "$(date): Tar file not found yet." | tee -a "${LOG}"
    fi
    sleep "${INTERVAL_SEC}"
done
