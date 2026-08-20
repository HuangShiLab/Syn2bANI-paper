#!/bin/bash
# submit_all.sh — bootstrap the mag_validation chain. The site caps
# MaxSubmit/MaxJobs per user (~50/45, array elements count), so only s0/s0b/s1
# + the first controller are submitted here; controller.sh jobs submit the rest
# as dependencies complete (with retry on quota). Run on the HPC login node
# from $WORK/scripts. Safe to re-run: stages already in jobs/jobs.tsv are
# skipped (delete a line there to force resubmission).
set -uo pipefail
cd "$(dirname "$0")"
source ./common.sh
mkdir -p "$WORK/logs" "$WORK/jobs"
JOBS=$WORK/jobs/jobs.tsv
touch "$JOBS"
HERE=$(pwd)

jid() { grep -P "^$1\t" "$JOBS" | tail -1 | cut -f2; }
have() { grep -qP "^$1\t\S" "$JOBS"; }
record() { echo -e "$1\t$2" >> "$JOBS"; }

submit() { # name dep script
    local name=$1 dep=$2 script=$3
    if have "$name"; then echo "[submit] $name already submitted: $(jid "$name")"; return 0; fi
    local id="" attempt
    for attempt in 1 2 3; do
        if [ -n "$dep" ]; then
            id=$(sbatch --parsable --export=NONE --dependency="afterok:$dep" "$script" 2>/dev/null) || id=""
        else
            id=$(sbatch --parsable --export=NONE "$script" 2>/dev/null) || id=""
        fi
        [ -n "$id" ] && break
        echo "[submit] $name attempt $attempt failed (quota?); sleeping 120s"
        sleep 120
    done
    if [ -z "$id" ]; then
        echo "[submit] FATAL: could not submit $name after 3 attempts. Re-run submit_all.sh later." >&2
        exit 1
    fi
    record "$name" "$id"
    echo "[submit] $name = $id (dep: ${dep:-none})"
}

submit s0_sourceprep "" s0_sourceprep.slurm
submit s0b_gtdb_prep "" s0b_gtdb_prep.slurm
submit s1_assemble "" s1_assemble.slurm
if ! have ctl_after_s1; then
    CID=$(sbatch --parsable --export=NONE --dependency="afterok:$(jid s1_assemble)" "$HERE/controller.sh" after_s1 2>/dev/null) || CID=""
    if [ -n "$CID" ]; then
        record ctl_after_s1 "$CID"
        echo "[submit] ctl_after_s1 = $CID"
    else
        echo "[submit] FATAL: could not submit controller; re-run submit_all.sh later." >&2
        exit 1
    fi
fi

echo "=== jobs.tsv ==="
cat "$JOBS"
