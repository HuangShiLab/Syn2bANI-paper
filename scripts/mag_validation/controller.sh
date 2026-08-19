#!/bin/bash
#SBATCH --job-name=mv_ctl
#SBATCH --partition=intel
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:10:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/ctl_%A.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/ctl_%A.err
# controller.sh <stage> — incremental submitter. The site caps MaxSubmit/MaxJobs
# per user (~50/45, array elements count), so the chain is submitted
# stage-by-stage: each controller fires after its dependency completes, submits
# the next stage(s), records job IDs in $WORK/jobs/jobs.tsv, and submits the
# next controller. Retries with backoff + delayed self-resubmit on submit-limit
# failures, so it survives other jobs occupying the user's quota.
set -uo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
cd "$SCRIPTS_HOME"
source ./common.sh
JOBS=$WORK/jobs/jobs.tsv
mkdir -p "$WORK/jobs" "$WORK/logs"
touch "$JOBS"
CUR_STAGE=$1
HERE=$SCRIPTS_HOME

jid() { grep -P "^$1\t" "$JOBS" | tail -1 | cut -f2; }
have() { grep -qP "^$1\t\S" "$JOBS"; }
record() { echo -e "$1\t$2" >> "$JOBS"; }

retry_later() { # schedule this same controller stage again in 20 min
    local rid
    rid=$(sbatch --parsable --export=NONE --begin=now+20min "$HERE/controller.sh" "$CUR_STAGE" 2>/dev/null) || rid=""
    if [ -n "$rid" ]; then
        record "ctl_retry_${CUR_STAGE}" "$rid"
        echo "[ctl] $CUR_STAGE: submit limit hit; retry scheduled in 20 min (job $rid)"
    else
        echo "[ctl] $CUR_STAGE: FATAL — could not even schedule a retry; resubmit manually: sbatch $HERE/controller.sh $CUR_STAGE" >&2
        exit 1
    fi
}

submit() { # name dep script [extra sbatch args...]
    local name=$1 dep=$2 script=$3; shift 3
    if have "$name"; then echo "[ctl] $name already submitted: $(jid "$name")"; return 0; fi
    local id="" attempt
    for attempt in 1 2 3; do
        if [ -n "$dep" ]; then
            id=$(sbatch --parsable --export=NONE --dependency="afterok:$dep" "$@" "$script" 2>/dev/null) || id=""
        else
            id=$(sbatch --parsable --export=NONE "$@" "$script" 2>/dev/null) || id=""
        fi
        [ -n "$id" ] && break
        echo "[ctl] $name submit attempt $attempt failed (quota?); sleeping 120s"
        sleep 120
    done
    [ -z "$id" ] && { retry_later; exit 0; }
    record "$name" "$id"
    echo "[ctl] submitted $name = $id (dep: ${dep:-none})"
}

ctl() { # dep next_stage
    local id
    id=$(sbatch --parsable --export=NONE --dependency="afterok:$1" "$HERE/controller.sh" "$2" 2>/dev/null) || id=""
    [ -z "$id" ] && { retry_later; exit 0; }
    record "ctl_$2" "$id"
    echo "[ctl] controller $2 = $id (dep: $1)"
}

case "$CUR_STAGE" in
    after_s1)
        submit s2_depth "$(jid s1_assemble)" s2_depth.slurm
        ctl "$(jid s2_depth)" after_s2
        ;;
    after_s2)
        submit s3_bin "$(jid s2_depth)" s3_bin.slurm
        ctl "$(jid s3_bin)" after_s3
        ;;
    after_s3)
        submit s4_checkm2 "$(jid s3_bin)" s4_checkm2.slurm
        submit s5_assign "$(jid s3_bin):$(jid s0_sourceprep)" s5_assign.slurm
        submit s7_repsearch "$(jid s4_checkm2):$(jid s0b_gtdb_prep)" s7_repsearch.slurm
        ctl "$(jid s4_checkm2):$(jid s5_assign)" after_s45
        ;;
    after_s45)
        submit s6_cohort "$(jid s4_checkm2):$(jid s5_assign)" s6_cohort.slurm
        ctl "$(jid s6_cohort):$(jid s7_repsearch)" after_s67
        ;;
    after_s67)
        submit s8_fasttools "$(jid s6_cohort):$(jid s7_repsearch)" s8_fasttools.slurm
        NPAIRS=$(( $(wc -l < "$WORK/pairs/pairs_anchor.tsv") - 1 ))
        NCH=$(( (NPAIRS + DNADIFF_CHUNK - 1) / DNADIFF_CHUNK ))
        [ "$NCH" -lt 1 ] && NCH=1
        echo "[ctl] truth: $NPAIRS pairs -> $NCH chunks"
        submit s9_truth "$(jid s6_cohort)" s9_truth.slurm --array=0-$((NCH-1))
        ctl "$(jid s8_fasttools):$(jid s9_truth)" after_s89
        ;;
    after_s89)
        submit s10_collect "$(jid s8_fasttools):$(jid s9_truth)" s10_collect.slurm
        echo "[ctl] full chain submitted; collect = $(jid s10_collect)"
        ;;
    *)
        echo "[ctl] unknown stage: $CUR_STAGE" >&2
        exit 1
        ;;
esac
