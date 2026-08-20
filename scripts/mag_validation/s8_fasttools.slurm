#!/bin/bash
#SBATCH --job-name=mv_fasttools
#SBATCH --partition=intel
#SBATCH --array=0-11
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/fasttools_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/fasttools_%A_%a.err
# s8: syn2bani / skani / fastANI on all pairs of each sample's bins
#     (anchor + top-3 same-species alternates + nearest GTDB r207 rep).
#     Also builds reusable syn2bani sketches of the bins. PACK=3 samples/task.
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
mkdir -p "$WORK/logs" "$WORK/fast/per_pair" "$WORK/fast/rl" "$WORK/sketches"
PACK=${PACK:-3}

mapfile -t ROWS < <(tail -n +2 "$SAMPLE_LIST")

for I in $(seq $((SLURM_ARRAY_TASK_ID * PACK)) $((SLURM_ARRAY_TASK_ID * PACK + PACK - 1))); do
    [ "$I" -ge "${#ROWS[@]}" ] && break
    IFS=$'\t' read -r DS SM _ <<< "${ROWS[$I]}"
    TAG="${DS}__${SM}"
    grep -P "^${TAG}__" "$WORK/bins_all/bin_stats.tsv" | cut -f1 > "$WORK/fast/rl/${TAG}.binids.txt"
    N=$(wc -l < "$WORK/fast/rl/${TAG}.binids.txt")
    echo "[fasttools] $TAG: $N bins"
    [ "$N" -eq 0 ] && continue

    # reusable syn2bani sketches for this sample's bins (idempotent)
    SKDIR=$WORK/sketches/$TAG
    if [ ! -d "$SKDIR" ]; then
        grep -P "^${TAG}__" "$WORK/bins_all/bin_stats.tsv" | cut -f7 > "$WORK/fast/rl/${TAG}.paths.txt"
        "$S2B" sketch -o "$SKDIR" $(cat "$WORK/fast/rl/${TAG}.paths.txt") -t "$SLURM_CPUS_PER_TASK" -p \
            || echo "[fasttools] $TAG: sketch step failed (non-fatal)"
    fi

    xargs -P "$SLURM_CPUS_PER_TASK" -n 1 "$SCRIPTS/run_fasttools_bin.sh" < "$WORK/fast/rl/${TAG}.binids.txt"
    echo "[fasttools] $TAG done"
done
