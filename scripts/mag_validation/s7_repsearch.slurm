#!/bin/bash
#SBATCH --job-name=mv_repsearch
#SBATCH --partition=intel
#SBATCH --array=0-11
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/repsearch_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/repsearch_%A_%a.err
# s7: skani search of each sample's staged bins against the GTDB r207 sketch DB;
#     nearest rep per bin -> repsearch/{tag}.map.tsv. PACK=3 samples per task.
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
mkdir -p "$WORK/logs" "$WORK/repsearch"
PACK=${PACK:-3}

mapfile -t ROWS < <(tail -n +2 "$SAMPLE_LIST")

for I in $(seq $((SLURM_ARRAY_TASK_ID * PACK)) $((SLURM_ARRAY_TASK_ID * PACK + PACK - 1))); do
    [ "$I" -ge "${#ROWS[@]}" ] && break
    IFS=$'\t' read -r DS SM _ <<< "${ROWS[$I]}"
    TAG="${DS}__${SM}"
    OUT=$WORK/repsearch/$TAG
    mkdir -p "$OUT"
    if [ -s "$OUT/map.tsv" ]; then
        echo "[repsearch] $TAG already done, skipping"
        continue
    fi
    grep -P "^${DS}__${SM}__" "$WORK/bins_all/bin_stats.tsv" | cut -f7 > "$OUT/ql.txt"
    N=$(wc -l < "$OUT/ql.txt")
    echo "[repsearch] $TAG: $N bins"
    if [ "$N" -eq 0 ]; then
        touch "$OUT/map.tsv"
        continue
    fi
    $SKANI search -d "$WORK/refs/gtdb/db" --ql "$OUT/ql.txt" -n 5 --min-af 0 \
        -o "$OUT/search.tsv" -t "$SLURM_CPUS_PER_TASK"
    $PY "$SCRIPTS/parse_repsearch.py" "$OUT/search.tsv" > "$OUT/map.tsv"
    echo "[repsearch] $TAG done"
done
