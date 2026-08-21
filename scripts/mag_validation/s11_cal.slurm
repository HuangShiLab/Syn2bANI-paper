#!/bin/bash
#SBATCH --job-name=mv_cal695
#SBATCH --partition=intel
#SBATCH --array=0-3
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/cal695_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/cal695_%A_%a.err
# s11: re-run the 695 anchor pairs with `ani --calibrate` (v5 embedded model,
#      binary @ fe0f36c) so calibrated MAG accuracy can be compared with the
#      raw gated estimates from s8.
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
OUT=$WORK/fast/per_pair_cal
mkdir -p "$WORK/logs" "$OUT"

# anchor pair list: binid <TAB> ref_path  (idempotent)
awk -F'\t' 'NR>1 && $4=="anchor" {print $1"\t"$2}' "$WORK/collect/ani_fast_tools.tsv" \
    > "$OUT/pairs_anchor.tsv"
mapfile -t ROWS < "$OUT/pairs_anchor.tsv"
N=${#ROWS[@]}
echo "[cal695] $N anchor pairs"

run_one() {
    local row=$1 bin ref
    IFS=$'\t' read -r bin ref <<< "$row"
    [ -s "$OUT/${bin}.tsv" ] && return 0
    "$S2B" ani "$WORK/bins_all/${bin}.fa" "$ref" --verbose --calibrate -t 1 \
        -o "$OUT/${bin}.tsv" 2> "$OUT/${bin}.err" \
        || echo "[cal695] FAIL $bin"
}
export -f run_one
export S2B OUT WORK

CHUNK=$(( (N + SLURM_ARRAY_TASK_COUNT - 1) / SLURM_ARRAY_TASK_COUNT ))
START=$((SLURM_ARRAY_TASK_ID * CHUNK))
END=$((START + CHUNK < N ? START + CHUNK : N))
printf '%s\n' "${ROWS[@]:$START:$((END - START))}" \
    | xargs -d '\n' -P "$SLURM_CPUS_PER_TASK" -I{} bash -c 'run_one "$@"' _ {}
echo "[cal695] task $SLURM_ARRAY_TASK_ID done: $(ls "$OUT"/*.tsv 2>/dev/null | grep -vc pairs_anchor)"
