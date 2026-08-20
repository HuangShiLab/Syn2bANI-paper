#!/bin/bash
#SBATCH --job-name=mv_metabat2
#SBATCH --partition=intel
#SBATCH --array=0-11
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/bin_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/bin_%A_%a.err
# s3: MetaBAT2 binning per sample. PACK=3 samples per task.
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
mkdir -p "$WORK/logs"
PACK=${PACK:-3}

mapfile -t ROWS < <(tail -n +2 "$SAMPLE_LIST")
export PATH="$MB2_BIN:$PATH"

for I in $(seq $((SLURM_ARRAY_TASK_ID * PACK)) $((SLURM_ARRAY_TASK_ID * PACK + PACK - 1))); do
    [ "$I" -ge "${#ROWS[@]}" ] && break
    IFS=$'\t' read -r DS SM _ <<< "${ROWS[$I]}"
    TAG="${DS}__${SM}"
    CTG=$WORK/asm/$TAG/final.contigs.fa
    OUT=$WORK/bins/$TAG
    if ls "$OUT"/bin.*.fa > /dev/null 2>&1; then
        echo "[bin] $TAG already binned, skipping"
        continue
    fi
    rm -rf "$OUT"
    mkdir -p "$OUT"
    metabat2 -i "$CTG" -a "$WORK/depth/$TAG/depth.txt" -o "$OUT/bin" \
        -m "$MB2_MIN_CONTIG" -t "$SLURM_CPUS_PER_TASK"
    echo "[bin] $TAG done: $(ls "$OUT"/bin.*.fa 2>/dev/null | wc -l) bins"
done
