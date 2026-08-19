#!/bin/bash
#SBATCH --job-name=mv_megahit
#SBATCH --partition=intel
#SBATCH --array=0-17
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/megahit_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/megahit_%A_%a.err
# s1: MEGAHIT assembly (interleaved paired reads). PACK=2 samples per task
# (site MaxJobs limit makes 35-task arrays infeasible alongside other jobs).
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
mkdir -p "$WORK/logs"
PACK=${PACK:-2}

mapfile -t ROWS < <(tail -n +2 "$SAMPLE_LIST")
export PATH="$(dirname "$MEGAHIT"):$PATH"

for I in $(seq $((SLURM_ARRAY_TASK_ID * PACK)) $((SLURM_ARRAY_TASK_ID * PACK + PACK - 1))); do
    [ "$I" -ge "${#ROWS[@]}" ] && break
    IFS=$'\t' read -r DS SM _ <<< "${ROWS[$I]}"
    READS=$(reads_path "$DS" "$SM")
    TAG="${DS}__${SM}"
    OUT=$WORK/asm/$TAG
    echo "[assemble] $TAG reads=$READS"
    if [ -s "$OUT/final.contigs.fa" ]; then
        echo "[assemble] $TAG already assembled, skipping"
        continue
    fi
    rm -rf "$OUT"
    mkdir -p "$OUT"
    $MEGAHIT --12 "$READS" -t "$SLURM_CPUS_PER_TASK" -o "$OUT" --min-contig-len 1000
    echo "[assemble] $TAG done: $(grep -c '^>' "$OUT/final.contigs.fa") contigs"
done
