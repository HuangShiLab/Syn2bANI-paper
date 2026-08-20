#!/bin/bash
#SBATCH --job-name=mv_depth
#SBATCH --partition=intel
#SBATCH --array=0-17
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=18:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/depth_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/depth_%A_%a.err
# s2: read depth per contig for MetaBAT2 (bowtie2 map interleaved reads -> jgi
#     depth). PACK=2 samples per task.
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
mkdir -p "$WORK/logs"
PACK=${PACK:-2}

mapfile -t ROWS < <(tail -n +2 "$SAMPLE_LIST")
export PATH="$BT2_BIN:$PATH"

for I in $(seq $((SLURM_ARRAY_TASK_ID * PACK)) $((SLURM_ARRAY_TASK_ID * PACK + PACK - 1))); do
    [ "$I" -ge "${#ROWS[@]}" ] && break
    IFS=$'\t' read -r DS SM _ <<< "${ROWS[$I]}"
    READS=$(reads_path "$DS" "$SM")
    TAG="${DS}__${SM}"
    CTG=$WORK/asm/$TAG/final.contigs.fa
    OUT=$WORK/depth/$TAG
    mkdir -p "$OUT"
    if [ -s "$OUT/depth.txt" ]; then
        echo "[depth] $TAG already done, skipping"
        continue
    fi
    if [ ! -s "$OUT/idx.1.bt2" ]; then
        bowtie2-build --threads "$SLURM_CPUS_PER_TASK" "$CTG" "$OUT/idx" > /dev/null
    fi
    bowtie2 -x "$OUT/idx" --interleaved "$READS" -p "$SLURM_CPUS_PER_TASK" --no-unal 2> "$OUT/bowtie2.log" \
        | samtools sort -@ 8 -m 2G -o "$OUT/sorted.bam" -
    "$MB2_BIN/jgi_summarize_bam_contig_depths" --outputDepth "$OUT/depth.txt" "$OUT/sorted.bam"
    rm -f "$OUT/sorted.bam"   # depth.txt is all MetaBAT2 needs; keep lustre clean
    echo "[depth] $TAG done"
done
