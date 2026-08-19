#!/bin/bash
#SBATCH --job-name=mv_assign
#SBATCH --partition=intel
#SBATCH --array=0-11
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/assign_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/assign_%A_%a.err
# s5: contig truth assignment per sample: bin contigs -> dataset source genomes
#     (minimap2 asm5 + CIGAR; >=80% contig coverage, >=95% identity). PACK=3.
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
mkdir -p "$WORK/logs" "$WORK/assign"
PACK=${PACK:-3}

mapfile -t ROWS < <(tail -n +2 "$SAMPLE_LIST")

for I in $(seq $((SLURM_ARRAY_TASK_ID * PACK)) $((SLURM_ARRAY_TASK_ID * PACK + PACK - 1))); do
    [ "$I" -ge "${#ROWS[@]}" ] && break
    IFS=$'\t' read -r DS SM _ <<< "${ROWS[$I]}"
    TAG="${DS}__${SM}"
    REF=$WORK/refs/$DS/source_all.fa
    OUT=$WORK/assign/$TAG
    mkdir -p "$OUT"
    if [ -s "$OUT/assign.tsv" ]; then
        echo "[assign] $TAG already done, skipping"
        continue
    fi
    $PY "$SCRIPTS/prep_sample_bins.py" "$WORK/bins/$TAG" "$DS" "$SM" "$OUT/sample_bins.fa" "$MIN_BIN_BP"
    $MINIMAP2 -x asm5 -c --eqx -t "$SLURM_CPUS_PER_TASK" "$REF" "$OUT/sample_bins.fa" > "$OUT/map.paf"
    ASSIGN_MIN_COV=$ASSIGN_MIN_COV ASSIGN_MIN_ID=$ASSIGN_MIN_ID \
        $PY "$SCRIPTS/parse_assign.py" "$OUT/sample_bins.fa" "$OUT/map.paf" \
        "$WORK/refs/$DS/genome_lengths.tsv" "$WORK/refs/$DS/species_clusters.tsv" "$OUT/assign.tsv"
    echo "[assign] $TAG done"
done
