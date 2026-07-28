#!/bin/bash
#SBATCH --job-name=anim
#SBATCH --array=1-32
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=anim_%A_%a.out
# dnadiff over the sampled pairs. ANIm is the independent reference.
set -euo pipefail
source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani
SAMPLE="/lustre1/g/aos_shihuang/Syn2bANI-paper/results/sample_anim_truth.tsv"
CHUNK=65
OUT="${OUT:-anim_results}"
mkdir -p "$OUT"
start=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 2 ))   # +2 skips the header
end=$(( start + CHUNK - 1 ))
sed -n "${start},${end}p" "$SAMPLE" | while IFS=$'\t' read -r q r band grp v qp rp; do
    [ -z "${q:-}" ] && continue
    tag="${q}__${r}"
    pre="$OUT/$tag"
    [ -s "$pre.report" ] && continue
    dnadiff -p "$pre" "$rp" "$qp" >/dev/null 2>&1 || { echo "FAILED $tag" >&2; continue; }
    # AvgIdentity from the 1-to-1 block is the ANIm value.
    ident=$(awk '/^AvgIdentity/ {print $2; exit}' "$pre.report")
    aln=$(awk '/^AlignedBases/ {print $2; exit}' "$pre.report")
    printf '%s\t%s\t%s\t%s\n' "$q" "$r" "${ident:-NA}" "${aln:-NA}" \
        >> "$OUT/anim_${SLURM_ARRAY_TASK_ID}.tsv"
    rm -f "$pre".{delta,1delta,mdelta,1coords,mcoords,qdiff,rdiff,snps,unref,unqry}
done
