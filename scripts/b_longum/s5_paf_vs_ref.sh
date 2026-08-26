#!/bin/bash
# s5_paf_vs_ref.sh <REF_ACC> — struct --paf of every B. longum genome vs the
# given reference (default FSHHK16M1_ctg). Output: $WORK/paf_vs_<REF>/.
set -uo pipefail
WORK=${B_LONGUM_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/b_longum_abfA}
GENOMES=${B_LONGUM_GENOMES:-/lustre1/g/aos_shihuang/Strain2b/data/JNU_genomes/genome2023/fna}
S2B=${S2B:-/lustre1/g/aos_shihuang/tools/syn2bani/syn2bani}
REF=${1:-FSHHK16M1_ctg}
OUT="$WORK/paf_vs_${REF}"
mkdir -p "$OUT"
tail -n +2 "$WORK/manifest.tsv" | while IFS=$'\t' read -r ACC FP; do
    [ "$ACC" = "$REF" ] && continue
    O="$OUT/${ACC}.paf"
    [ -s "$O" ] && continue
    "$S2B" struct "$FP" "$GENOMES/${REF}.fna" --paf -o "$O" 2>/dev/null
done
echo "$OUT: $(ls "$OUT" | wc -l) files"
