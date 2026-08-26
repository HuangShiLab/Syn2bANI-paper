#!/bin/bash
# s4_struct_vs_ref.sh — struct every B. longum genome vs the abfA+ reference
# FSHHK16M1_ctg (abfA cluster: contig10 8,546-37,075, 28.5 kb).
# Output: $WORK/struct_vs_ref/<acc>.struct.tsv
set -uo pipefail

WORK=${B_LONGUM_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/b_longum_abfA}
GENOMES=${B_LONGUM_GENOMES:-/lustre1/g/aos_shihuang/Strain2b/data/JNU_genomes/genome2023/fna}
REF=${ABFA_REF:-$GENOMES/FSHHK16M1_ctg.fna}
S2B=${S2B:-/lustre1/g/aos_shihuang/tools/syn2bani/syn2bani}
OUT="$WORK/struct_vs_ref"
mkdir -p "$OUT"

tail -n +2 "$WORK/manifest.tsv" | while IFS=$'\t' read -r ACC FP; do
    [ "$FP" = "$REF" ] && continue
    O="$OUT/${ACC}.struct.tsv"
    [ -s "$O" ] && continue
    "$S2B" struct "$FP" "$REF" -o "$O" 2>/dev/null || echo "FAIL $ACC"
done
echo "struct_vs_ref done: $(ls "$OUT" | wc -l) genomes"
