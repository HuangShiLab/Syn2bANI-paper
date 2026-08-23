#!/bin/bash
# run_skani_high_ani_final_slice.sh <slice_id> — skani dist for final high-ANI pairs.
# Output: $WORK/skani_high_ani_final/{pairid}.txt
set -uo pipefail

WORK=${SYN2BANI_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
GENOMES=${SYN2BANI_GENOMES:-$WORK/genomes_high_ani}
SKANI=${SKANI:-/lustre1/g/aos_shihuang/tools/skani-conda/bin/skani}
NSLICES=${NSLICES:-45}

SLICE=$1
PAIRS=$WORK/high_ani_pairs_final.tsv
N=$(tail -n +2 "$PAIRS" | wc -l)
CHUNK=$(( (N + NSLICES - 1) / NSLICES ))
START=$(( SLICE * CHUNK + 1 ))
END=$(( (SLICE + 1) * CHUNK ))

mkdir -p "$WORK/skani_high_ani_final"

tail -n +2 "$PAIRS" | sed -n "${START},${END}p" | while IFS=$'\t' read -r PID QA RA _; do
    OUT="$WORK/skani_high_ani_final/${PID}.txt"
    [ -s "$OUT" ] && continue
    QF="$GENOMES/${QA}.fna"
    RF="$GENOMES/${RA}.fna"
    [ -s "$QF" ] && [ -s "$RF" ] || continue
    "$SKANI" dist "$QF" "$RF" > "$OUT" 2>/dev/null || true
done

echo "[high_ani_final skani] slice $SLICE done"
