#!/bin/bash
# run_skani_high_ani_slice.sh <slice_id> — skani dist for high-ANI pairs.
# Output: $WORK/skani_high_ani/{pairid}.txt
set -uo pipefail

WORK=${SYN2BANI_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
GENOMES=${SYN2BANI_GENOMES:-$WORK/genomes_high_ani}
SKANI=${SKANI:-/group/aos_shihuang/conda/envs/skani/bin/skani}
NSLICES=${NSLICES:-80}

SLICE=$1
PAIRS=$WORK/high_ani_candidates.tsv
N=$(tail -n +2 "$PAIRS" | wc -l)
CHUNK=$(( (N + NSLICES - 1) / NSLICES ))
START=$(( SLICE * CHUNK + 1 ))
END=$(( (SLICE + 1) * CHUNK ))

mkdir -p "$WORK/skani_high_ani"

tail -n +2 "$PAIRS" | sed -n "${START},${END}p" | while IFS=$'\t' read -r PID QA RA _; do
    OUT="$WORK/skani_high_ani/${PID}.txt"
    [ -s "$OUT" ] && continue
    QF="$GENOMES/${QA}.fna"
    RF="$GENOMES/${RA}.fna"
    [ -s "$QF" ] && [ -s "$RF" ] || continue
    "$SKANI" dist "$QF" "$RF" > "$OUT" 2>/dev/null || true
done

echo "[high_ani skani] slice $SLICE done"
