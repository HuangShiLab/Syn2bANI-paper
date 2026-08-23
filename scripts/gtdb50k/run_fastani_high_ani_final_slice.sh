#!/bin/bash
# run_fastani_high_ani_final_slice.sh <slice_id> — FastANI for final high-ANI pairs.
# Output: $WORK/fastani_high_ani_final/{pairid}.txt
set -uo pipefail

WORK=${SYN2BANI_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
GENOMES=${SYN2BANI_GENOMES:-$WORK/genomes_high_ani}
FASTANI=${FASTANI:-/group/aos_shihuang/conda/envs/fastani/bin/fastANI}
NSLICES=${NSLICES:-45}

SLICE=$1
PAIRS=$WORK/high_ani_pairs_final.tsv
N=$(tail -n +2 "$PAIRS" | wc -l)
CHUNK=$(( (N + NSLICES - 1) / NSLICES ))
START=$(( SLICE * CHUNK + 1 ))
END=$(( (SLICE + 1) * CHUNK ))

mkdir -p "$WORK/fastani_high_ani_final"

tail -n +2 "$PAIRS" | sed -n "${START},${END}p" | while IFS=$'\t' read -r PID QA RA _; do
    OUT="$WORK/fastani_high_ani_final/${PID}.txt"
    [ -s "$OUT" ] && continue
    QF="$GENOMES/${QA}.fna"
    RF="$GENOMES/${RA}.fna"
    [ -s "$QF" ] && [ -s "$RF" ] || continue
    "$FASTANI" -q "$QF" -r "$RF" -o "$OUT" > /dev/null 2>&1 || true
done

echo "[high_ani_final fastani] slice $SLICE done"
