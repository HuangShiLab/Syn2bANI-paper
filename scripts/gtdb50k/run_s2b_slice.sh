#!/bin/bash
# run_s2b_slice.sh <slice_id> — syn2bani ani --verbose --calibrate for a
# contiguous slice of pairs_50k.tsv. Idempotent: pairs with an existing
# output row are skipped. Same slicing arithmetic as run_dnadiff_slice.sh.
set -uo pipefail

WORK=${GTDB50K_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
GENOMES=${GTDB50K_GENOMES:-/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all}
S2B=${S2B:-/lustre1/g/aos_shihuang/Syn2bANI-hi95/target/release/syn2bani}
NSLICES=${NSLICES:-190}

SLICE=$1
PAIRS=$WORK/pairs_50k.tsv
N=$(tail -n +2 "$PAIRS" | wc -l)
CHUNK=$(( (N + NSLICES - 1) / NSLICES ))
START=$(( SLICE * CHUNK + 1 ))
END=$(( (SLICE + 1) * CHUNK ))

OUT=$WORK/s2b_out
mkdir -p "$OUT"

tail -n +2 "$PAIRS" | sed -n "${START},${END}p" | while IFS=$'\t' read -r QA RA _; do
    PID="${QA}__${RA}"
    [ -s "$OUT/${PID}.tsv" ] && continue
    "$S2B" ani "$GENOMES/${QA}.fna" "$GENOMES/${RA}.fna" \
        --verbose --calibrate -t 1 -o "$OUT/${PID}.tsv" \
        2> "$OUT/${PID}.err" || echo "[gtdb50k-s2b] FAIL $PID"
done
echo "[gtdb50k-s2b] slice $SLICE done"
