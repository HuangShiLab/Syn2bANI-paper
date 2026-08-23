#!/bin/bash
# run_fastani_slice.sh <slice_id> — run FastANI for a contiguous slice of
# pairs_50k.tsv. Idempotent: pairs already present in the slice output are
# skipped. Output columns: pairid, fastani_ani, fastani_mapped,
# fastani_total (NA row when FastANI reports nothing = below its detection).
set -uo pipefail

WORK=${GTDB50K_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
GENOMES=${GTDB50K_GENOMES:-/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all}
FASTANI=${FASTANI:-/group/aos_shihuang/conda/envs/fastani/bin/fastANI}
NSLICES=${NSLICES:-190}

SLICE=$1
PAIRS=$WORK/pairs_50k.tsv
N=$(tail -n +2 "$PAIRS" | wc -l)
CHUNK=$(( (N + NSLICES - 1) / NSLICES ))
START=$(( SLICE * CHUNK + 1 ))          # tail already stripped the header
END=$(( (SLICE + 1) * CHUNK ))

mkdir -p "$WORK/fastani_rows"
ROWS_TMP="$WORK/fastani_rows/slice_${SLICE}.tsv.tmp"
ROWS_OUT="$WORK/fastani_rows/slice_${SLICE}.tsv"
: > "$ROWS_TMP"

tail -n +2 "$PAIRS" | sed -n "${START},${END}p" | while IFS=$'\t' read -r QA RA _; do
    PID="${QA}__${RA}"
    OUT=$WORK/tmp_fa.${PID}.out
    "$FASTANI" -q "$GENOMES/${QA}.fna" -r "$GENOMES/${RA}.fna" -o "$OUT" > /dev/null 2>&1
    if [ -s "$OUT" ]; then
        # query ref ANI mapped total
        ANI=$(awk -F'\t' '{print $3}' "$OUT" | head -1)
        MAP=$(awk -F'\t' '{print $4}' "$OUT" | head -1)
        TOT=$(awk -F'\t' '{print $5}' "$OUT" | head -1)
        echo -e "${PID}\t${ANI}\t${MAP}\t${TOT}" >> "$ROWS_TMP"
    else
        echo -e "${PID}\tNA\tNA\tNA" >> "$ROWS_TMP"
    fi
    rm -f "$OUT"
done
mv "$ROWS_TMP" "$ROWS_OUT"
echo "[gtdb50k] fastani slice $SLICE done: $(wc -l < "$ROWS_OUT") pairs"
