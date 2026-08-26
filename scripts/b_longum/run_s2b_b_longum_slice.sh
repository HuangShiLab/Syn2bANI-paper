#!/bin/bash
# run_s2b_b_longum_slice.sh <slice_id> — run Syn2bANI ani for a slice of
# B. longum all-vs-all pairs. Idempotent: pairs with an existing .tsv are skipped.
set -uo pipefail

WORK=${B_LONGUM_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/b_longum_abfA}
TOOLS=${TOOLS:-/lustre1/g/aos_shihuang/tools}
S2B=${S2B:-$TOOLS/syn2bani/syn2bani}
NSLICES=${NSLICES:-20}

SLICE=$1
PAIRS=$WORK/pairs_all_vs_all.tsv
N=$(wc -l < "$PAIRS")
CHUNK=$(( (N + NSLICES - 1) / NSLICES ))
START=$(( SLICE * CHUNK + 1 ))
END=$(( (SLICE + 1) * CHUNK ))

OUT=$WORK/s2b_out
mkdir -p "$OUT" "$WORK/logs"

# Build genome path lookup
declare -A PATHS
while IFS=$'\t' read -r ACC FP; do
    PATHS[$ACC]=$FP
done < "$WORK/manifest.tsv"

tail -n +2 "$PAIRS" | sed -n "${START},${END}p" | while IFS=$'\t' read -r QA RA; do
    PID="${QA}__${RA}"
    [ -s "$OUT/${PID}.tsv" ] && continue
    QF=${PATHS[$QA]:-}
    RF=${PATHS[$RA]:-}
    if [ -z "$QF" ] || [ -z "$RF" ]; then
        echo "[b_longum-s2b] missing FASTA for $PID"
        continue
    fi
    "$S2B" ani "$QF" "$RF" --verbose --calibrate -t 1 -o "$OUT/${PID}.tsv" \
        2> "$OUT/${PID}.err" || echo "[b_longum-s2b] FAIL $PID"
done
echo "[b_longum-s2b] slice $SLICE done"
