#!/bin/bash
# run_minimap2_slice.sh <slice_id> — run minimap2 for a contiguous slice of
# pairs_50k.tsv, parse SV metrics on the fly, and delete PAF immediately to
# save disk space. Idempotent: pairs already present in the output TSV are
# skipped.
set -uo pipefail

WORK=${GTDB50K_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
SCRIPTS=${GTDB50K_SCRIPTS:-$WORK/scripts}
GENOMES=${GTDB50K_GENOMES:-/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all}
# software lives in /lustre1/g/aos_shihuang/tools to save shared HPC space
TOOLS=${TOOLS:-/lustre1/g/aos_shihuang/tools}
MM2=${MM2:-$TOOLS/minimap2/minimap2}
PY=${PY:-$TOOLS/python3/bin/python3}
[ -x "$PY" ] || PY=$TOOLS/anaconda3/bin/python3
NSLICES=${NSLICES:-190}

SLICE=$1
PAIRS=$WORK/pairs_50k.tsv
N=$(tail -n +2 "$PAIRS" | wc -l)
CHUNK=$(( (N + NSLICES - 1) / NSLICES ))
START=$(( SLICE * CHUNK + 1 ))
END=$(( (SLICE + 1) * CHUNK ))

mkdir -p "$WORK/minimap2_rows" "$WORK/logs"
ROWS_TMP="$WORK/minimap2_rows/slice_${SLICE}.tsv.tmp"
ROWS_OUT="$WORK/minimap2_rows/slice_${SLICE}.tsv"
# Append to tmp so a killed job can resume; final sort/de-dupe happens at the end.
[ -f "$ROWS_TMP" ] || : > "$ROWS_TMP"

tail -n +2 "$PAIRS" | sed -n "${START},${END}p" | while IFS=$'\t' read -r QA RA _; do
    PID="${QA}__${RA}"
    # Skip if already finalized or already in the current tmp (resume-friendly).
    if [ -s "$ROWS_OUT" ] && grep -qm1 "^${PID}" "$ROWS_OUT"; then
        continue
    fi
    if grep -qm1 "^${PID}" "$ROWS_TMP" 2>/dev/null; then
        continue
    fi
    QF=$GENOMES/${QA}.fna
    RF=$GENOMES/${RA}.fna
    if [ ! -s "$QF" ] || [ ! -s "$RF" ]; then
        echo -e "${PID}\tNA\tNA\tNA\tmissing_fasta" >> "$ROWS_TMP"
        continue
    fi
    TMP=$(mktemp -d "$WORK/tmp.mm2.${PID}.XXXXXX")
    if "$MM2" -x asm5 -c --eqx "$RF" "$QF" > "$TMP/mm2.paf" 2> /dev/null; then
        $PY "$SCRIPTS/parse_minimap2_sv.py" "$PID" "$TMP/mm2.paf" >> "$ROWS_TMP"
    else
        echo -e "${PID}\tNA\tNA\tNA\tminimap2_failed" >> "$ROWS_TMP"
    fi
    rm -rf "$TMP"
done

# merge with existing output if present (idempotent) and de-duplicate
cat "$ROWS_TMP" > "$WORK/minimap2_rows/slice_${SLICE}.tsv.all"
if [ -s "$ROWS_OUT" ]; then
    cat "$ROWS_OUT" >> "$WORK/minimap2_rows/slice_${SLICE}.tsv.all"
fi
sort -u "$WORK/minimap2_rows/slice_${SLICE}.tsv.all" > "$ROWS_OUT"
rm -f "$ROWS_TMP" "$WORK/minimap2_rows/slice_${SLICE}.tsv.all"
echo "[minimap2] slice $SLICE done: $(wc -l < "$ROWS_OUT") pairs"
