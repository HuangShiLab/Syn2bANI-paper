#!/bin/bash
# run_dnadiff_high_ani_slice.sh <slice_id> — dnadiff ANIm truth for high-ANI
# candidate pairs. Idempotent: pairs with existing dd.report are skipped.
set -uo pipefail

WORK=${SYN2BANI_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
SCRIPTS=${SYN2BANI_SCRIPTS:-$WORK/scripts}
GENOMES=${SYN2BANI_GENOMES:-$WORK/genomes_high_ani}
ANVIO_BIN=${ANVIO_BIN:-/group/aos_shihuang/conda/envs/anvio/bin}
PY=${PY:-/group/aos_shihuang/conda/bin/python3}
NSLICES=${NSLICES:-80}

SLICE=$1
PAIRS=$WORK/high_ani_pairs_ready.tsv
N=$(tail -n +2 "$PAIRS" | wc -l)
CHUNK=$(( (N + NSLICES - 1) / NSLICES ))
START=$(( SLICE * CHUNK + 1 ))
END=$(( (SLICE + 1) * CHUNK ))

mkdir -p "$WORK/out_high_ani" "$WORK/rows_high_ani"
ROWS_TMP="$WORK/rows_high_ani/slice_${SLICE}.tsv.tmp"
ROWS_OUT="$WORK/rows_high_ani/slice_${SLICE}.tsv"
export PATH="$ANVIO_BIN:$PATH"
: > "$ROWS_TMP"

tail -n +2 "$PAIRS" | sed -n "${START},${END}p" | while IFS=$'\t' read -r PID QA RA _; do
    DD=$WORK/out_high_ani/$PID
    if [ ! -s "$DD/dd.report" ]; then
        TMP=$(mktemp -d "$WORK/tmp_ha.${PID}.XXXXXX")
        for ACC in "$QA" "$RA"; do
            SRC=$GENOMES/${ACC}.fna
            if [ ! -s "$SRC" ]; then
                echo -e "${PID}\tNA\tNA\tNA\tNA\tmissing_genome" >> "$ROWS_TMP"
                rm -rf "$TMP"
                continue 2
            fi
            if grep -qm1 -P "\s" "$SRC"; then
                awk '/^>/{sub(/[ \t].*$/,"");print;next}{gsub(/[ \t]/,"");print}' \
                    "$SRC" > "$TMP/${ACC}.fna"
            else
                ln -s "$SRC" "$TMP/${ACC}.fna"
            fi
        done
        mkdir -p "$DD"
        (cd "$DD" && dnadiff -p dd "$TMP/${RA}.fna" "$TMP/${QA}.fna" > /dev/null 2>&1) \
            || echo "[high_ani] $PID: dnadiff failed"
        rm -rf "$TMP"
        rm -f "$DD"/dd.delta "$DD"/dd.mdelta "$DD"/dd.1delta "$DD"/dd.mcoords \
              "$DD"/dd.qdiff "$DD"/dd.rdiff "$DD"/dd.snps 2>/dev/null
    fi
    if [ -s "$DD/dd.report" ]; then
        $PY "$SCRIPTS/parse_dnadiff_pair.py" "$PID" "$DD/dd.report" >> "$ROWS_TMP"
    else
        echo -e "${PID}\tNA\tNA\tNA\tNA\tdnadiff_failed" >> "$ROWS_TMP"
    fi
done
mv "$ROWS_TMP" "$ROWS_OUT"
echo "[high_ani dnadiff] slice $SLICE done: $(wc -l < "$ROWS_OUT") pairs"
