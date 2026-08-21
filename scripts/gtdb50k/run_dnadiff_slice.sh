#!/bin/bash
# run_dnadiff_slice.sh <slice_id> — run dnadiff for a contiguous slice of
# pairs_50k.tsv. Idempotent: pairs with an existing dd.report are skipped.
# Per-pair scratch under $WORK/tmp (sanitized FASTA copies), deleted after.
# Keeps dd.report and dd.1coords; deletes the large delta files.
set -uo pipefail

WORK=${GTDB50K_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k}
SCRIPTS=${GTDB50K_SCRIPTS:-$WORK/scripts}
GENOMES=${GTDB50K_GENOMES:-/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all}
ANVIO_BIN=${ANVIO_BIN:-/group/aos_shihuang/conda/envs/anvio/bin}
PY=${PY:-/group/aos_shihuang/conda/bin/python3}
NSLICES=${NSLICES:-49}

SLICE=$1
PAIRS=$WORK/pairs_50k.tsv
N=$(tail -n +2 "$PAIRS" | wc -l)
CHUNK=$(( (N + NSLICES - 1) / NSLICES ))
START=$(( SLICE * CHUNK + 2 ))          # +2: skip header, 1-based sed
END=$(( (SLICE + 1) * CHUNK + 1 ))

mkdir -p "$WORK/out" "$WORK/rows"
ROWS_TMP="$WORK/rows/slice_${SLICE}.tsv.tmp"
ROWS_OUT="$WORK/rows/slice_${SLICE}.tsv"
export PATH="$ANVIO_BIN:$PATH"
: > "$ROWS_TMP"

tail -n +2 "$PAIRS" | sed -n "${START},${END}p" | while IFS=$'\t' read -r QA RA _; do
    PID="${QA}__${RA}"
    DD=$WORK/out/$PID
    if [ ! -s "$DD/dd.report" ]; then
        TMP=$(mktemp -d "$WORK/tmp.${PID}.XXXXXX")
        for ACC in "$QA" "$RA"; do
            SRC=$GENOMES/${ACC}.fna
            if grep -qm1 -P "\s" "$SRC"; then
                awk '/^>/{sub(/[ \t].*$/,"");print;next}{gsub(/[ \t]/,"");print}' \
                    "$SRC" > "$TMP/${ACC}.fna"
            else
                ln -s "$SRC" "$TMP/${ACC}.fna"
            fi
        done
        mkdir -p "$DD"
        (cd "$DD" && dnadiff -p dd "$TMP/${RA}.fna" "$TMP/${QA}.fna" > /dev/null 2>&1) \
            || echo "[gtdb50k] $PID: dnadiff failed"
        rm -rf "$TMP"
        # keep only report + 1coords (delta files are large)
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
echo "[gtdb50k] slice $SLICE done: $(wc -l < "$ROWS_OUT") pairs"
