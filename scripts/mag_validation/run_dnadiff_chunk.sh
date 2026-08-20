#!/bin/bash
# run_dnadiff_chunk.sh <chunk_id> — run dnadiff for a chunk of pairs_anchor.tsv
# (CHUNK pairs per task, set by DNADIFF_CHUNK in common.sh). Idempotent: pairs
# with an existing parsed report are skipped. REF=majority source genome, QRY=MAG.
set -uo pipefail
source "$(dirname "$0")/common.sh"

CHUNK_ID=$1
START=$(( CHUNK_ID * DNADIFF_CHUNK + 1 ))
END=$(( (CHUNK_ID + 1) * DNADIFF_CHUNK ))
OUTDIR=$WORK/truth
mkdir -p "$OUTDIR/out" "$OUTDIR/rows"
ROWS_TMP="$OUTDIR/rows/chunk_${CHUNK_ID}.tsv.tmp"
ROWS_OUT="$OUTDIR/rows/chunk_${CHUNK_ID}.tsv"

export PATH="$ANVIO_BIN:$PATH"
: > "$ROWS_TMP"
tail -n +2 "$WORK/pairs/pairs_anchor.tsv" | sed -n "${START},${END}p" | while IFS=$'\t' read -r BINID REF; do
    BINFA=$WORK/bins_all/${BINID}.fa
    DD=$OUTDIR/out/$BINID
    if [ ! -s "$DD/dd.report" ]; then
        # dnadiff (mummer 3.23) aborts on ANY whitespace in a FASTA (headers
        # AND sequence lines); CAMI2 source genomes have both. Use a fully
        # whitespace-stripped cached copy for those refs.
        REFUSE=$REF
        if grep -qm1 -P "\s" "$REF"; then
            mkdir -p "$OUTDIR/refs_san"
            KEY=$(basename "$REF" | tr -c 'A-Za-z0-9._-' '_')
            REFUSE=$OUTDIR/refs_san/$KEY
            [ -s "$REFUSE" ] || awk '/^>/{sub(/[ \t].*$/, ""); print; next} {gsub(/[ \t]/, ""); print}' "$REF" > "$REFUSE"
        fi
        mkdir -p "$DD"
        (cd "$DD" && dnadiff -p dd "$REFUSE" "$BINFA" > /dev/null 2>&1) \
            || echo "[truth] $BINID: dnadiff failed"
    fi
    if [ -s "$DD/dd.report" ]; then
        $PY "$SCRIPTS/parse_dnadiff.py" "$BINID" "$DD/dd.report" >> "$ROWS_TMP"
    else
        echo -e "${BINID}\tNA\tNA\tNA\tNA\tdnadiff_failed" >> "$ROWS_TMP"
    fi
done
mv "$ROWS_TMP" "$ROWS_OUT"
echo "[truth] chunk $CHUNK_ID done: $(wc -l < "$ROWS_OUT") pairs"
