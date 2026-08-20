#!/bin/bash
# run_fasttools_bin.sh <binid> — run syn2bani dist / skani dist / fastANI for one
# bin against its anchor + <=3 same-species alternates + nearest GTDB r207 rep.
# Outputs: $WORK/fast/per_pair/{binid}.{syn2bani,skani,fastani}.tsv and {binid}.refs.tsv
set -uo pipefail
source "$(dirname "$0")/common.sh"

BINID=$1
BINFA=$WORK/bins_all/${BINID}.fa
PP=$WORK/fast/per_pair
RLD=$WORK/fast/rl
mkdir -p "$PP" "$RLD"

[ -s "$PP/${BINID}.syn2bani.tsv" ] && [ -s "$PP/${BINID}.skani.tsv" ] && [ -s "$PP/${BINID}.fastani.tsv" ] && exit 0

# ref list with roles (anchor, alt1..3 from pairs.tsv; gtdb_rep from repsearch)
awk -F'\t' -v b="$BINID" '$1==b {print $3"\t"$4}' "$WORK/pairs/pairs.tsv" > "$RLD/${BINID}.refs.tsv"
grep -hP "^${BINID}\t" "$WORK"/repsearch/*/map.tsv 2>/dev/null | head -1 \
    | awk -F'\t' '{print $2"\tgtdb_rep"}' >> "$RLD/${BINID}.refs.tsv"
cut -f1 "$RLD/${BINID}.refs.tsv" > "$RLD/${BINID}.rl.txt"
printf '%s\n' "$BINFA" > "$RLD/${BINID}.ql.txt"
NREF=$(wc -l < "$RLD/${BINID}.rl.txt")
[ "$NREF" -eq 0 ] && { echo "[fasttools] $BINID: no refs, skip"; exit 0; }

if [ ! -s "$PP/${BINID}.syn2bani.tsv" ]; then
    "$S2B" dist --ql "$RLD/${BINID}.ql.txt" --rl "$RLD/${BINID}.rl.txt" \
        --verbose -t 1 -o "$PP/${BINID}.syn2bani.tsv" 2> "$PP/${BINID}.syn2bani.err" \
        || echo "[fasttools] $BINID: syn2bani failed"
fi
if [ ! -s "$PP/${BINID}.skani.tsv" ]; then
    "$SKANI" dist --ql "$RLD/${BINID}.ql.txt" --rl "$RLD/${BINID}.rl.txt" \
        --min-af 0 -t 1 -o "$PP/${BINID}.skani.tsv" 2> "$PP/${BINID}.skani.err" \
        || echo "[fasttools] $BINID: skani failed"
fi
if [ ! -s "$PP/${BINID}.fastani.tsv" ]; then
    "$FASTANI" -q "$BINFA" --rl "$RLD/${BINID}.rl.txt" -t 1 \
        -o "$PP/${BINID}.fastani.tsv" 2> "$PP/${BINID}.fastani.err" \
        || echo "[fasttools] $BINID: fastani failed"
fi
