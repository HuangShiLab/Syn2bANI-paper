#!/bin/bash
# Search: 100 held-out queries vs 5000-genome sketch DB. syn2bani search vs skani search.
# Requires sketches/s2ba_n5000 and sketches/skani_n5000 (from bench_sketch.sh 5000).
# Usage: bench_search.sh <reps>
source "$(dirname "$0")/common.sh"
REPS=${1:-1}
TSV="$BASE/search_scaling.tsv"
[ -s "$TSV" ] || printf 'tool\tphase\tn\trep\twall_s\tmax_rss_kb\texit_status\ttimestamp_utc\tout_rows\n' > "$TSV"
QLIST="$BASE/lists/queries100.txt"
S2B_DB="$BASE/sketches/s2ba_n5000"
SKANI_DB="$BASE/sketches/skani_n5000"
[ -d "$S2B_DB" ] || { echo "missing $S2B_DB"; exit 1; }
[ -d "$SKANI_DB" ] || { echo "missing $SKANI_DB"; exit 1; }

mapfile -t Q < "$QLIST"
echo "[$(date -u +%H:%M:%S)] search bench: 100 queries vs n=5000 DB, reps=$REPS threads=$THREADS"
for rep in $(seq 1 "$REPS"); do
    OUT_S="$BASE/out/search_s2b_r${rep}.tsv"
    run_timed "$TSV.tmp" syn2bani search "100x5000" "$rep" "$BASE/logs/search_s2b_r${rep}.log" \
        "$S2B" search "${Q[@]}" "$S2B_DB" -o "$OUT_S" -p -t "$THREADS" --min-ani 0.8 || true
    rows=$(($(wc -l < "$OUT_S" 2>/dev/null || echo 1) - 1))
    sed -i "\$ s/\$/\t$rows/" "$TSV.tmp"

    OUT_K="$BASE/out/search_skani_r${rep}.tsv"
    run_timed "$TSV.tmp" skani search "100x5000" "$rep" "$BASE/logs/search_skani_r${rep}.log" \
        "$SKANI" search --ql "$QLIST" -d "$SKANI_DB" -t "$THREADS" -o "$OUT_K" || true
    rows=$(($(wc -l < "$OUT_K" 2>/dev/null || echo 1) - 1))
    sed -i "\$ s/\$/\t$rows/" "$TSV.tmp"
done
cat "$TSV.tmp" >> "$TSV"; rm -f "$TSV.tmp"
echo "[$(date -u +%H:%M:%S)] search bench done"
