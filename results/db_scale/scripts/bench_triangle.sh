#!/bin/bash
# All-vs-all: syn2bani triangle (edge list, BcgI-only hardcoded) vs skani triangle (full matrix).
# Usage: bench_triangle.sh <n> <reps>
source "$(dirname "$0")/common.sh"
N=$1; REPS=$2
LIST="$BASE/lists/n${N}.txt"
TSV="$BASE/triangle_scaling.tsv"
[ -s "$TSV" ] || printf 'tool\tphase\tn\trep\twall_s\tmax_rss_kb\texit_status\ttimestamp_utc\tout_rows\n' > "$TSV"

mapfile -t FILES < "$LIST"
echo "[$(date -u +%H:%M:%S)] triangle bench n=$N reps=$REPS threads=$THREADS"
for rep in $(seq 1 "$REPS"); do
    OUT_S="$BASE/out/triangle_s2b_n${N}_r${rep}.tsv"
    run_timed "$TSV.tmp" syn2bani triangle "$N" "$rep" "$BASE/logs/triangle_s2b_n${N}_r${rep}.log" \
        "$S2B" triangle "${FILES[@]}" --edge-list -p -t "$THREADS" -o "$OUT_S" || true
    rows=$(($(wc -l < "$OUT_S" 2>/dev/null || echo 1) - 1))
    sed -i "\$ s/\$/\t$rows/" "$TSV.tmp"

    OUT_K="$BASE/out/triangle_skani_n${N}_r${rep}.tsv"
    run_timed "$TSV.tmp" skani triangle "$N" "$rep" "$BASE/logs/triangle_skani_n${N}_r${rep}.log" \
        "$SKANI" triangle -l "$LIST" -t "$THREADS" -o "$OUT_K" || true
    rows=$(($(wc -l < "$OUT_K" 2>/dev/null || echo 1) - 1))
    sed -i "\$ s/\$/\t$rows/" "$TSV.tmp"
done
cat "$TSV.tmp" >> "$TSV"; rm -f "$TSV.tmp"
echo "[$(date -u +%H:%M:%S)] triangle bench n=$N done"
