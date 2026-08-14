#!/bin/bash
# Sketch scaling: syn2bani sketch (4-enzyme panel) vs skani sketch.
# Usage: bench_sketch.sh <n> <reps>
source "$(dirname "$0")/common.sh"
N=$1; REPS=$2
LIST="$BASE/lists/n${N}.txt"
TSV="$BASE/sketch_scaling.tsv"
[ -s "$TSV" ] || printf 'tool\tphase\tn\trep\twall_s\tmax_rss_kb\texit_status\ttimestamp_utc\tstore_bytes\n' > "$TSV"

mapfile -t FILES < "$LIST"
echo "[$(date -u +%H:%M:%S)] sketch bench n=$N reps=$REPS threads=$THREADS"
for rep in $(seq 1 "$REPS"); do
    rm -rf "$BASE/sketches/s2ba_n${N}"
    run_timed "$TSV.tmp" syn2bani sketch "$N" "$rep" "$BASE/logs/sketch_s2b_n${N}_r${rep}.log" \
        "$S2B" sketch "${FILES[@]}" -o "$BASE/sketches/s2ba_n${N}" --enzymes "$ENZ4" -p -t "$THREADS" || true
    sz=$(du -sb "$BASE/sketches/s2ba_n${N}" | cut -f1)
    sed -i "\$ s/\$/\t$sz/" "$TSV.tmp"

    rm -rf "$BASE/sketches/skani_n${N}"
    run_timed "$TSV.tmp" skani sketch "$N" "$rep" "$BASE/logs/sketch_skani_n${N}_r${rep}.log" \
        "$SKANI" sketch -l "$LIST" -o "$BASE/sketches/skani_n${N}" -t "$THREADS" || true
    sz=$(du -sb "$BASE/sketches/skani_n${N}" | cut -f1)
    sed -i "\$ s/\$/\t$sz/" "$TSV.tmp"
done
cat "$TSV.tmp" >> "$TSV"; rm -f "$TSV.tmp"
# keep the largest stores for downstream phases; remove smaller ones to save quota
if [ "$N" -lt 5000 ]; then rm -rf "$BASE/sketches/s2ba_n${N}" "$BASE/sketches/skani_n${N}"; fi
echo "[$(date -u +%H:%M:%S)] sketch bench n=$N done"
