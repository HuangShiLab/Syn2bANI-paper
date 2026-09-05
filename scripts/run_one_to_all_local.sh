#!/bin/bash
set -u
OUT=/Users/macstudio/Downloads/Syn2bANI-paper/results/efficiency_v8
PY=/Users/macstudio/Downloads/Syn2bANI-paper/scripts/bench_one_to_all_local.py
RES=$OUT/one_to_all_scaling.tsv
mkdir -p "$OUT"

if [ ! -f "$RES" ]; then
    printf "tool\tmode\tn_genomes\trep\twall_s\tpeak_rss_mb\n" > "$RES"
fi

for n in 2 5 10 15 21; do
  FLIST=$OUT/lists/fasta_n${n}_one_to_all.txt
  for rep in 1 2 3; do
    python3 "$PY" "$FLIST" "$n" "$rep" "$RES"
  done
done

echo "ALL DONE"
