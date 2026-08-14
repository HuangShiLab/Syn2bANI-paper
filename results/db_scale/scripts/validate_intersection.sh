#!/bin/bash
# Validate the 10 search-intersection pairs with the `ani` MLE path.
B=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/db_scale
S2B=/lustre1/g/aos_shihuang/Syn2bANI-bench/target/release/syn2bani
G=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
OUT="$B/out/search_intersection_ani.tsv"
printf 'query\tref\tani_mle\tani_gated\taf_q\taf_r\n' > "$OUT"
for pair in \
  "GCA_002862645.1 GCA_009886265.1" \
  "GCA_002862745.1 GCA_002708985.1" \
  "GCA_013361685.1 GCA_013362315.1" \
  "GCA_013911625.1 GCA_011525035.1" \
  "GCA_013911625.1 GCA_012269745.1" \
  "GCA_014190035.1 GCA_000428345.2" \
  "GCA_014190155.1 GCA_013205455.1"; do
    set -- $pair
    row=$("$S2B" ani "$G/$1.fna" "$G/$2.fna" 2>/dev/null | tail -1)
    ani=$(echo "$row" | cut -f3); gated=$(echo "$row" | awk -F'\t' '{print $(NF-1)}')
    afq=$(echo "$row" | cut -f5); afr=$(echo "$row" | cut -f6)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$ani" "$gated" "$afq" "$afr" >> "$OUT"
done
cat "$OUT"
