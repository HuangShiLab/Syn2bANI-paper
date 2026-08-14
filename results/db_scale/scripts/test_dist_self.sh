#!/bin/bash
B=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/db_scale
S2B=/lustre1/g/aos_shihuang/Syn2bANI-bench/target/release/syn2bani
P=$(awk -F'\t' '$3>95{print; exit}' "$B/lists/accuracy_pairs.tsv")
q=$(echo "$P" | cut -f1); r=$(echo "$P" | cut -f2)
echo "=== dist q q (self, ENZ4, min-af 0)"
"$S2B" dist "$q" "$q" --enzymes BcgI,AlfI,AloI,FalI --min-af 0 -p -t 4 2>/dev/null
echo "=== three more 95-100 pairs: dist ENZ4 min-af 0"
awk -F'\t' '$3>95{c++; if(c>1 && c<=4) print}' "$B/lists/accuracy_pairs.tsv" | while IFS=$'\t' read -r qq rr tt; do
    echo "--- $(basename "$qq" .fna) vs $(basename "$rr" .fna) truth=$tt"
    "$S2B" dist "$qq" "$rr" --enzymes BcgI,AlfI,AloI,FalI --min-af 0 -p -t 4 2>/dev/null | tail -1
done
echo "=== same three pairs: ani rerun"
awk -F'\t' '$3>95{c++; if(c>1 && c<=4) print}' "$B/lists/accuracy_pairs.tsv" | while IFS=$'\t' read -r qq rr tt; do
    "$S2B" ani "$qq" "$rr" 2>/dev/null | tail -1 | cut -f1-3
done
