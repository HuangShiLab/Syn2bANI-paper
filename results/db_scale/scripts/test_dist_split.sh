#!/bin/bash
B=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/db_scale
S2B=/lustre1/g/aos_shihuang/Syn2bANI-bench/target/release/syn2bani
P=$(awk -F'\t' '$3>95{print; exit}' "$B/lists/accuracy_pairs.tsv")
q=$(echo "$P" | cut -f1); r=$(echo "$P" | cut -f2); g3=$(head -1 "$B/lists/n100.txt")
echo "q=$(basename "$q") r=$(basename "$r") g3=$(basename "$g3")"
echo "--- dist q r (expect 1 row q-r)"
"$S2B" dist "$q" "$r" --enzymes BcgI,AlfI,AloI,FalI -p -t 4
echo "--- dist q r g3 (which pairs reported?)"
"$S2B" dist "$q" "$r" "$g3" --enzymes BcgI,AlfI,AloI,FalI -p -t 4
