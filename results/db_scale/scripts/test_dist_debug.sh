#!/bin/bash
B=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/db_scale
S2B=/lustre1/g/aos_shihuang/Syn2bANI-bench/target/release/syn2bani
P=$(awk -F'\t' '$3>95{print; exit}' "$B/lists/accuracy_pairs.tsv")
q=$(echo "$P" | cut -f1); r=$(echo "$P" | cut -f2)
echo "q=$(basename "$q") r=$(basename "$r") truth=$(echo "$P" | cut -f3)"
echo "--- ani (validated)"
"$S2B" ani "$q" "$r"
echo "--- dist ENZ4 --min-af 0"
"$S2B" dist "$q" "$r" --enzymes BcgI,AlfI,AloI,FalI --min-af 0 -p -t 4
echo "--- dist ENZ4 --min-af 0 --mash-ani"
"$S2B" dist "$q" "$r" --enzymes BcgI,AlfI,AloI,FalI --min-af 0 --mash-ani -p -t 4
echo "--- dist default (AloI,BslFI) --min-af 0"
"$S2B" dist "$q" "$r" --min-af 0 -p -t 4
echo "--- dist single -e BcgI --min-af 0"
"$S2B" dist "$q" "$r" -e BcgI --min-af 0 -p -t 4
