#!/bin/bash
# Run syn2bani struct on all ANI>=99.9 E. coli pairs.
# ref = alphabetically first accession, query = second.
set -u
cd /Users/macstudio/Downloads/Syn2bANI-paper/case_studies/fda_argos
SYN=/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani
mkdir -p struct_pairs99.9 struct_logs
n=0
tail -n +2 pairs99.9_ecoli.tsv | while IFS=$'\t' read -r acc1 s1 acc2 s2 ani; do
  if [[ "$acc1" < "$acc2" ]]; then ref=$acc1; qry=$acc2; else ref=$acc2; qry=$acc1; fi
  out="struct_pairs99.9/${ref}__${qry}.bed"
  log="struct_logs/${ref}__${qry}.log"
  if [ -s "$out" ] && [ -s "$log" ]; then continue; fi
  contigs=$(grep '^>' "genomes/ecoli/${ref}.fna" | sed 's/^>//; s/ .*//' | paste -sd, -)
  "$SYN" struct "genomes/ecoli/${qry}.fna" "genomes/ecoli/${ref}.fna" \
    --bed --circular "$contigs" -o "$out" > "$log" 2>&1
  rc=$?
  n=$((n+1))
  echo "[$n] $ref vs $qry rc=$rc $(tail -1 "$log" | cut -c1-100)"
done
echo "ALL DONE: $(ls struct_pairs99.9 | wc -l) outputs"
