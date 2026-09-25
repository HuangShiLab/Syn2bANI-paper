#!/bin/bash
# Run minimap2 asm20 for all 103 struct pairs: ref vs query (same orientation as struct runs).
set -u
cd /Users/macstudio/Downloads/Syn2bANI-paper/case_studies/fda_argos
MM=/opt/homebrew/bin/minimap2
mkdir -p validate/pafs
n=0
for bed in struct_pairs99.9/*.bed; do
  base=$(basename "$bed" .bed)
  ref=${base%%__*}
  qry=${base##*__}
  out="validate/pafs/${base}.paf"
  [ -s "$out" ] && continue
  "$MM" -cx asm20 "genomes/ecoli/${ref}.fna" "genomes/ecoli/${qry}.fna" > "$out" 2> "validate/pafs/${base}.log"
  rc=$?
  n=$((n+1))
  echo "[$n] $base rc=$rc lines=$(wc -l < "$out")"
done
echo "ALL DONE: $(ls validate/pafs/*.paf | wc -l) PAFs"
