#!/bin/bash
set -euo pipefail
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
S=$L/Syn2b-paper/results/census/gtdb
TABLE=$S/gtdb_r207_effective_normalized_dedup.tsv.gz
QUAL=$W/gtdb_genome_quality.tsv
OUT=$S/effective_quality_sensitivity.tsv
gzip -cd "$TABLE" | awk -v QUAL="$QUAL" '
BEGIN { FS="\t"; OFS="\t" }
NR==FNR {
  if (FNR==1) next
  comp[$1]=$4+0; cont[$1]=$5+0; n50[$1]=$7+0
  next
}
FNR==1 { next }
{
  a=$1; b=$2; bp=$3+0
  ca=(a in comp ? comp[a] : "NA"); cb=(b in comp ? comp[b] : "NA")
  xa=(a in cont ? cont[a] : "NA"); xb=(b in cont ? cont[b] : "NA")
  na=(a in n50 ? n50[a] : "NA"); nb=(b in n50 ? n50[b] : "NA")
  mc="NA"; xc="NA"; mn="NA"
  if (ca!="NA" && cb!="NA") mc=(ca<cb?ca:cb)
  if (xa!="NA" && xb!="NA") xc=(xa>xb?xa:xb)
  if (na!="NA" && nb!="NA") mn=(na<nb?na:nb)
  add("all",bp)
  if (mc!="NA" && mc>=90) add("min_completeness_ge90",bp)
  if (mc!="NA" && mc<70) add("min_completeness_lt70",bp)
  if (mc!="NA" && mc>=70 && mc<90) add("min_completeness_70_90",bp)
  if (xc!="NA" && xc<=5) add("max_contamination_le5",bp)
  if (xc!="NA" && xc>5) add("max_contamination_gt5",bp)
  if (mc!="NA" && mc>=90 && xc!="NA" && xc<=5) add("HQ_comp90_contam5",bp)
  if (mn!="NA" && mn>=10000) add("min_N50_ge10kb",bp)
  if (mn!="NA" && mn>=100000) add("min_N50_ge100kb",bp)
  if (mc!="NA" && mc>=90 && xc!="NA" && xc<=5 && mn!="NA" && mn>=10000) add("HQ_N50_ge10kb",bp)
}
function add(s,bp) { pairs[s]++; if (bp>=2) ge2[s]++ }
END {
  print "stratum\tpairs\tpairs_ge2\tfraction_ge2"
  n=asorti(pairs,idx)
  for(i=1;i<=n;i++) { s=idx[i]; print s,pairs[s],ge2[s],ge2[s]/pairs[s] }
}
' > "$OUT"
echo GTDB-QUALITY-SENSITIVITY-DONE
