#!/bin/bash
#SBATCH --job-name=hrom-repair
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2b-paper/results/census/hrom/quality_stratified/repair_screen_%j.log
set -euo pipefail
E=/lustre1/g/aos_shihuang/databases/HROM/data/protein_catalog/HROP-100_eggnog.mapper.tsv
C=/tmp/hrom_ani_candidate_clusters.txt
O=/lustre1/g/aos_shihuang/Syn2b-paper/results/census/hrom/quality_stratified/hrom_ani_enriched_repair_genes.tsv
awk -F '\t' 'NR==FNR {p[$1]=1;next} FNR==1 {print;next} {for(x in p) if(index($1,x)==1){ if($0 ~ /[Mm]ut[SsLlHhYyTtMm]|dnaQ|uvrD|ogt|ada|alkB|mutM|mutY|mutT|MMR|mismatch repair/) print; break}}' "$C" "$E" > "$O"
echo REPAIR_SCREEN_DONE
