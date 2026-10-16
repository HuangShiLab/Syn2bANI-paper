#!/bin/bash
#SBATCH --job-name=gtdb-tri-convert
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-4
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/ani_triangle_convert_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/ani_triangle_convert_%A_%a.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
R=$L/Syn2b-paper/results/census/gtdb
cd "$L/Syn2b-paper"
case "$SLURM_ARRAY_TASK_ID" in
  1) SAFE="Escherichia coli"; CLUSTER="Escherichia coli" ;;
  2) SAFE="Salmonella enterica"; CLUSTER="Salmonella enterica" ;;
  3) SAFE="Klebsiella pneumoniae"; CLUSTER="Klebsiella pneumoniae" ;;
  4) SAFE="Staphylococcus aureus"; CLUSTER="Staphylococcus aureus" ;;
  *) echo "invalid array task" >&2; exit 2 ;;
esac
python3 scripts/gtdb_census/convert_skani_triangle_raw.py \
  --raw "$W/ani_triangle_tmp/$SAFE.triangle.raw" \
  --accessions "$W/effective_cluster_accessions/$SAFE.txt" \
  --cluster "$CLUSTER" \
  --output "$W/ani_triangle/$SAFE.ani.tsv.gz" \
  --report "$W/ani_triangle/$SAFE.ani.tsv.gz.report.json" \
  --pigz-threads 8
echo GTDB-TRIANGLE-CONVERT-DONE
