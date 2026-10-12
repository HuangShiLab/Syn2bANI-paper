#!/bin/bash
#SBATCH --job-name=gtdb-triangle
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --array=1-4
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/ani_triangle_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/ani_triangle_%A_%a.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
R=$L/Syn2b-paper/results/census/gtdb
case "$SLURM_ARRAY_TASK_ID" in
  1) CLUSTER="Escherichia coli" ;;
  2) CLUSTER="Salmonella enterica" ;;
  3) CLUSTER="Klebsiella pneumoniae" ;;
  4) CLUSTER="Staphylococcus aureus" ;;
  *) echo "invalid array task" >&2; exit 2 ;;
esac
SAFE=${CLUSTER//\//_}
mkdir -p "$W/ani_triangle" "$W/ani_triangle_tmp" "$W/logs"
cd "$L/Syn2b-paper"
python3 scripts/gtdb_census/skani_ani_triangle_cluster.py \
  --workdir "$W" \
  --skani "$L/tools/bin/skani" \
  --cluster "$CLUSTER" \
  --output "$W/ani_triangle/$SAFE.ani.tsv.gz" \
  --tmpdir "$W/ani_triangle_tmp" \
  --threads 60 \
  --pigz-threads 4 \
  --min-af 0 \
  --screen 0
