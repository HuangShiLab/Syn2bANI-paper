#!/bin/bash
#SBATCH --job-name=census-ani
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/ani_stream_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/ani_stream_%j.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
cd "$L/Syn2b-paper"
python3 scripts/gtdb_census/skani_ani_pass.py \
  --workdir "$L/data/syn2b_census" \
  --skani "$L/tools/bin/skani" \
  --workers 8 --large-threshold 1000 --large-workers 2 \
  --large-triangle-threads 8 --small-triangle-threads 2
echo ANI-STREAM-PASS-DONE
