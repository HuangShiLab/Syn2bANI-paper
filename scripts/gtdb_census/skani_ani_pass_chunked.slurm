#!/bin/bash
#SBATCH --job-name=gtdb-ani-chunk
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/ani_chunked_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/ani_chunked_%j.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
R=$L/Syn2b-paper/results/census/gtdb
mkdir -p "$W/ani_chunk_norm" "$W/ani_chunk_tmp" "$R"
cd "$L/Syn2b-paper"
python3 scripts/gtdb_census/skani_ani_pass_chunked.py \
  --workdir "$W" --skani "$L/tools/bin/skani" \
  --output "$R/gtdb_ani_remaining_chunked.tsv.gz" \
  --norm-dir "$W/ani_chunk_norm" --tmpdir "$W/ani_chunk_tmp" \
  --workers 4 --search-threads 4 --chunk-size 5000 \
  --min-af 0 --screen 0
echo GTDB-ANI-CHUNKED-DONE
