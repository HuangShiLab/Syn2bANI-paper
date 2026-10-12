#!/bin/bash
set -euo pipefail
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
cd "$L/Syn2b-paper"
mkdir -p "$W/logs"
# The old merge dependency referred to the cancelled 16-core search job.
scancel 4258680 2>/dev/null || true
python3 scripts/gtdb_census/promote_gtdb_ani_layers.py \
  --workdir "$W" \
  --output "$W/ani_layer_promotion_report.json"
CHUNK=$(sbatch --parsable scripts/gtdb_census/skani_ani_pass_chunked.slurm)
echo "resume_chunk_job=$CHUNK"
MERGE=$(sbatch --parsable --dependency=afterok:"$CHUNK" scripts/gtdb_census/merge_gtdb_ani.sh)
echo "merge_job=$MERGE"
SUMMARY=$(sbatch --parsable --dependency=afterok:"$MERGE" scripts/gtdb_census/gtdb_ani_summary_join.slurm)
echo "summary_join_job=$SUMMARY"
