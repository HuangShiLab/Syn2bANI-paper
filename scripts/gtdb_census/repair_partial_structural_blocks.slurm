#!/bin/bash
#SBATCH --job-name=gtdb-struct-repair
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=72:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_struct_repair_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_struct_repair_%j.err
set -euo pipefail
export RAYON_NUM_THREADS=1
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
R=$L/Syn2b-paper/results/census/gtdb
cd "$L/Syn2b-paper"
mkdir -p "$W/logs" /tmp/$SLURM_JOB_ID
TASKS="$W/ecoli_structural_repair_tasks.jsonl"
if [[ ! -s "$TASKS" ]]; then
  python3 scripts/gtdb_census/repair_partial_structural_blocks.py \
    --coverage "$R/ecoli_blockpair_coverage.tsv" \
    --cluster "Escherichia coli" \
    --workdir "$W" \
    --archive-dir "$W/structural_repair_quarantine_4259011" \
    --output "$TASKS" \
    --execute
fi
python3 scripts/gtdb_census/census_worker.py \
  --workdir "$W" \
  --tasks "$TASKS" \
  --syn2b "$L/Syn2b/target/release/syn2b" \
  --tmpdir /tmp/$SLURM_JOB_ID \
  --workers "$SLURM_CPUS_PER_TASK" --max-mega 12 \
  --soft-gb 250 --hard-gb 320 --min-free-gb 50 \
  --drain
echo GTDB-STRUCTURAL-REPAIR-DONE
