#!/bin/bash
#SBATCH --job-name=gtdb-struct-repair2
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=36:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_struct_repair_pass2_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_struct_repair_pass2_%j.err
set -euo pipefail
export RAYON_NUM_THREADS=1
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
cd "$L/Syn2b-paper"
mkdir -p "$W/logs" /tmp/$SLURM_JOB_ID
python3 scripts/gtdb_census/census_worker.py \
  --workdir "$W" \
  --tasks "$W/ecoli_structural_repair_pass2_tasks.jsonl" \
  --syn2b "$L/Syn2b/target/release/syn2b" \
  --tmpdir /tmp/$SLURM_JOB_ID \
  --workers "$SLURM_CPUS_PER_TASK" --max-mega 12 \
  --soft-gb 250 --hard-gb 320 --min-free-gb 50 \
  --drain
echo GTDB-STRUCTURAL-REPAIR-PASS2-DONE
