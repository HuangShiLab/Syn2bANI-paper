#!/bin/bash
#SBATCH --job-name=gtdb-ani-summary
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_ani_summary_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_ani_summary_%j.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
R=$L/Syn2b-paper/results/census/gtdb
cd "$L/Syn2b-paper"
python3 scripts/gtdb_census/summarize_gtdb_ani_by_species.py \
  --ani "$R/gtdb_r207_ani_merged_sorted.tsv.gz" \
  --workdir "$W" \
  --quality "$W/gtdb_genome_quality.tsv" \
  --output "$R/gtdb_species_ani_summary.tsv" \
  --report "$R/gtdb_ani_merge_validation.json"
python3 scripts/gtdb_census/join_gtdb_species_tables.py \
  --structural "$R/effective_repaired/effective_species_summary.tsv" \
  --ani "$R/gtdb_species_ani_summary.tsv" \
  --output "$R/gtdb_species_ani_sv_joint.tsv"
echo GTDB-ANI-SUMMARY-JOIN-DONE
