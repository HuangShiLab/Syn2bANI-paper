#!/bin/bash
#SBATCH --job-name=gtdb-qual
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_quality_sensitivity_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_quality_sensitivity_%j.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
cd "$L/Syn2b-paper"
bash scripts/gtdb_census/gtdb_quality_sensitivity.sh
