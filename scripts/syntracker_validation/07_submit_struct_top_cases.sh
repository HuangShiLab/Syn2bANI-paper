#!/usr/bin/env bash
# Submit syn2bani struct job for selected top-discordant Syntracker pairs.
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
mkdir -p "${BASE}/logs" "${BASE}/struct_top_cases"

PAIR_LIST="${1:-${BASE}/struct_top_cases/pairs_to_struct.tsv}"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=syntracker_struct
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=${BASE}/logs/struct_top_cases_%j.out
#SBATCH --error=${BASE}/logs/struct_top_cases_%j.err

/lustre1/g/aos_shihuang/data/syntracker_validation/scripts/07_run_struct_top_cases.sh "${PAIR_LIST}"
EOF

echo "Submitted struct job. Monitor with: squeue -u \$USER"
