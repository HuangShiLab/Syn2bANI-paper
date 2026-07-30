#!/usr/bin/env bash
# Submit skani comparison job.
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
mkdir -p "${BASE}/logs"

sbatch <<'EOF'
#!/bin/bash
#SBATCH --job-name=syntracker_skani
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=${BASE}/logs/skani_%j.out
#SBATCH --error=${BASE}/logs/skani_%j.err

/lustre1/g/aos_shihuang/data/syntracker_validation/scripts/05_run_skani.sh
EOF

echo "Submitted skani job. Monitor with: squeue -u \$USER"
