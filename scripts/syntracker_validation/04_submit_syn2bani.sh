#!/usr/bin/env bash
# Submit Syn2bANI comparison job for all four species.
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
mkdir -p "${BASE}/logs"

sbatch <<'EOF'
#!/bin/bash
#SBATCH --job-name=syntracker_syn2bani
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=${BASE}/logs/syn2bani_%j.out
#SBATCH --error=${BASE}/logs/syn2bani_%j.err

/lustre1/g/aos_shihuang/data/syntracker_validation/scripts/04_run_syn2bani.sh
EOF

echo "Submitted Syn2bANI job. Monitor with: squeue -u \$USER"
