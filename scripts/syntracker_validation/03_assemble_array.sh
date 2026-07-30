#!/usr/bin/env bash
# SLURM array driver for assembling all SynTracker isolate FASTQs.
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
SAMPLE_DIR="${BASE}/samples"
ASM_DIR="${BASE}/assemblies"
LIST="${BASE}/manifests/all_isolates.tsv"
mkdir -p "${ASM_DIR}" "${BASE}/manifests"

# Build one task list from the four sample sheets
> "${LIST}"
for tsv in "${SAMPLE_DIR}"/samples_*.tsv; do
  [ -e "${tsv}" ] || continue
  tail -n +2 "${tsv}" | awk -F'\t' -v OFS='\t' '{print $2, $1}' >> "${LIST}"
done

N=$(wc -l < "${LIST}")
echo "Prepared ${N} assembly tasks -> ${LIST}"

# Submit array job
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=syntracker_asm
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-${N}%20
#SBATCH --output=${BASE}/logs/asm_%A_%a.out
#SBATCH --error=${BASE}/logs/asm_%A_%a.err

set -euo pipefail
LINE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${LIST}")
ISOLATE=\$(echo "\${LINE}" | cut -f1)
SPECIES=\$(echo "\${LINE}" | cut -f2)

"${BASE}/scripts/03_assemble_one.sh" "\${ISOLATE}" "${BASE}/reads" "${ASM_DIR}"
EOF

echo "Submitted assembly array. Monitor with: squeue -u \$USER"
