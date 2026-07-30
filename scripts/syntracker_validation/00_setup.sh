#!/usr/bin/env bash
# Copy metadata and scripts to the HPC working tree.
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
REPO="/lustre1/g/aos_shihuang/Syn2bANI-paper"

mkdir -p "${BASE}/"{samples,manifests,references,reads,assemblies,syn2bani,skani,figures,logs,scripts}

cp "${REPO}/data/syntracker/samples_"*.tsv "${BASE}/samples/"
cp "${REPO}/data/syntracker/references.tsv" "${BASE}/samples/"
cp "${REPO}/data/syntracker/fastq_manifest_"*.tsv "${BASE}/manifests/"
cp "${REPO}/scripts/syntracker_validation/"*.sh "${BASE}/scripts/"
cp "${REPO}/scripts/syntracker_validation/"*.py "${BASE}/scripts/"

echo "Setup complete in ${BASE}"
