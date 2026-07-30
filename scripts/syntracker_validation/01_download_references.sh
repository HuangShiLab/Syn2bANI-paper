#!/usr/bin/env bash
# Download reference genomes used by SynTracker (Fig. 3).
# Run on the I/O node or login node; light I/O task.
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
OUTDIR="${BASE}/references"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Accession -> desired filename
# E. coli K-12 MG1655 is GCF_000005845.2 (SynTracker paper: NC_000913.3)
# H. pylori reference CP032479.1 corresponds to assembly GCF_006337385.1
declare -A ACC=(
  [Neisseria_gonorrhoeae]=GCF_900087635.2
  [Escherichia_coli]=GCF_000005845.2
  [Helicobacter_pylori]=GCF_006337385.1
  [Streptomyces_rimosus]=GCF_000331185.2
)

module load ncbi-datasets-cli 2>/dev/null || true

for sp in "${!ACC[@]}"; do
  acc="${ACC[$sp]}"
  echo "[ref] ${sp} : ${acc}"
  rm -rf "tmp_${sp}"
  datasets download genome accession "${acc}" --filename "${sp}.zip"
  unzip -q "${sp}.zip" -d "tmp_${sp}"
  fna=$(find "tmp_${sp}/ncbi_dataset/data/${acc}" -maxdepth 1 -name '*_genomic.fna' | head -1)
  cp "${fna}" "${OUTDIR}/${sp}.fna"
  rm -rf "tmp_${sp}" "${sp}.zip"
done

echo "References ready in ${OUTDIR}"
ls -lh "${OUTDIR}"
