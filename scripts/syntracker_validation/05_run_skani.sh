#!/usr/bin/env bash
# Run skani dist on all isolate assemblies per species.
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
ASM_DIR="${BASE}/assemblies"
OUT_DIR="${BASE}/skani"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
mkdir -p "${OUT_DIR}"

SPECIES_LIST=(
  Neisseria_gonorrhoeae
  Escherichia_coli_hypermutator
  Helicobacter_pylori
  Streptomyces_rimosus
)

for sp in "${SPECIES_LIST[@]}"; do
  sample_tsv="${BASE}/samples/samples_${sp}.tsv"
  genomes_list="${OUT_DIR}/${sp}_genomes.txt"
  tail -n +2 "${sample_tsv}" | while IFS=$'\t' read -r species isolate rest; do
    asm="${ASM_DIR}/${isolate}.fna"
    [ -s "${asm}" ] && echo "${asm}"
  done > "${genomes_list}"

  [ -s "${genomes_list}" ] || continue
  out_tsv="${OUT_DIR}/skani_${sp}.tsv"
  echo "[skani] ${sp}: $(wc -l < "${genomes_list}") genomes -> ${out_tsv}"
  skani dist -q "${genomes_list}" -r "${genomes_list}" -t "${THREADS}" -o "${out_tsv}"
done

echo "skani results in ${OUT_DIR}"
