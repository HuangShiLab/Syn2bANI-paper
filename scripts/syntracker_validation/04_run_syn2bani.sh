#!/usr/bin/env bash
# Run Syn2bANI ani (default AloI,BslFI panel) on all isolate assemblies per species.
# Submit as a SLURM job.
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
ASM_DIR="${BASE}/assemblies"
OUT_DIR="${BASE}/syn2bani"
SYN2BANI="${SYN2BANI:-/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani}"
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
    if [ ! -s "${asm}" ]; then
      echo "WARN: missing assembly ${asm}" >&2
      continue
    fi
    echo "${asm}"
  done > "${genomes_list}"

  if [ ! -s "${genomes_list}" ]; then
    echo "WARN: no assemblies for ${sp}, skipping"
    continue
  fi

  out_tsv="${OUT_DIR}/syn2bani_${sp}.tsv"
  echo "[syn2bani] ${sp}: $(wc -l < "${genomes_list}") genomes -> ${out_tsv}"
  "${SYN2BANI}" ani --ql "${genomes_list}" --rl "${genomes_list}" \
    --calibrate -t "${THREADS}" -p -o "${out_tsv}"
done

echo "Syn2bANI results in ${OUT_DIR}"
