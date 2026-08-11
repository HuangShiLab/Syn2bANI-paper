#!/usr/bin/env bash
# Run skani dist on all isolate assemblies per species.
set -euo pipefail

# skani is expected in the 'syn2bani' conda environment.
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate syn2bani 2>/dev/null || true
fi

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
ASM_DIR="${BASE}/assemblies"
OUT_DIR="${BASE}/skani"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
mkdir -p "${OUT_DIR}"

if ! command -v skani >/dev/null 2>&1; then
  echo "ERROR: skani not found. Install it first (see scripts/syntracker_validation/install_skani.sh)." >&2
  exit 1
fi

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
  skani dist --ql "${genomes_list}" --rl "${genomes_list}" -t "${THREADS}" -o "${out_tsv}"
done

echo "skani results in ${OUT_DIR}"
