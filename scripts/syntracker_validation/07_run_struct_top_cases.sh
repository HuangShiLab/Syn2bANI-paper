#!/usr/bin/env bash
# Run syn2bani struct on selected top-discordant Syntracker pairs.
# Input: a TSV with columns species, query, reference (no header).
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/syntracker_validation"
ASM_DIR="${BASE}/assemblies"
OUT_DIR="${BASE}/struct_top_cases"
SYN2BANI="${SYN2BANI:-/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

PAIR_LIST="${1:-${BASE}/struct_top_cases/pairs_to_struct.tsv}"

mkdir -p "${OUT_DIR}"

while IFS=$'\t' read -r species query reference; do
  q_asm="${ASM_DIR}/${query}.fna"
  r_asm="${ASM_DIR}/${reference}.fna"
  if [ ! -s "${q_asm}" ]; then
    echo "WARN: missing query assembly ${q_asm}" >&2
    continue
  fi
  if [ ! -s "${r_asm}" ]; then
    echo "WARN: missing reference assembly ${r_asm}" >&2
    continue
  fi

  case_name="${species}_${query}_vs_${reference}"
  out_paf="${OUT_DIR}/${case_name}.tsv"
  echo "[struct] ${case_name}"
  "${SYN2BANI}" struct -q "${q_asm}" -r "${r_asm}" -t "${THREADS}" -p -o "${out_paf}"
done < "${PAIR_LIST}"

echo "Struct results in ${OUT_DIR}"
