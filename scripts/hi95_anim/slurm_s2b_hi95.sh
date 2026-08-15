#!/bin/bash -l
#SBATCH --job-name=hi95_s2b
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work/logs/s2b_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work/logs/s2b_%j.err

set -euo pipefail

SYN2BANI=/lustre1/g/aos_shihuang/Syn2bANI-hi95/target/release/syn2bani
GDIR=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work
PAIRS=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/anim_truth_hi95.tsv
mkdir -p "${WORK}/s2b_out" "${WORK}/logs/s2b_pair"

[ -x "${SYN2BANI}" ] || { echo "ERROR: ${SYN2BANI} missing/not executable" >&2; exit 1; }

run_pair() {
  q="$1"; r="$2"
  SYN2BANI=/lustre1/g/aos_shihuang/Syn2bANI-hi95/target/release/syn2bani
  GDIR=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
  WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work
  out="${WORK}/s2b_out/${q}__${r}.tsv"
  [ -s "${out}" ] && return 0
  # default --enzymes is the standard BcgI,AlfI,AloI,FalI panel
  "${SYN2BANI}" ani "${GDIR}/${q}.fna" "${GDIR}/${r}.fna" --verbose -t 1 \
    2> "${WORK}/logs/s2b_pair/${q}__${r}.err" \
    | awk -F'\t' 'NR==1 && $1=="query" {next} {print}' > "${out}.tmp" \
    || true
  if [ -s "${out}.tmp" ]; then mv "${out}.tmp" "${out}"; else
    rm -f "${out}.tmp"; echo "${q}	${r}" >> "${WORK}/s2b_fail.txt"
  fi
}
export -f run_pair

tail -n +2 "${PAIRS}" | cut -f1,2 | xargs -P 31 -L1 bash -c 'run_pair "$@"' _

echo "fragments: $(ls "${WORK}"/s2b_out/*.tsv | wc -l) / pairs: $(($(wc -l < "${PAIRS}") - 1))"
[ -f "${WORK}/s2b_fail.txt" ] && { echo "FAILURES: $(wc -l < "${WORK}/s2b_fail.txt")"; cat "${WORK}/s2b_fail.txt"; } || echo "no failures"
