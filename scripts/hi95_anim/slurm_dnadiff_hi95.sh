#!/bin/bash -l
#SBATCH --job-name=hi95_dnadiff
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work/logs/dnadiff_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work/logs/dnadiff_%j.err

set -euo pipefail

export PATH=/group/aos_shihuang/conda/envs/anvio/bin:$PATH
GDIR=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work
mkdir -p "${WORK}/dnadiff" "${WORK}/logs/dnadiff_pair"
command -v dnadiff >/dev/null || { echo "ERROR: dnadiff not in PATH" >&2; exit 1; }

run_pair() {
  q="$1"; r="$2"
  GDIR=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
  WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work
  prefix="${WORK}/dnadiff/${q}__${r}"
  [ -s "${prefix}.report" ] && return 0
  dnadiff -p "${prefix}" "${GDIR}/${r}.fna" "${GDIR}/${q}.fna" \
    > "${WORK}/logs/dnadiff_pair/${q}__${r}.log" 2>&1 \
    || echo "${q}	${r}" >> "${WORK}/dnadiff_fail.txt"
}
export -f run_pair

tail -n +2 "${WORK}/pairs_dnadiff.tsv" | cut -f1,2 | \
  xargs -P 31 -L1 bash -c 'run_pair "$@"' _

echo "reports: $(ls "${WORK}"/dnadiff/*.report | wc -l) / pairs: $(($(wc -l < "${WORK}/pairs_dnadiff.tsv") - 1))"
[ -f "${WORK}/dnadiff_fail.txt" ] && { echo "FAILURES: $(wc -l < "${WORK}/dnadiff_fail.txt")"; } || echo "no dnadiff failures"
