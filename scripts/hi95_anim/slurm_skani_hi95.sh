#!/bin/bash -l
#SBATCH --job-name=hi95_skani
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work/logs/skani_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work/logs/skani_%j.err

set -euo pipefail

SKANI=/group/aos_shihuang/conda/envs/gtdbtk310/bin/skani
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work
mkdir -p "${WORK}/skani_out" "${WORK}/logs"

run_batch() {
  list="$1"
  b=$(basename "${list}" .list)
  out="/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work/skani_out/${b}.tsv"
  [ -s "${out}" ] && { echo "[skip] ${b}"; return 0; }
  /group/aos_shihuang/conda/envs/gtdbtk310/bin/skani dist \
    --ql "${list}" --rl "${list}" -t 4 --min-af 15 -o "${out}" \
    2> "/lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work/logs/${b}.err" \
    || echo "FAIL ${b}" >> /lustre1/g/aos_shihuang/Syn2bANI-paper/results/hi95_work/skani_fail.txt
}
export -f run_batch

ls "${WORK}"/batches/*.list | xargs -P 8 -I{} bash -c 'run_batch "$@"' _ {}

echo "batches done: $(ls "${WORK}"/skani_out/*.tsv | wc -l)"
[ -f "${WORK}/skani_fail.txt" ] && { echo "FAILURES:"; cat "${WORK}/skani_fail.txt"; } || echo "no batch failures"
