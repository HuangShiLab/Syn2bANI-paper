#!/usr/bin/env bash
# Build ANIm (dnadiff/nucmer) ground truth + Syn2bANI v8 estimates for the
# 15 mid-ANI (85-95%) oral/gut validation pairs.
# Runs on HPC2021 as a single SLURM job.
set -euo pipefail

BASE="/lustre1/g/aos_shihuang/data/validation_mid_ani"
PAIRS="${BASE}/mid_ani_pairs_85_95.tsv"
GENOMES="${BASE}/genomes"
OUT="${BASE}/anim"
SYN2BANI="${SYN2BANI:-/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate syn2bani 2>/dev/null || true
fi
command -v dnadiff >/dev/null || { echo "ERROR: dnadiff not found" >&2; exit 1; }

mkdir -p "${OUT}/dnadiff"

# --- 1. Syn2bANI v8, all-vs-all over the union of genomes -------------------
tail -n +2 "${PAIRS}" | cut -f1,2 | tr '\t' '\n' | sort -u | \
  sed "s|^|${GENOMES}/|; s|$|.fna|" > "${OUT}/genome_list.txt"
echo "[syn2bani] $(wc -l < "${OUT}/genome_list.txt") genomes, all-vs-all"
"${SYN2BANI}" ani --ql "${OUT}/genome_list.txt" --rl "${OUT}/genome_list.txt" \
  --verbose -t "${THREADS}" -p -o "${OUT}/syn2bani_v8.tsv"

# --- 2. skani on the same genome set (for a like-for-like row) ---------------
if command -v skani >/dev/null 2>&1; then
  skani dist --ql "${OUT}/genome_list.txt" --rl "${OUT}/genome_list.txt" \
    -t "${THREADS}" -o "${OUT}/skani.tsv" || echo "WARN: skani failed" >&2
fi

# --- 3. dnadiff (ANIm) per pair ----------------------------------------------
tail -n +2 "${PAIRS}" | while IFS=$'\t' read -r q r rest; do
  prefix="${OUT}/dnadiff/${q}__${r}"
  if [ -s "${prefix}.report" ]; then
    echo "[dnadiff] ${q} vs ${r} exists, skipping"
    continue
  fi
  dnadiff -p "${prefix}" "${GENOMES}/${r}.fna" "${GENOMES}/${q}.fna" >/dev/null 2>&1 \
    || echo "WARN: dnadiff failed for ${q} vs ${r}" >&2
done

# --- 4. Extract ANIm summary table -------------------------------------------
out_tsv="${OUT}/anim_truth.tsv"
printf "query\treference\tanim_ani\tanim_aligned_query\tanim_aligned_ref\n" > "${out_tsv}"
tail -n +2 "${PAIRS}" | while IFS=$'\t' read -r q r rest; do
  rep="${OUT}/dnadiff/${q}__${r}.report"
  [ -s "${rep}" ] || continue
  python3 - "${rep}" "${q}" "${r}" >> "${out_tsv}" <<'PY'
import re, sys
rep, q, r = sys.argv[1], sys.argv[2], sys.argv[3]
txt = open(rep).read()
# The first (1-to-1) block: AlignedBases and AvgIdentity lines carry
# "ref_value query_value" pairs.
def grab(label):
    m = re.search(rf"^{label}\s+(.*)$", txt, re.M)
    return m.group(1).split() if m else []
ab = grab("AlignedBases")       # e.g. ['123(45.67%)', '124(46.78%)']
ai = grab("AvgIdentity")        # e.g. ['87.12', '87.23']
try:
    aq = float(re.search(r"\(([\d.]+)%\)", ab[0]).group(1))
    ar = float(re.search(r"\(([\d.]+)%\)", ab[1]).group(1))
    ident = (float(ai[0]) + float(ai[1])) / 2.0
    print(f"{q}\t{r}\t{ident:.4f}\t{aq:.2f}\t{ar:.2f}")
except Exception:
    sys.stderr.write(f"WARN: could not parse {rep}\n")
PY
done

echo "Done. Results in ${OUT}"
