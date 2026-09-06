#!/bin/bash
#SBATCH --job-name=s2b_fmh1582
#SBATCH --partition=intel
#SBATCH --array=0-49
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/s2b_fmh1582_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/s2b_fmh1582_%A_%a.err
set -euo pipefail

ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
GENOMES=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
SYN2B=/lustre1/g/aos_shihuang/Syn2b/target/release/syn2b
NSLICES=50

mkdir -p "$WORK/logs" "$WORK/syn2b_tgts_cache_fmh1582" "$WORK/syn2b_inverted_fraction_50k_fmh1582_parts"

SLICE=$SLURM_ARRAY_TASK_ID
OUT="$WORK/syn2b_inverted_fraction_50k_fmh1582_parts/part_${SLICE}.tsv"
[ -s "$OUT" ] && exit 0

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate checkm2

python3 "$ROOT/scripts/run_syn2b_inverted_fraction.py" \
    --pairs      "$WORK/pairs_50k.tsv" \
    --genome-dir "$GENOMES" \
    --syn2b      "$SYN2B" \
    --mode       fracminhash \
    --kmer       31 \
    --scale      1582 \
    --tgt-dir    "$WORK/syn2b_tgts_cache_fmh1582" \
    --out        "$OUT" \
    --slice-id   "$SLICE" \
    --n-slices   "$NSLICES" \
    --workers    2

echo "[fmh1582] slice $SLICE done"
