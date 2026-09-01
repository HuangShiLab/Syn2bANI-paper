#!/bin/bash
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
PAIRS=$ROOT/results/gtdb50k/high_ani_pairs_ready.tsv
GENOME_DIR=$ROOT/results/gtdb50k/genomes_high_ani
SYN2B=/lustre1/g/aos_shihuang/Syn2b/target/release/syn2b
TGT_DIR=$ROOT/results/gtdb50k/syn2b_tgts_cache_high_ani
OUT=$ROOT/results/gtdb50k/syn2b_inverted_fraction_high_ani_all.tsv
LOG=$ROOT/results/gtdb50k/logs/syn2b_invfrac_high_ani_all.log

mkdir -p "$TGT_DIR"
mkdir -p "$ROOT/results/gtdb50k/logs"

cd "$ROOT"
export RAYON_NUM_THREADS=1
python3 "$ROOT/scripts/run_syn2b_inverted_fraction.py" \
  --pairs "$PAIRS" \
  --genome-dir "$GENOME_DIR" \
  --syn2b "$SYN2B" \
  --tgt-dir "$TGT_DIR" \
  --out "$OUT" \
  --workers 2 \
  > "$LOG" 2>&1
