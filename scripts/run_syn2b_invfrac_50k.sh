#!/bin/bash
# Run Syn2b inverted_fraction on all 43,334 held-out GTDB pairs
set -euo pipefail

ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
PAIRS=$ROOT/results/gtdb50k/pairs_50k.tsv
GENOME_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all
SYN2B=/lustre1/g/aos_shihuang/Syn2b/target/release/syn2b
TGT_DIR=$ROOT/results/gtdb50k/syn2b_tgts_cache
OUT=$ROOT/results/gtdb50k/syn2b_inverted_fraction_50k.tsv
LOG=$ROOT/results/gtdb50k/logs/syn2b_invfrac_50k.log

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
  --workers 4 \
  > "$LOG" 2>&1
