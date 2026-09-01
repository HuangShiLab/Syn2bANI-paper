#!/bin/bash
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
python3 "$ROOT/scripts/compute_dnadiff_inverted_fraction.py" \
  --pairs "$ROOT/results/gtdb50k/high_ani_pairs_ready.tsv" \
  --outdir "$ROOT/results/gtdb50k/out_high_ani" \
  --outfile "$ROOT/results/gtdb50k/dnadiff_inverted_fraction_high_ani_all.tsv" \
  --workers 16
