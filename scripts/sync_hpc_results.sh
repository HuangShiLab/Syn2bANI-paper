#!/bin/bash
# Sync key results from HPC to local Mac
set -euo pipefail
HPC=shihuang@hpc2021.hku.hk
HPC_ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
LOCAL_ROOT=/Users/macstudio/Downloads/Syn2bANI-paper/results/gtdb50k

mkdir -p "$LOCAL_ROOT"

for f in \
  dnadiff_inverted_fraction.tsv \
  dnadiff_inverted_fraction_high_ani.tsv \
  dnadiff_inverted_fraction_high_ani_all.tsv \
  syn2b_inverted_fraction_50k.tsv \
  syn2b_inverted_fraction_high_ani.tsv \
  syn2b_inverted_fraction_high_ani_all.tsv \
  inverted_fraction_comparison_report.md; do
    scp "$HPC:$HPC_ROOT/$f" "$LOCAL_ROOT/$f" 2>/dev/null || echo "skip $f (not ready)"
done
