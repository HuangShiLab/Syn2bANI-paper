# Calibration provenance for the benchmark tables

## Two models

- **v5** (`results/panel_by_band/linear_cal_v5.json`) is the deployed model in
  the Syn2bANI 0.1.0 software snapshot. It was trained on 2,520 ANIm-truth
  finite pairs (2,053 usable pairs from the nominal 2,074-pair post-rescue set
  plus the 467-pair 95--99.5% expansion); its 2,074-pair band-holdout CV
  result is MAE 0.731. The other 21 source pairs lacked finite training
  features and were excluded before training.
- **v6** (`results/panel_by_band/linear_cal_v6.json`) is the benchmark model
  used by `ani_cal` in `results/gtdb50k/s2b_50k.tsv`,
  `results/gtdb50k/high_ani_results.tsv`, the GTDB held-out benchmark, and the
  unified benchmark. It was trained on the 2,520 v5 pairs plus 1,615 high-ANI
  training pairs (`training_n = 4,134`, after one of 1,615 sampled pairs was
  dropped as non-finite), then evaluated on the 43,334 representative-genome hold-out and
  the 727-pair high-ANI test split.

## Leakage control

The high-ANI split excluded, in both directions, genomes in the v5 training
set and the 43,334-pair held-out set (`scripts/gtdb50k/sample_high_ani_pairs.py`).
Genome-level train/test assignment further dropped mixed pairs. Therefore the
43,334-pair benchmark is independent of v6 training. The additional 727
high-ANI test pairs are also independent of both the v5 and high-ANI training
splits.

## Numerical attribution

The previously generated reports labelled the held-out `ani_cal` column as
"v5"; those labels were stale. The numerical results are v6 results and are now
labelled as v6. No numerical benchmark value has been changed by this
correction.

## Threshold provenance

The 98% hybrid-rule threshold was selected by MAE on the unified benchmark.
Accordingly, unified hybrid rows are operational-threshold sensitivity analyses;
the confirmatory hold-out claims are restricted to the raw gated and
benchmark-v6-calibrated estimates.
