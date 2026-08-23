# Supplementary Information — Syn2b-ANI: Strain-level ANI estimation and structural comparison via fixed restriction-site anchors

This document collects the supplementary figure and tables referenced from
the main text. All underlying data and analysis scripts are released at
https://github.com/HuangShiLab/Syn2bANI-paper.

## Supplementary Figure S1 — Synteny benchmark: exact-truth inversion ladder

(`figures/report/fig_synteny_ladder.png`)

E. coli MG1655 evolved to ANI 95.00/98.00 (counted substitutions) with 0–32
non-overlapping inversions (100–400 kb; 2 true breakpoints each).
(a) syn2bani breakpoint_count vs truth: exact on all 14 rungs (0–64).
(b) synteny_score vs inversion count: monotone decline 1.0000 → 0.9835.
(c) ANI estimate vs inversion count: invariant (94.95–95.00 at true 95.00;
98.00–98.01 at true 98.00). Wall time 50–70 ms per pair. Specificity on 695
real MAG pairs (dnadiff-derived truth, median 0 breaks): median
synteny_score 0.9997, no spurious structural flags. Data and report:
`results/synteny_bench/`.

## Supplementary Figure S2 — Held-out GTDB-R207 benchmark

(`figures/report/fig_gtdb50k_heldout.png`)

43,334 same-genus pairs sampled from GTDB R207 representatives (24,831
genomes, 74 phyla), with every genome of the 2,520-pair calibration training
matrix excluded in both directions; dnadiff/ANIm truth for all pairs. The
gate, consistency flag, and calibration v5 were frozen before evaluation.
Panels show the scored subset (n = 39,903; 3,431 BELOW_DETECTION pairs
unscored by design). (a) Calibrated Syn2bANI v5 vs ANIm (MAE 0.619, bias
-0.12, r = 0.962). (b) skani vs ANIm on the same pairs (MAE 0.957, bias
-0.80, r = 0.973). (c) Per-band MAE. (d) Signed-error distributions.
Per-pair truth, tool output, per-band and per-phylum metrics:
`results/gtdb50k/` (report `GTDB50K_HELDOUT_REPORT.md`, metrics
`gtdb50k_metrics.tsv`); pipeline `scripts/gtdb50k/`.

## Supplementary Table S1 — MAG accuracy by CheckM2 quality tier

695 CAMI2 bins vs dnadiff ANIm truth (raw gated syn2bani; skani; FastANI).
Full per-group metrics: `results/mag_validation/MAG_METRICS.tsv`.

| Tier | n | syn2bani MAE | skani MAE | FastANI MAE | syn2bani within 0.5 |
|---|---|---|---|---|---|
| HQ | 200 | 0.041 | 0.024 | 0.030 | 100% |
| MQ | 210 | 0.105 | 0.067 | 0.131 | 99.0% |
| LQ | 285 | 0.291 | 0.160 | 0.426 | 86.7% |

## Supplementary Table S2 — MAG accuracy by contamination class

Class from CAMI2 ground-truth contig assignment (clean 305 / strain-mixed
253 / cross-species 137).

| Class | syn2bani MAE | syn2bani bias | skani MAE | FastANI MAE | FastANI bias |
|---|---|---|---|---|---|
| clean | 0.084 | +0.083 | 0.030 | 0.055 | +0.044 |
| strain-mixed | 0.170 | +0.159 | 0.111 | 0.165 | −0.081 |
| cross-species | 0.326 | +0.314 | 0.196 | 0.656 | −0.605 |

## Supplementary Table S3 — MAG accuracy by dnadiff alignment-fraction tier

| AF tier | n | syn2bani MAE | syn2bani within 0.5 |
|---|---|---|---|
| strict (≥ 60%) | 638 | 0.129 | 96.7% |
| low-AF (30–60%) | 40 | 0.461 | 72.5% |
| verylow-AF (< 30%) | 17 | 0.744 | 52.9% |

## Supplementary Note 1 — Calibration is input-regime-specific (MAG test)

The 695 MAG anchor pairs were re-run with `ani --calibrate` (deployed v5
ridge model trained on complete GTDB genomes). Calibrated MAE is 1.261
(bias −1.243) versus 0.163 for the raw gated estimate; the degradation is
worst on low-quality (2.307) and cross-species (2.091) bins, and 43
near-clonal pairs are pushed above 100. The model learned the +2.02 raw
bias of complete-genome pairs, while raw bias on MAGs is only +0.16, so the
correction overshoots. This is the direct evidence for the design rule that
calibrated output is a separate, labelled column and never a silent
replacement; on draft/MAG inputs the raw gated estimate is recommended.
Data: `results/mag_validation/collect/ani_calibrated.tsv`.

## Supplementary Note 2 — Calibration learning curve

Band-holdout CV MAE of the v5 ridge model as a function of training-set
size (subsampled from the 2,520-pair ANIm training matrix):

| Training pairs | CV MAE |
|---|---|
| 492 | 0.952 |
| 984 | 0.902 |
| 1,476 | 0.887 |
| 1,969 | 0.893 |

The linear model plateaus at ~1,500 pairs; expanding the training set
further is not expected to improve it (main text, Methods 4.6).

## Supplementary Note 3 — Pipeline failure modes and fixes (MAG benchmark)

The CAMI2 MAG pipeline required six engineering fixes that bear on
reproducibility of metagenome benchmarks generally: (i) dnadiff/nucmer
(MUMmer 3.23) rejects FASTA files containing any embedded whitespace, and
the CAMI2 source genomes carry spaces in both headers and sequence lines —
inputs were sanitized (headers truncated at first space, sequence
whitespace stripped) with cached copies; (ii) the FastANI CLI flag for the
query is `-q`, not `--q`; (iii) conda environments on shared clusters can
be shadowed by user-site packages (numpy ABI mismatch in CheckM2), fixed
with PYTHONNOUSERSITE=1. Full incident log:
`results/mag_validation/RUNBOOK.md`.
