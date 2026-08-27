# Supplementary Information — Syn2b-ANI: Strain-level ANI estimation and structural comparison via fixed restriction-site anchors

This document collects the supplementary figures and tables referenced from
the main text. All underlying data and analysis scripts are released at
https://github.com/HuangShiLab/Syn2bANI-paper.

## Supplementary Figure S1 — Synteny benchmark: exact-truth inversion ladder

(`figures/report/fig_synteny_ladder.png`)

E. coli MG1655 evolved to ANI 95.00/98.00 (counted substitutions) with 0–32
non-overlapping inversions (100–400 kb; 2 true breakpoints each).
(a) syn2bani breakpoint_count vs truth: exact on all 14 rungs (0–64).
(b) anchor_adjacency vs inversion count: monotone decline 1.0000 → 0.9835.
(c) ANI estimate vs inversion count: invariant (94.95–95.00 at true 95.00;
98.00–98.01 at true 98.00). Wall time 50–70 ms per pair. Specificity on 695
real MAG pairs (dnadiff-derived truth, median 0 breaks): median
anchor_adjacency 0.9997, no spurious structural flags. Data and report:
`results/synteny_bench/`.

## Supplementary Figure S2 — Held-out GTDB-R207 benchmark

(`figures/report/fig_gtdb50k_heldout.png`)

43,334 same-genus pairs sampled from GTDB R207 representatives (24,831
genomes, 74 phyla), with every genome of the 2,520-pair calibration training
matrix excluded in both directions; dnadiff/ANIm truth for all pairs. The
gate, consistency flag, and calibration v5 were frozen before evaluation.
Panels (a–c) Calibrated Syn2bANI v5, skani, and FastANI vs ANIm (MAE 0.619, 0.958, 0.977 respectively). (d) Per-band MAE. (e) Signed-error distributions: skani and FastANI carry a large negative bias at divergent bands; the calibrated estimator is centered (bias −0.12). Per-pair truth, tool output, per-band and per-phylum metrics:
`results/gtdb50k/` (report `GTDB50K_HELDOUT_REPORT.md`, metrics
`gtdb50k_metrics.tsv`); pipeline `scripts/gtdb50k/`.

## Supplementary Figure S3 — Unified GTDB-R207 80–100% benchmark

(`figures/gtdb_r207_unified_benchmark.png`)

43,334 held-out same-genus pairs (80–95% plus the 95–100% stratum split into
95–97 and 97–100 by truth) combined with 727 high-ANI test pairs sampled from
non-representative GTDB-R207 genomes. (a) Estimated vs ANIm truth for Syn2bANI
raw gated, Syn2bANI v6-calibrated, Syn2bANI hybrid (raw ≥98% / calibrated
<98%), skani, and FastANI. (b) Per-band MAE. (c) Signed-error distributions.
(d) High-ANI zoom (95–100%). On the full 80–100% range the Syn2bANI hybrid
estimator attains the lowest overall MAE (**0.615**), while recovering the
accurate raw gated estimate in the 97–100% sub-band (MAE 0.23).

## Supplementary Figure S4 — Genome quality does not drive ANI error on GTDB-R207

(`figures/report/gtdb_quality_vs_mae_combined.png`)

39,903 held-out same-genus pairs with GTDB R207 metadata (completeness,
contamination, contig count, mean contig length) and dnadiff/ANIm truth.
Calibrated Syn2bANI v5 overall MAE = 0.619, Pearson r = 0.962. (a) MAE by
minimum pair completeness (n = 3,972 below 70%, MAE 0.82; n = 22,259 at
95–100%, MAE 0.59). (b) MAE by maximum pair contamination. (c) MAE by maximum
contig count. (d) MAE by minimum mean contig length (<5 kb: MAE 0.81; ≥10 kb:
MAE 0.59–0.66). Error bars are 95% confidence intervals by bootstrap. Only the
most fragmented genomes (<5 kb mean contig length) show a modest accuracy
degradation; typical MAG-level completeness/contamination variation has little
effect on ANI accuracy.

## Supplementary Figure S5 — Exact-truth indel ladder and indel sweep

(`figures/report/fig_s5_simulation_indel.png`)

(a) ANI ladder evolved from *E. coli* K-12 MG1655 at 85.0–99.9% true ANI with
a 400 kb inversion plus ~1 deletion per 100 kb; Syn2bANI (4-enzyme panel) MAE
= 0.073 ANI points. (b) Indel sweep at fixed true ANI 95.000 with 0–4 deletions
per 100 kb; error stays between +0.032 and +0.184 (MAE 0.081). Gamma and
uniform estimates are identical on uniform-rate data.

## Supplementary Figure S6 — GC coverage ladder

(`figures/report/fig_s6_simulation_gc.png`)

Substitution ladders evolved on five real template genomes spanning GC
27.2–72.1%. Syn2bANI error is 0.074–0.356 ANI points and is *not* monotone in
GC; the worst genome (*B. longum*, GC 60.1%) sits mid-range. The GC sweep
source genomes are not present in the repository, so this panel reproduces the
historical 4-enzyme sweep reported in `ALGORITHM_MLE.md` §4.8.

## Supplementary Figure S7 — Simulated fragmentation

(`figures/report/fig_s7_simulation_fragment.png`)

A 95% ANI *E. coli* query was fragmented into 20–201 contigs (~50%
reverse-complemented, order shuffled, no sequence lost). Syn2bANI error is
0.018–0.179 (MAE 0.070), confirming that random fragmentation per se does not
bias the chain-restricted MLE. Real assemblies below ~20 kb N50 behave
differently because they lose sequence at repeat boundaries; see Results 2.3
and Fig. 10c.

## Supplementary Figure S8 — Accessory-content confound

(`figures/report/fig_s8_simulation_accessory.png`)

Core ANI fixed at 95.000 with 0–50% of the genome replaced by
composition-preserving, homology-destroying shuffled blocks. (a) ANI estimate
remains flat across accessory fractions (+0.044 to +0.250, MAE 0.112) because
only tags inside chains contribute to the likelihood. (b) `af_query` tracks
the true shared fraction (1 − F) to within 0.004, demonstrating the decoupling
of divergence and shared-content estimates.

## Supplementary Figure S9 — Mosaic/rate-heterogeneity family

(`figures/report/fig_s9_simulation_mosaic.png`)

Per-block divergence rates sampled from Gamma(α, α) with α = 0.5/1/2 (mean
ANI 90–98) and deliberately misspecified bimodal 50/70%-core cases. The gamma
estimator reduces error relative to the uniform estimator (e.g., gamma MAE
1.35 vs uniform 2.75 on gamma regimes), but residual bias remains on bimodal
misspecification (gamma MAE 2.43 vs uniform 3.42), motivating the real-genome
ridge calibration of Results 2.5.

## Supplementary Figure S10 — Breakpoint counts scale independently of ANI in high-ANI isolate collections

(`figures/syntracker_validation/syntracker_breakpoints_vs_ani.png`)

Syn2bANI ANI versus breakpoint count for (a) *E. coli* hypermutator (23
isolates, 253 pairs) and (b) *H. pylori* (77 isolates, 2,926 pairs). In both
collections ANI is pinned near 100% while breakpoint counts vary over more
than an order of magnitude, illustrating that rearrangement burden is not
captured by the divergence scalar. Annotated points are the top discordant
cases by ANI–synteny rank difference.

## Supplementary Figure S11 — Near-clonal *Neisseria gonorrhoeae* and *Streptomyces rimosus* isolates also show ANI–anchor-adjacency decoupling

(`figures/syntracker_validation/syntracker_supp_ngonorrhoeae_srimosus.png`)

Left column: Syn2bANI ANI vs anchor adjacency; right column: skani ANI vs the
same anchor adjacency. (a,b) 12 *N. gonorrhoeae* isolates (66 pairs); pairs with
ANI >99.9% have anchor adjacencies 0.85–0.89 and up to 376 breakpoints. (c,d) 20
*S. rimosus* isolates (190 pairs); near-clonal pairs (ANI 99.99–100%) have
anchor adjacencies 0.89–0.96 and an average of ~1,980 breakpoints. The most
discordant cases (NG_05 vs NG_10 and SR_06 vs SR_10) were re-analyzed with
`syn2bani struct`; counts are in Table 4 and
`results/syntracker_validation/struct_top_cases_summary.tsv`.

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
near-clonal pairs are pushed above 100. The model learned the +2.02 raw bias
of complete-genome pairs, while raw bias on MAGs is only +0.16, so the
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
