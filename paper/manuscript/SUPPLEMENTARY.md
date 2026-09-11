# Supplementary Information — Syn2b-ANI: Strain-level ANI estimation and structural comparison via fixed restriction-site anchors

This document collects the supplementary figures and tables referenced from
the main text. All underlying data and analysis scripts are released at
https://github.com/HuangShiLab/Syn2bANI-paper.

## Supplementary Figure S1 — Synteny benchmark: exact-truth inversion ladder

(`paper/figures/supplementary/fig_s1_synteny_ladder.png`)

E. coli MG1655 evolved to ANI 95.00/98.00 (counted substitutions) with 0–32
non-overlapping inversions (100–400 kb; 2 true breakpoints each).
(a) syn2bani breakpoint_count vs truth: exact on all 14 rungs (0–64).
(b) anchor_adjacency vs inversion count: monotone decline 1.0000 → 0.9835.
(c) ANI estimate vs inversion count: invariant (94.95–95.00 at true 95.00;
98.00–98.01 at true 98.00). Wall time 50–70 ms per pair. Specificity on 695
real MAG pairs (dnadiff-derived truth, median 0 breaks): median
anchor_adjacency 0.9997, no spurious structural flags. Data and report:
`results/synteny_bench/`.

## Supplementary Figure S2 — Held-out GTDB-R207 benchmark (v5 detail)

(`paper/figures/supplementary/fig_s2_gtdb50k_heldout.png`)

43,334 same-genus pairs sampled from GTDB R207 representatives (24,831
genomes, 74 phyla), with every genome of the 2,520-pair calibration training
matrix excluded in both directions; dnadiff/ANIm truth for all pairs. The
gate, consistency flag, and calibration v5 were frozen before evaluation.
Panels (a–c) Calibrated Syn2bANI v5, skani, and FastANI vs ANIm (MAE 0.619, 0.958, 0.977 respectively). (d) Per-band MAE. (e) Signed-error distributions: skani and FastANI carry a large negative bias at divergent bands; the calibrated estimator is centered (bias −0.12). Per-pair truth, tool output, per-band and per-phylum metrics:
`results/gtdb50k/` (report `GTDB50K_HELDOUT_REPORT.md`, metrics
`gtdb50k_metrics.tsv`); pipeline `scripts/gtdb50k/`.

## Supplementary Figure S3 — Per-band ANI comparison

(`paper/figures/supplementary/fig_s3_anim_by_band.png`)

Band-holdout cross-validation on the 2,074-pair GTDB-R207 ANIm benchmark:
MAE by ANI band for raw gamma, calibrated v5, skani, and FastANI. This is the
per-band detail behind the unified benchmark of Fig. 3.

## Supplementary Figure S4 — Genome quality does not drive ANI error on GTDB-R207

(`paper/figures/supplementary/fig_s4_gtdb_quality_vs_mae.png`)

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

## Supplementary Figure S5 — GC coverage ladder

(`paper/figures/supplementary/fig_s5_simulation_gc.png`)

Substitution ladders evolved on five real template genomes spanning GC
27.2–72.1%. Syn2bANI error is 0.074–0.356 ANI points and is *not* monotone in
GC; the worst genome (*B. longum*, GC 60.1%) sits mid-range.

## Supplementary Figure S6 — Mosaic/rate-heterogeneity family

(`paper/figures/supplementary/fig_s6_simulation_mosaic.png`)

Per-block divergence rates sampled from Gamma(α, α) with α = 0.5/1/2 (mean
ANI 90–98) and deliberately misspecified bimodal 50/70%-core cases. The gamma
estimator reduces error relative to the uniform estimator (e.g., gamma MAE
1.35 vs uniform 2.75 on gamma regimes), but residual bias remains on bimodal
misspecification (gamma MAE 2.43 vs uniform 3.42), motivating the real-genome
ridge calibration of Results 2.5.

## Supplementary Figure S7 — Accuracy on 695 CAMI2 MAGs

(`paper/figures/supplementary/fig_s7_mag_validation.png`)

695 CAMI2 bins vs dnadiff ANIm truth (raw gated syn2bani; skani; FastANI).
(a) Syn2bANI raw gated estimate vs dnadiff ANIm truth, colored by contamination
class. (b) Absolute-error distributions by tool. (c) Syn2bANI error by CheckM2
quality tier.

## Supplementary Figure S8 — High-ANI isolate collections show extensive rearrangements

(`paper/figures/supplementary/fig_s8_syntracker_breakpoints.png`)

skani ANI versus Syn2bANI `breakpoint_count` for four published near-clonal
isolate collections: (a) *E. coli* hypermutator (253 pairs), (b) *H. pylori*
(2,926 pairs), (c) *N. gonorrhoeae* (66 pairs), and (d) *S. rimosus* (190
pairs). In all four collections ANI is pinned near 100% while breakpoint counts
vary over orders of magnitude.

## Supplementary Figure S9 — High-ANI *E. coli* O157:H7 pairs carry hundreds of breakpoints

(`paper/figures/supplementary/fig_s9_ecoli_o157_breakpoints.png`)

Seventy-four genomes from Fitzgerald et al. (2021), 2,701 non-self pairs.
(a) ANI vs breakpoint count colored by lineage (I/II, II, Ia, Ic). (b) Same
data colored by host category (bovine, human, other/unknown). All pairwise ANIs
exceed 99.886% yet breakpoints range from 171 to >1,100.

## Supplementary Figure S10 — High-ANI FDA-ARGOS *Staphylococcus aureus* pairs show wide breakpoint variation

(`paper/figures/supplementary/fig_s10_saureus_breakpoints.png`)

One hundred and twenty-two genomes, 7,381 refined pairs. (a) ANI vs breakpoint
count colored by country. (b) ANI vs breakpoint count colored by isolation
source. Pairs at 100% ANI still carry >150 breakpoints.

## Supplementary Figure S11 — cagPAI extended state by country and phylogenetic population

(`paper/figures/supplementary/fig_s11_cagpai_country_population.png`)

Stacked-bar distribution of extended cagPAI states (after circular-origin
filtering) across countries of isolation and phylogenetic population for the
528 *H. pylori* cohort genomes.

## Supplementary Figure S12 — Circular-origin artifact filtering and sensitivity in *H. pylori* cagPAI

(`paper/figures/supplementary/fig_s12_circular_origin_filtering.png`)

(a, b) Before and after counts of the four extended cagPAI states. Filtering
reclassifies 133 complete-marker genomes from `complete_rearranged` to
`complete_collinear` because their only cagPAI-overlapping SV is a
genome-spanning call consistent with a different arbitrary start coordinate on
the circular *H. pylori* chromosome. (c) Number of `complete_rearranged` genomes
and number of excluded calls as a function of the artifact-span threshold.
(d) Crude association p-value for cagPAI rearrangement (GC vs NAG) across the
same threshold range; the association remains non-significant throughout.
Source data: `case_studies/h_pylori_cagpai/results/cagpai_artifact_threshold_sensitivity.tsv`.

## Supplementary Figure S13 — Lineage-stratified association forest plot for *H. pylori* metadata tests

(`case_studies/h_pylori_cagpai/results/metadata_assoc/fig_forest_stratified.png`)

Crude (circle) and FastBAPS-lineage-stratified Mantel–Haenszel (triangle) odds
ratios with 95% CIs for the headline association tests: cagPAI presence vs cagA
detection and vacA s1/i1/m1 (top, strong and lineage-robust), and SV burden vs
macrolide resistance (bottom, crude OR ≈ 2.8 collapsing to ≈ 1.3 after
stratification). The macrolide rows straddle OR = 1, visually demonstrating the
lineage confounding described in the main text.
Source data: `case_studies/h_pylori_cagpai/results/metadata_assoc/association_tests.tsv`.



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

## Supplementary Table S4 — Crude and lineage-stratified cagPAI–disease associations

Five hundred and twenty-eight *H. pylori* isolates from Song et al. (2026) were
classified into extended cagPAI states after circular-origin filtering (empty
85, partial 11, complete_collinear 145, complete_rearranged 287). The raw
Pearson χ² test for disease stage is significant (χ² = 24.60, df = 9, p = 0.0034),
but the association is confounded by FastBAPS lineage. The table below reports
crude odds ratios (Woolf 95% CI), χ² p-values, Cochran–Mantel–Haenszel (CMH)
p-values and Mantel–Haenszel odds ratios, and Breslow–Day tests of homogeneity
across FastBAPS strata.

| Contrast | n | Case | Control | Outcome | crude OR (95% CI) | crude p | CMH χ² | CMH p | OR_MH | Breslow–Day p |
|---|---:|---|---|---|---|---:|---:|---:|---:|---:|---:|
| cagPAI presence | 285 | GC | NAG | complete vs empty/partial | 2.16 (1.13–4.13) | 0.018 | 2.406 | 0.121 | 1.70 | 0.005 |
| cagPAI rearrangement | 230 | GC | NAG | complete_rearranged vs complete_collinear | 1.11 (0.64–1.93) | 0.700 | 0.315 | 0.574 | 1.17 | 0.153 |
| cagPAI presence (advanced vs early) | 528 | GC/IM | AG/NAG | complete vs empty/partial | 1.14 (0.73–1.78) | 0.576 | 1.864 | 0.172 | 1.39 | 0.021 |
| cagPAI rearrangement (advanced vs early) | 432 | GC/IM | AG/NAG | complete_rearranged vs complete_collinear | 1.02 (0.68–1.52) | 0.928 | 0.004 | 0.951 | 1.01 | 0.780 |

Only the first contrast (presence, GC vs NAG) showed a nominally significant
crude association. After FastBAPS stratification the CMH p-value is 0.12 and the
Mantel–Haenszel OR is 1.70. The Breslow–Day test is significant (p = 0.005),
indicating that the OR is heterogeneous across lineages. The other three
contrasts were not significant even in the crude analysis. Discordant-cell
counts are small: only 15 cagPAI-negative gastric-cancer genomes exist across
the whole cohort, distributed as 4 / 9 / 1 / 1 / 0 across FastBAPS L2–L6. Thus
the presence association is attenuated but not abolished by lineage
stratification, and the current design is under-powered to resolve it.

Rearrangement-state associations were also evaluated across artifact-span
thresholds from 0.20 to 0.50; p-values remained non-significant throughout
(`case_studies/h_pylori_cagpai/results/cagpai_artifact_threshold_sensitivity.tsv`).
Full per-stratum counts and CMH diagnostics:
`case_studies/h_pylori_cagpai/results/cagpai_association_crude_and_stratified.tsv`.

## Supplementary Table S5 — cagPAI cross-validation against independent cagA/vacA typing and SV–resistance associations

Five hundred and twenty-eight isolates were merged with collaborator-generated
metadata (cagA detection with EPIYA genotyping and genome-level rescue, vacA
s/i/m typing, CARD resistance hits; merge losses: none). Crude odds ratios
(Woolf 95% CI), Fisher p, Cochran–Mantel–Haenszel (CMH) stratified OR and
p, Breslow–Day homogeneity (B–D p), and Benjamini–Hochberg q across the
executed test family are shown; q ≈ 0 denotes underflow.

| Test | n | crude OR (95% CI) | crude p | stratified OR (95% CI) | CMH p | B–D p | q |
|---|---:|---|---|---|---|---|---:|---:|
| cagPAI complete vs cagA detected | 524 | 782.9 (179.2–3420.8) | 7.2 × 10⁻⁶⁸ | 3986 (263.3–60342.6) | < 10⁻¹⁶ | 0.956 | ≈ 0 |
| cagPAI complete vs cagA EPIYA genotype (4-level) | 450 | — | 0.515 (χ²) | — | — | — | 0.57 |
| cagPAI complete vs vacA s1 | 462 | 556.1 (74.0–4180.4) | 1.4 × 10⁻³⁹ | 510.9 (74.1–3520.8) | < 10⁻¹⁶ | 1 × 10⁻⁶ | ≈ 0 |
| cagPAI complete vs vacA i1 | 477 | 38.0 (19.6–73.7) | 4.6 × 10⁻³² | 27.0 (13.6–53.8) | < 10⁻¹⁶ | 0.77 | ≈ 0 |
| cagPAI complete vs vacA m1 | 519 | 11.8 (6.6–21.4) | 1.4 × 10⁻²¹ | 11.4 (6.4–20.4) | < 10⁻¹⁶ | 0.018 | ≈ 0 |
| SV burden (median split) vs macrolide resistance | 528 | 2.81 (1.71–4.62) | 3.0 × 10⁻⁵ | 1.33 (0.73–2.40) | 0.35 | 0.90 | 8.5 × 10⁻⁵ |
| cagPAI complete vs GC (vs NAG) | 285 | 2.16 (1.13–4.13) | 0.022 | 1.70 (0.87–3.32) | 0.121 | 0.005 | 0.044 |

Contingency for the headline test: cagPAI-complete → 429/431 cagA-positive;
partial → 10/11; empty → 10/82. Fluoroquinolone, tetracycline, and
nitroimidazole classes were excluded pre-test (intrinsic hp1181 efflux pump at
99.8% prevalence); beta-lactam was excluded (0.6% prevalence). The 26 genomes
discordant between marker state and cagA typing were reviewed at assembly level:
17 showed loss of the island backbone with cagA sequence retained at the
evacuated site (8 near full-length, 9 fragmentary), 2 carried a precise cagA
excision on an otherwise intact island (both FastBAPS L3, Latin America), 5
were threshold-level partial calls in high-divergence backgrounds, and 2 were
clean fully assembled island deletions; none was attributable to assembly
incompleteness.
Source data: `case_studies/h_pylori_cagpai/results/metadata_assoc/association_tests.tsv`;
assembly-level review: `case_studies/h_pylori_cagpai/results/metadata_assoc/discordant_review/`.

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

## Supplementary Note 4 — Chain-restricted stratified MLE, gating, and flags

For a query tag inside a chained region, with tag length `k`, recognition-site
length `s`, mutable body `b = k − s`, and mismatch budget `tol`, the outcome
probabilities under per-base identity `a` are:

- **Found with `m` body mismatches** (for `m ≤ tol`):
  `P_m(a) = C(b, m) · (1 − a)^m · a^(k−m)`.
- **Not found** (miss): `P_miss(a) = 1 − Σ_{m=0}^{tol} P_m(a)`.

Mismatches are observable only in the mutable body — a mutation inside the
recognition site deletes the tag entirely — so the site contributes a factor
`a^s` to `P_found` but never appears in the histogram. This also sets a hard
retention ceiling (for example, ~0.72 at 95% ANI and ~0.50 at 90% with the
default budget). The per-enzyme strata differ in `k` and `b`; the total
log-likelihood is the sum over strata and is maximized over the scalar `a` by
golden-section search.

Two partial estimators are computed as diagnostics:

- `ani_from_loss` — uses `P_miss(a)` only (the miss-rate channel).
- `ani_from_hist` — uses the mismatch histogram renormalized by `P_found(a)`.

When the data are well behaved these two estimates agree; large divergence
signals either model misspecification or sparse data.

**Heterogeneous (gamma) model.** To account for rate heterogeneity, each
chained region receives a per-region rate multiplier `r ~ Gamma(α, α)` (mean 1,
shape α). Because rate variation acts at kilobase scale while a tag is ~30 bp,
all sites within one tag share the same `r`. Integrating `r` out converts the
mismatch count distribution into a negative binomial with two parameters: mean
divergence `d` and shape `α`. These are identifiable from the three histogram
degrees of freedom plus the miss count. Two guardrails are applied: the second
parameter is used only when the likelihood-ratio statistic exceeds 3.841, and
`α` is clamped to `[0.1, 200]`.

At the identifiability boundary with few surviving tags, the gamma shape and
mean couple and can overshoot by 4–10 ANI points in the mid-ANI range. The
shipped point estimate is therefore **gated**: report the homogeneous (uniform)
fit when `|ani_from_loss − ani_from_hist| > 5` ANI points, otherwise report the
gamma fit. The threshold is an effect size, not a significance level; threshold
sweeps on the GTDB-ANIm training matrix show a flat optimum between 4.5 and 6
points. The gate fires on 5.5% of GTDB pairs, on all 15 mid-ANI pairs (where it
reduces gamma MAE from 4.48 to 1.42), and never on uniform-rate simulations,
mosaic simulations, or oral/gut same-species pairs, preserving gamma's advantage
exactly where it is real.

**Flags.** The `flag` column is hierarchical:

- `BELOW_DETECTION` — expected retention < 0.20 (takes precedence).
- `INCONSISTENT` — the gate fell back to the uniform model, or the chains carry
  > 0.5 rearrangement breakpoints per anchor (a structural statistic independent
  of the chain-restricted likelihood denominator).
- `ok` — otherwise.

The recalibrated flag never inverts its error ranking on any validation set
(kept vs flagged MAE: GTDB-ANIm 2.415 vs 4.153; oral/gut 0.514 vs 4.269;
mid-ANI 0.343 vs 1.113), at a disclosed cost in near-clonal sensitivity.

## Supplementary Note 5 — The spatial-model negative result (identifiability analysis)

A mechanistic spatial-rate model was evaluated as an alternative to the ridge
calibration layer. The prototype was written in Python, tested on 19 exact-truth
mosaic simulations, and then evaluated on the GTDB-ANIm matrix.

The residual bias decomposes into two terms. The first is an in-chain
distribution-family term: the asymptotic in-chain MAE is 2.25 for gamma and 1.25
for a capped grid NPMLE, so a nonparametric in-chain fit can fix at most
~0.1–0.2 ANI points of the GTDB MAE. The second is a coverage term: the
divergence of the unchained fraction is not identifiable from tag data, because
the out-of-chain anchor residual is approximately zero within multi-match noise
whether the unchained mass is saturated-divergent or accessory. Even an oracle
given the exact identity of the chained sample still errs by MAE 1.36 against
whole-genome truth.

Candidates tested included AF-weighted mixtures, discrete mixtures,
ascertainment-aware tilts, capped NPMLE with and without likelihood-ratio gating,
and model averaging. None beats the gated homogeneous/gamma baseline while
holding all simulation gates. The raw estimator is therefore at its
identifiability floor at 4-enzyme tag density; this negative result is why the
shipped pipeline uses calibration, rather than a more elaborate likelihood, as
the real-genome correction layer.

## Supplementary Note 6 — Ridge calibration protocol details

The v5 calibration model is ridge regression on nine Syn2bANI-internal features:
`ani_gated`, `ani_uniform`, `af_query`, `af_reference`, `std_err`, `retention`,
`n_anchors`, `n_chains`, and `n_tags_in_chains`. Features are median-imputed and
standardized; regularization strength `α = 10` is selected by inner
cross-validation with scikit-learn29.

**Training set.** 2,520 ANIm-truth pairs (80–99.5% ANIm). This combines 2,053
finite pairs from a stratified 2,074-pair GTDB benchmark and 467 targeted
95.0–99.15% ANIm pairs selected by screening all 2.26 M same-genus GTDB-R207
representative pairs with skani and running dnadiff on 650 candidates. The
expanded set spans 607 species and 26 phyla and has zero genome-level overlap
with the original 2,074-pair set.

**Validation.** Band-holdout cross-validation: each of the four ANI bands
(80–85, 85–90, 90–95, 95–99) is predicted by a model trained on the other three
only. Training rows are shuffled with a fixed seed before CV because the inner
KFold of the ridge CV is order-sensitive.

**Model comparisons.** An expanded 18-feature variant and a gradient-boosted
regressor on the same features were evaluated and rejected: nonlinearity buys
nothing at `n ≈ 2,500`, and the expanded feature set degrades out-of-distribution
on near-clonal pairs. Seed spread over 5 seeds is ≤ 0.0015 MAE.

**Learning curve.** Subsampling the 2,520-pair training matrix at 492/984/1,476
and 1,969 pairs gives band-holdout CV MAE of 0.952, 0.902, 0.887, and 0.893,
respectively (Supplementary Note 2); the linear model plateaus at ~1,500 pairs.

**Deployment.** The model ships as a versioned JSON file containing ridge
weights, scaling parameters, and feature order; it is embedded in the binary and
evaluated at runtime in < 1 µs. All feature matrices were regenerated with the
post-rescue binary before training, so there is no version skew between the
shipped estimator and the deployed model. Calibrated output is a separate column
(`--calibrate`); it returns no value on `BELOW_DETECTION` pairs rather than
extrapolating.

**Hybrid estimator threshold.** The raw gated estimate is more accurate than the
calibrated estimate in the near-clonal regime, whereas calibration is essential
at lower ANI. A hybrid rule was therefore evaluated: report the raw gated
estimate when it exceeds a threshold `t`, otherwise report the calibrated
estimate. The threshold was chosen by overall MAE on the unified 80–100% ANIm
benchmark (43,334 held-out GTDB-R207 pairs plus 727 high-ANI test pairs,
n = 40,629 pairs with non-missing raw and calibrated values):

| Threshold `t` (%) | n | MAE | Bias |
|---|---:|---:|---:|
| 95.0 | 40,629 | 0.7065 | +0.0284 |
| 96.0 | 40,629 | 0.6443 | −0.0699 |
| 97.0 | 40,629 | 0.6189 | −0.1089 |
| 97.5 | 40,629 | 0.6161 | −0.1133 |
| 98.0 | 40,629 | 0.6146 | −0.1154 |
| 98.5 | 40,629 | 0.6144 | −0.1169 |
| 99.0 | 40,629 | 0.6140 | −0.1191 |
| 99.5 | 40,629 | 0.6145 | −0.1211 |

MAE decreases sharply from 95% to 97%, then plateaus; 98% was selected as the
operating threshold because it lies on the plateau while keeping the calibrated
estimate active for the broadest range of strain-level pairs. The hybrid rule is
not sensitive to the exact threshold in the 97–99% window.
