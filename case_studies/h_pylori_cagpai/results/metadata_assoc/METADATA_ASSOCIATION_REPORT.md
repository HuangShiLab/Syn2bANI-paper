# H. pylori SV × isolate-metadata association analysis

**Date:** 2026-09-11
**Driver script:** `case_studies/h_pylori_cagpai/scripts/run_metadata_association.py`
**Data:** `/Volumes/MoneyCat/Data/song_2026_hpylori/` (read-only)
**Cohort:** 528 *Helicobacter pylori* genome assemblies (keyed by `GCA_accession`)

This is an **exploratory, hypothesis-generating** analysis. Its two main goals were:

1. **Cross-validate our cagPAI marker panel against independent cagA molecular
   typing** (new collaborator data), plus vacA typing.
2. **Test whether structural-variation (SV) burden is associated with
   antibiotic-resistance genotype** from CARD (new data).

Disease-outcome (`group`) tests are secondary and included only for completeness.

---

## 1. Data inventory

Merged per-genome table: `merged_genome_table.tsv` — **528 genomes × 55 columns**.
**Merge losses: none.** All 528 GCAs have SV calls, cagPAI status, metadata,
cagA and vacA typing; CARD long format covers all 528 genomes (618 hit rows).

| Variable | n (non-missing) | Notes |
|---|---|---|
| cagPAI status (28-marker panel) | 528 | complete 432 (81.8%), empty 85 (16.1%), partial 11 (2.1%) |
| cagPAI extended state | 528 | complete_collinear 145, complete_rearranged 287, empty 85, partial 11 |
| cagA result | 528 | Detected 449, Not_detected 75, Indeterminate 2, Detected_partial 1, Detected_disrupted 1 |
| cagA genotype (EPIYA) | 450 | ABC 177, ABD 171, ABCC 34, AABD 21, AB 9, other/multi 48; 78 missing |
| East-West type | 439 | Western_like 241, EastAsian_like 198; 78 NaN + 11 "atypical" excluded |
| vacA s / i / m | 462 / 477 / 519 | s1 418, s2 44; i1 409, i2 68; m1 311, m2 208 |
| fastbaps lineage | 528 | L2 197, L3 149, L5 113, L4 48, L6 21 (no L1 in cohort) |
| Disease group | 528 | GC 118, AG 133, NAG 167, IM 110 |
| Country | 515 | 13 missing |
| CARD macrolide resistance (23S rRNA A2146G/A2147G) | 85/528 (16.1%) | Perfect/Strict hits only |
| CARD beta-lactam (TEM-116/-229) | 3/528 (0.6%) | very rare |
| CARD fluoroquinolone / tetracycline / nitroimidazole | 527/528 (99.8%) | **intrinsic hp1181 MFS efflux pump — near-fixation, non-informative** |
| SV count vs 26695 | 528 | median 36, IQR 32–40, range 18–50 |
| SV span vs 26695 | 528 | median ~1.10 Mb |
| SV >50 kb count | 528 | median 6, IQR 4–7 (527/528 have ≥1 large SV) |

**Important CARD caveat.** The only discriminatory resistance signal in this
dataset is the macrolide class (23S rRNA mutation, n=85). The
fluoroquinolone/tetracycline/nitroimidazole "resistance" calls are the
chromosomally intrinsic hp1181 MFS efflux pump (529 hits over 528 genomes)
annotated by CARD against those drug classes; prevalence 99.8% makes them
untestable as phenotypic resistance proxies. No gyrA/rdxA mutational FQ or
nitroimidazole calls were detected. Beta-lactam (TEM) hits are too rare
(n=3) to test. **Consequently the only executable SV×resistance comparison is
macrolide (clarithromycin proxy), for which genotype–phenotype concordance is
literature-confirmed per the collaborator.**

## 2. Methods

- **Crude 2×2:** Fisher exact; Woolf OR with Haldane–Anscombe 0.5 correction
  for zero cells.
- **Multi-level tables:** Pearson chi-square (low-count cells flagged).
- **Lineage stratification:** Mantel–Haenszel (CMH) pooled OR + conditional
  independence test across fastbaps strata (statsmodels `StratifiedTable`,
  Robins–Breslow–Greenland CI); Breslow–Day homogeneity computed manually
  (MH pooled OR, df = k−1). For continuous SV outcomes: lineage-adjusted
  logistic regression of above-median burden (fastbaps dummies).
- **FDR:** Benjamini–Hochberg across the full family of 20 executed tests
  (4 resistance classes skipped pre-test for non-informative prevalence are
  not counted).
- Cells with n<5 flagged `low_count`; several headline tables have zero cells
  by construction (e.g. only 2 cagPAI-complete genomes are cagA-negative), so
  CIs are wide even when p is extreme.

## 3. Results

Full machine-readable table: `association_tests.tsv` (24 rows: 20 executed +
4 skipped). q = BH-FDR across executed tests.

### 3.1 cagPAI × cagA — cross-validation of the marker panel (headline)

| Test | n | Crude OR (95% CI) | p | Lineage-stratified OR (95% CI) | CMH p | B–D p | q_FDR |
|---|---|---|---|---|---|---|---|
| cagPAI complete vs cagA detected | 524 | **782.9 (179.2–3420.8)** | 7.2e-68 | **3986 (263.3–60342.6)** | <1e-16† | 0.956 | ≈0 |
| cagPAI status 3-level vs cagA detected | 524 | — | 7.9e-94 (chi²) | — | — | — | 3.1e-93 |
| cagPAI complete vs cagA genotype (ABC/ABD/ABCC/multi) | 450 | — | 0.515 (chi²) | — | — | — | 0.57 |
| cagPAI complete vs East-Asian type | 439 | 0.81 (0.33–2.00) | 0.654 | 2.53 (0.002–3173)‡ | 0.847 | 0.906 | 0.65 |
| cagPAI rearranged vs collinear (among complete) vs cagA | 431 | 0.39 (0.019–8.20) | 0.552 | — (degenerate) | — | — | 0.58 |
| cagPAI rearranged vs collinear vs East-Asian | 419 | 0.85 (0.57–1.27) | 0.469 | 0.21 (0.016–2.76) | 0.246 | 0.437 | 0.38 |

† p underflows to 0.0 in floating point; treat as p < 1e-16.
‡ uninformative: near-empty strata make the stratified CI explode; direction
  and crude estimate both null.

Contingency (cagPAI status × cagA result): complete → 429/431 cagA-detected;
partial → 10/11 detected; empty → 10/82 detected. The two cagPAI-complete but
cagA-negative genomes and the 10 empty-but-cagA-detected genomes are the
discordant sets worth manual review (possible cagA sequence divergence below
the typing threshold, or marker-panel false positives/negatives).

**Interpretation.** The 28-marker cagPAI panel reproduces the independent cagA
PCR/in-silico typing with near-perfect concordance, and the association is
**homogeneous across lineages** (Breslow-Day p=0.96) — it is not driven by any
single lineage. cagPAI presence does **not** predict cagA EPIYA genotype
structure (p=0.52), as expected: the panel assays island presence, not the
toxin allele architecture. No East-West association is detectable in this cohort.

### 3.2 cagPAI × vacA

| Test | n | Crude OR (95% CI) | p | Stratified OR (95% CI) | CMH p | B–D p | q_FDR |
|---|---|---|---|---|---|---|---|
| cagPAI complete vs vacA s1 | 462 | **556.1 (74.0–4180.4)** | 1.4e-39 | **510.9 (74.1–3520.8)** | <1e-16† | 1e-6 (heterog.) | ≈0 |
| cagPAI complete vs vacA i1 | 477 | **38.0 (19.6–73.7)** | 4.6e-32 | **27.0 (13.6–53.8)** | <1e-16† | 0.77 | ≈0 |
| cagPAI complete vs vacA m1 | 519 | **11.8 (6.6–21.4)** | 1.4e-21 | **11.4 (6.4–20.4)** | <1e-16† | 0.018 (heterog.) | ≈0 |

Only 1 of 389 cagPAI-complete genomes is vacA s2, vs 30/73 s1 among
cagPAI-incomplete. These strong, lineage-robust links are consistent with the
known tight co-inheritance of cagPAI/cagA with vacA s1/i1 (the toxigenic
strain genotype); the Breslow-Day signals reflect magnitude (not direction)
differences among lineages. These are sanity confirmations, not novel findings.

### 3.3 SV burden × resistance (the genuinely new test)

| Test | n | Crude | p | q_FDR | Lineage-stratified | Strat. p |
|---|---|---|---|---|---|---|
| SV burden median-split vs macrolide R | 528 | OR 2.81 (1.71–4.62) | **3.0e-5** | **8.5e-5** | OR 1.33 (0.73–2.40), CMH | 0.35 |
| SV count continuous vs macrolide R | 528 | median 38 vs 36 (+2) | 6.5e-5 (MWU) | 1.5e-4 | OR 1.34 (0.73–2.45), logistic | 0.34 |
| SV span continuous vs macrolide R | 528 | median 1.165 vs 1.049 Mb | 0.060 | 0.10 | OR 1.18 (0.72–1.93), logistic | 0.51 |
| large-SV (>50 kb) count median-split vs macrolide R | 528 | OR 1.22 (0.75–2.00) | 0.44 | 0.47 | OR 1.25 (0.75–2.09), CMH | 0.40 |
| fluoroquinolone / tetracycline / nitroimidazole | — | **skipped (prevalence 99.8%, intrinsic hp1181)** | — | — | — | — |
| beta-lactam | — | **skipped (prevalence 0.6%)** | — | — | — | — |

**Interpretation (candid).** The crude SV-burden–macrolide association is
significant and survives FDR, **but it does not survive lineage
stratification**: CMH OR collapses from 2.8 to 1.33 (95% CI 0.73–2.40,
p=0.35), with no inter-lineage heterogeneity (B–D p=0.90). The pattern is
classic confounding: fastbaps L2 (hspEAsia, mostly Chinese isolates) carries
both the highest macrolide-resistance prevalence (55/197 = 28%) and the
highest SV burden (166/197 above the cohort median), while L3/L5 carry little
of either. The same caution that applied to the earlier cagPAI–disease-stage
analysis applies here: **unstratified SV-resistance correlations in this
cohort reflect lineage structure, not biology of resistance.** There is no
evidence that high SV burden per se marks clarithromycin-resistant isolates
once lineage is accounted for.

### 3.4 Disease outcome (secondary)

| Test | n | Crude OR (95% CI) | p | Stratified | Strat. p | q_FDR |
|---|---|---|---|---|---|---|
| cagPAI complete vs GC (vs NAG) | 285 | 2.16 (1.13–4.13) | 0.022 | 1.70 (0.87–3.32) | 0.121 | 0.044 |
| macrolide R vs GC (vs NAG) | 285 | 1.05 (0.50–2.19) | 1.0 | 0.70 (0.33–1.50) | 0.351 | 0.47 |

The crude cagPAI–GC signal weakens after lineage adjustment (p=0.12), fully
consistent with the previous analysis — we report it as such and do not claim
a disease association.

## 4. What failed or is limited

- **No executable fluoroquinolone (levofloxacin) SV test**: CARD FQ signal is
  the intrinsic hp1181 efflux pump at 99.8% prevalence. A gyrA-based test
  would require calling gyrA QRDR mutations directly (not in the collaborator
  data). The headline "SV×levofloxacin" question remains **unanswered**, not
  negative.
- **Beta-lactam untestable** (n=3 TEM hits; for H. pylori these may also be
  assembly contamination — worth a glance before any claim).
- **East-West stratified estimates are unstable** (empty/sparse strata);
  treat the null there as "no detectable signal in this cohort", not proof of
  no association.
- **12 discordant genomes** between cagPAI panel and cagA typing (2
  complete/cagA-negative, 10 empty/cagA-detected) — candidates for manual
  curation; they do not threaten the headline concordance.
- Low-count cells in several tables (flagged in `association_tests.tsv`);
  CIs for the stratified cagA/vacA ORs are extremely wide even though
  direction is unanimous across strata — the point estimate is best read
  together with the contingency table.

## 5. Figures (300 dpi PNG)

- `fig_cagpai_by_cagA_eastwest.png` — cagPAI status stacked bars by East-West
  type and cagA EPIYA genotype group.
- `fig_sv_burden_by_resistance.png` — SV-count distributions by macrolide
  (23S mutation) and by the (non-informative) FQ-class indicator, which
  visually demonstrates the hp1181 near-fixation problem.
- `fig_forest_stratified.png` — lineage-stratified ORs for all headline
  tests; note macrolide→SV-burden rows centered on OR≈1.3 straddling 1.

## 6. Reproducibility

```
python3 case_studies/h_pylori_cagpai/scripts/run_metadata_association.py
```

Regenerates `merged_genome_table.tsv`, `association_tests.tsv`, and all
figures from the raw collaborator files. Environment: python3 with pandas
2.3.3, scipy 1.13.1, statsmodels 0.14.6, matplotlib, openpyxl 3.1.5.
