# Single-enzyme BcgI validation against the four-enzyme panel and dnadiff truth

**Date:** 2026-09-01  
**Reproduce:**
```bash
python3 scripts/gtdb50k/validate_bcgI_error_model.py \
  --four results/gtdb50k/syn2b_inverted_fraction_50k.tsv \
  --bcgi results/gtdb50k/syn2b_inverted_fraction_50k_bcgI.tsv \
  --out results/gtdb50k/bcgI_error_validation.tsv \
  --plot results/gtdb50k/fig_bcgI_error_validation.png

python3 scripts/gtdb50k/validate_inverted_fraction_truth.py \
  --truth results/gtdb50k/dnadiff_inverted_fraction.tsv \
  --syn2b results/gtdb50k/syn2b_inverted_fraction_50k_bcgI.tsv \
  --label BcgI \
  --out results/gtdb50k/inverted_fraction_truth_bcgI.tsv
```

All comparisons use the **fixed-reference `raw_inverted_fraction`**, which is the
metric that agrees with dnadiff across the full [0,1] range. The majority-frame
`inverted_fraction` saturates above 0.5 and is not comparable to dnadiff there.

---

## 1. Tag yield: one enzyme gives ~47% of the four-enzyme landmarks

| panel | median shared tags | mean shared tags | vs four-enzyme |
|---|---:|---:|---:|
| four-enzyme | 337 | 585 | 1.00 |
| BcgI | 142 | 295 | 0.47 (median), 0.42 (mean) |

The four-enzyme panel does **not** quadruple tag count relative to BcgI; the
median ratio is 0.47. This is expected because the four enzymes have overlapping
recognition logic and because shared sites are conditional on both genomes
retaining the site.

## 2. BcgI vs four-enzyme: the panel is internally consistent

- n = 41,485 pairs with both estimates
- Pearson r = **0.897**
- `BcgI = 1.012 × four-enzyme − 0.006` (R² = 0.805)
- Mean difference = −0.0001; SD of difference = **0.0763**
- Predicted SE of difference from the four-enzyme error model = **0.0835**

The observed spread is **slightly smaller** than the conservative model predicts
(z SD = 0.667, 1.4% |z| > 2). The two panels agree to within the error model.

## 3. BcgI vs dnadiff truth: accurate but noisier than the panel

| metric | four-enzyme | BcgI |
|---|---:|---:|
| n | 43,312 | 41,485 |
| Pearson r | **0.9355** | **0.8303** |
| slope | 1.0039 | 0.9989 |
| intercept | −0.0024 | −0.0000 |
| RMSE | 0.0555 | **0.0961** |
| MAE | 0.0370 | 0.0576 |
| bias | −0.0005 | −0.0006 |

BcgI is **unbiased** (slope ≈ 1, intercept ≈ 0) but roughly **1.7× noisier**
than the four-enzyme panel.

## 4. The error model is panel-specific, not transferable

The fitted model
```
SE = sqrt(1.504 × p(1−p) / m + 0.0205²)
```
was trained on the four-enzyme panel. Applied to BcgI:

| panel | mean predicted SE | observed RMSE | z SD |
|---|---:|---:|---:|
| four-enzyme | 0.0463 | 0.0555 | **1.08** |
| BcgI | 0.0697 | 0.0961 | **2.88** |

For the four-enzyme panel the model is well calibrated (z SD ≈ 1). For BcgI it
is **too optimistic**: the actual error is ~1.4× the predicted SE on average, and
the z-score distribution is over-dispersed.

The under-prediction is strongest at very low shared-tag counts (median m = 8:
observed RMSE 0.236 vs predicted 0.173). At high tag counts the model is closer.
This suggests that a single-enzyme library not only samples fewer landmarks but
also has a higher **detection floor** for small events, which the current model
captures only through m.

## 5. Implications for the paper

1. **The four-enzyme panel is justified.** It does not merely increase tag count;
   it reduces the method floor relative to a single-enzyme design.
2. **The error model is useful within a panel but should not be used to predict
   the performance of a different enzyme set without recalibration.** A design
   rule that relies only on m under-predicts single-enzyme error.
3. **BcgI alone is still informative.** At r = 0.83, unbiased, it gives a
   length-weighted SV proxy at a fraction of the wet-lab cost. It is a viable
   option when only one enzyme is available, but with ~70% higher RMSE.

## 6. Next step

Refit the error model separately to BcgI (and ideally to each enzyme) to obtain
panel-specific constants `c_panel` and `sigma0_panel`. That would turn the
formula into a genuine design rule for arbitrary enzyme panels.
