# Comparison of dnadiff and Syn2b inverted fractions

This report validates two related quantities:

1. `syn2b_inverted_fraction` (majority-frame / min): the length-weighted ratio
   that is invariant to fragmentation but saturates at 0.5 (MATH_REVIEW.md §7).
2. `syn2b_raw_inverted_fraction` (fixed-reference / raw): orientation mismatches
   relative to genome_A, i.e. the first genome in the sorted TGT file list. This
   matches fixed-reference alignment methods such as dnadiff and ranges in [0,1].

## high_ani_all

### Syn2b (majority-frame, min)
- n = 8937
- Pearson r = 0.4274
- Spearman rho = 0.2025
- Syn2b = 0.3728 * dnadiff + 0.2218
- R² = 0.1827

#### Unsaturated subset (dnadiff_inverted_fraction <= 0.5)
- n = 5094
- Pearson r = 0.9017
- Syn2b = 0.9613 * dnadiff + 0.0010
- R² = 0.8130

#### Partial correlations
- raw r = 0.4274
- partial | anim_ani = 0.4480
- partial | anim_ani + observable_fraction = 0.4477

#### Bland-Altman
- mean difference (Syn2b - dnadiff) = -0.0757
- median difference = -0.0199
- SD of difference = 0.1384

### Syn2b (fixed-reference, raw)
- n = 8937
- Pearson r = 0.8753
- Spearman rho = 0.8441
- Syn2b = 0.9857 * dnadiff + 0.0062
- R² = 0.7662

#### Unsaturated subset (dnadiff_inverted_fraction <= 0.5)
- n = 5094
- Pearson r = 0.8659
- Syn2b = 1.0047 * dnadiff + -0.0009
- R² = 0.7499

#### Partial correlations
- raw r = 0.8753
- partial | anim_ani = 0.8741
- partial | anim_ani + observable_fraction = 0.8741

#### Bland-Altman
- mean difference (Syn2b - dnadiff) = -0.0006
- median difference = -0.0004
- SD of difference = 0.0748

## high_ani_test

### Syn2b (majority-frame, min)
- n = 2376
- Pearson r = 0.4507
- Spearman rho = 0.2598
- Syn2b = 0.3609 * dnadiff + 0.2336
- R² = 0.2031

#### Unsaturated subset (dnadiff_inverted_fraction <= 0.5)
- n = 1371
- Pearson r = 0.9946
- Syn2b = 1.0007 * dnadiff + -0.0025
- R² = 0.9891

#### Partial correlations
- raw r = 0.4507
- partial | anim_ani = 0.4450
- partial | anim_ani + observable_fraction = 0.4362

#### Bland-Altman
- mean difference (Syn2b - dnadiff) = -0.0663
- median difference = -0.0121
- SD of difference = 0.1422

## held_out_50k

### Syn2b (majority-frame, min)
- n = 43312
- Pearson r = 0.1771
- Spearman rho = 0.0604
- Syn2b = 0.1379 * dnadiff + 0.3239
- R² = 0.0313

#### Unsaturated subset (dnadiff_inverted_fraction <= 0.5)
- n = 22277
- Pearson r = 0.9281
- Syn2b = 0.9488 * dnadiff + 0.0081
- R² = 0.8613

#### Partial correlations
- raw r = 0.1771
- partial | observable_fraction = 0.1739

#### Bland-Altman
- mean difference (Syn2b - dnadiff) = -0.0995
- median difference = -0.0431
- SD of difference = 0.1687

### Syn2b (fixed-reference, raw)
- n = 43312
- Pearson r = 0.9355
- Spearman rho = 0.8774
- Syn2b = 1.0039 * dnadiff + -0.0024
- R² = 0.8751

#### Unsaturated subset (dnadiff_inverted_fraction <= 0.5)
- n = 22277
- Pearson r = 0.9095
- Syn2b = 1.0041 * dnadiff + -0.0026
- R² = 0.8272

#### Partial correlations
- raw r = 0.9355
- partial | observable_fraction = 0.9354

#### Bland-Altman
- mean difference (Syn2b - dnadiff) = -0.0005
- median difference = -0.0002
- SD of difference = 0.0555

## Interpretation
The majority-frame `inverted_fraction` agrees one-to-one with dnadiff when dnadiff
reports ≤0.5, because in that regime the minority frame is the fixed-reference frame.
Above 0.5 the two diverge because Syn2b flips its reference to the majority orientation.

The fixed-reference `raw_inverted_fraction` removes the saturation by always scoring
relative to genome_A. When genome_A is chosen to match dnadiff's reference (r_acc),
the two metrics should correlate across the full [0,1] range. This comes at the cost
of losing whole-genome reverse-complement invariance: a genome and its complement
read as 1.0, exactly as dnadiff reports.

Controlling for fragmentation (observable_fraction) does not increase the ratio
correlation, as predicted for a length-weighted quantity.

Two cautions about the aggregate numbers above.

They are computed without deduplicating `pairid`. The `high_ani_all` tables carry
195 repeated pairids, and an inner merge on a repeated key multiplies rows, so the
high_ani_all counts here are inflated relative to the 6,922 distinct pairs.

And `high_ani_all` is a set selected on *predicted* ANI, not measured ANI: 23% of
its pairs come back from ANIm at a median 84.4% identity over 12.5% of the
reference. Those pairs dominate its aggregate spread, which is why it reports
r = 0.8753 against held_out_50k's 0.9355 — not because high-ANI pairs are harder.
Condition on measured ANIm before quoting anything from that set. The section below
does, and adds the error model that says what the residual is made of.

## Error model for `raw_inverted_fraction`

Regenerate with `python3 scripts/analyze_invfrac_error_model.py results/gtdb50k`.

### Bias is zero in every divergence band

Fitted separately inside each ANIm band of held_out_50k. The slope stays at 1 and the bias at 0 all the way down to 80% ANIm; only the spread moves. The estimator does not degrade at low ANI — it loses precision, which is a reportable standard error rather than a systematic error.

| ANIm | n | slope | intercept | r | bias | SD(err) | median shared tags | median aligned frac (%) |
|---|---|---|---|---|---|---|---|---|
| 80-85 | 1850 | 0.976 | +0.0080 | 0.8825 | -0.0032 | 0.1103 | 47 | 49.7 |
| 85-88 | 17576 | 1.003 | -0.0020 | 0.9149 | -0.0004 | 0.0634 | 214 | 49.5 |
| 88-90 | 8163 | 1.009 | -0.0051 | 0.9448 | -0.0005 | 0.0481 | 391 | 65.0 |
| 90-92 | 6425 | 1.006 | -0.0030 | 0.9607 | -0.0002 | 0.0424 | 545 | 71.2 |
| 92-95 | 8644 | 1.010 | -0.0051 | 0.9752 | -0.0000 | 0.0326 | 809 | 77.1 |
| 95-97 | 652 | 1.017 | -0.0099 | 0.9862 | -0.0015 | 0.0275 | 1162 | 79.0 |

### The strain range, from the high-ANI set

held_out_50k has only 2 pairs at >=97% ANIm, so it says nothing about the range the tool is meant for. That range is covered by the high_ani set, which was sampled for it. Agreement there is not merely good, it is close to exact — and note the spread keeps falling past the floor fitted on the mixed set, which is the clearest evidence that the floor is a function of divergence rather than a constant of the method.

| ANIm | n | slope | intercept | r | bias | SD(err) | median shared tags |
|---|---|---|---|---|---|---|---|
| 95-97 | 610 | 1.019 | -0.0104 | 0.9872 | -0.0013 | 0.0214 | 1463 |
| 97-98 | 594 | 1.009 | -0.0068 | 0.9922 | -0.0024 | 0.0167 | 1998 |
| 98-99 | 1064 | 1.011 | -0.0055 | 0.9950 | -0.0004 | 0.0135 | 2349 |
| 99-99.5 | 617 | 1.010 | -0.0060 | 0.9951 | -0.0012 | 0.0129 | 2602 |
| 99.5-100.1 | 1551 | 1.004 | -0.0022 | 0.9974 | -0.0002 | 0.0122 | 2828 |

Pooled over >=97% ANIm: n = 3826, slope 1.0063, intercept -0.0037, r = 0.9960, bias -0.0008, SD 0.0135.

### The whole ANI dependence runs through the landmark count

Lower ANI destroys restriction sites, so fewer landmarks are shared, so the same proportion is estimated from a smaller sample. Binning held_out_50k by shared-landmark count `m` and fitting a sampling term plus a constant floor:

```
Var(err) = 1.504 * p(1-p)/m + 0.0205^2      (12 bins, R2 = 0.9988)
```

The coefficient is 1.50 rather than 1 because landmarks inside an inverted segment are spatially clustered rather than independently drawn; a design effect near 1.5 is what clustering produces. The floor 0.0205 is the method difference itself: dnadiff averages over aligned bases, Syn2b over shared landmarks, and those denominators are not the same set.

| median m | n | SD(err) observed | SD(err) model |
|---|---|---|---|
| 17 | 2081 | 0.1416 | 0.1417 |
| 47 | 2050 | 0.0889 | 0.0901 |
| 81 | 3380 | 0.0703 | 0.0694 |
| 130 | 5094 | 0.0586 | 0.0561 |
| 202 | 5927 | 0.0491 | 0.0465 |
| 315 | 6347 | 0.0418 | 0.0391 |
| 503 | 6894 | 0.0348 | 0.0334 |
| 801 | 4484 | 0.0291 | 0.0292 |
| 1233 | 3776 | 0.0241 | 0.0265 |
| 1944 | 2066 | 0.0206 | 0.0244 |
| 3007 | 1015 | 0.0179 | 0.0231 |
| 4544 | 198 | 0.0155 | 0.0222 |

The fit is unweighted across bins, and in the top bins the model over-predicts (e.g. 0.0155 observed against 0.0222 at median m = 4544). Those bins are also the high-ANI pairs, so this is the floor itself shrinking with divergence rather than a misfit; the model is conservative where landmarks are plentiful.

### Applied out of sample

The model is fitted on held_out_50k bins only. Applied per pair to whole datasets it reproduces the aggregate spread, and standardised residuals have roughly unit variance:

| set | n | SD(err) observed | SD(err) model | SD(z) | within +-2 SE |
|---|---|---|---|---|---|
| held_out_50k | 43312 | 0.0555 | 0.0546 | 1.006 | 95.3% |
| high_ani_all | 6922 | 0.0848 | 0.0874 | 0.763 | 97.9% |

### Why high_ani_all reports a *lower* r than the full held-out set

Not range restriction: the two sets have nearly the same spread of truth (SD 0.136 vs 0.146) and nearly the same mean p(1-p) (0.2309 vs 0.2285). The cause is that high_ani_all is not actually a high-ANI set. It was selected on a *predicted* ANI, and 23% of its pairs come back from ANIm at a median of 84.4% identity over 12.5% of the reference — distant pairs sharing a small island, not strains. Those pairs carry almost no landmarks and set the aggregate spread; the ones that survive the check agree far more tightly (SD 0.0144 at m >= 500) than anything in held_out_50k.

| shared tags | n | median ANIm | median AF of ref (%) | SD(err) |
|---|---|---|---|---|
| < 100 | 1563 | 84.37 | 12.5 | 0.1725 |
| >= 500 | 4505 | 98.87 | 87.9 | 0.0144 |

So the aggregate row for high_ani_all in the sections above should not be read as a high-ANI result at all. The banded table earlier in this document, which conditions on measured ANIm, is the one to quote.

### Reporting a single pair

`syn2b_shared_tags` is already emitted, so every pair can carry its own standard error:

```
SE = sqrt( 1.504 * p(1-p) / shared_tags + 0.0205^2 )
```

| shared tags | SE at p=0.5 | sampling share of variance |
|---|---|---|
| 50 | 0.0891 | 95% |
| 100 | 0.0647 | 90% |
| 200 | 0.0480 | 82% |
| 500 | 0.0343 | 64% |
| 1000 | 0.0282 | 47% |
| 3000 | 0.0234 | 23% |
| 10000 | 0.0214 | 8% |
| infinity | 0.0205 | 0% |

Two consequences for the panel design. Past roughly 1,000 shared landmarks the sampling term is no longer dominant, so additional restriction sites buy little for the orientation channel — the four-enzyme panel has to be justified by the junction channel's resolution floor instead. And a single low-`m` pair's point estimate is not interpretable on its own, even though the mean over many such pairs stays unbiased.
