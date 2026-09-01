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
