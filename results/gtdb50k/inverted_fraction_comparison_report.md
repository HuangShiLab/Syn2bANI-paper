# Comparison of dnadiff and Syn2b inverted fractions

This report validates the length-weighted ratio invariance argument in MATH_REVIEW.md §7.

## high_ani_all
- n = 8937
- Pearson r = 0.4274
- Spearman rho = 0.2025
- Syn2b = 0.3728 * dnadiff + 0.2218
- R² = 0.1827

### Unsaturated subset (dnadiff_inverted_fraction <= 0.5)
- n = 5094
- Pearson r = 0.9017
- Syn2b = 0.9613 * dnadiff + 0.0010
- R² = 0.8130

### Partial correlations
- raw r = 0.4274
- partial | anim_ani = 0.4480
- partial | anim_ani + observable_fraction = 0.4479

### Bland-Altman
- mean difference (Syn2b - dnadiff) = -0.0757
- median difference = -0.0199
- SD of difference = 0.1384

## high_ani_test
- n = 2376
- Pearson r = 0.4507
- Spearman rho = 0.2598
- Syn2b = 0.3609 * dnadiff + 0.2336
- R² = 0.2031

### Unsaturated subset (dnadiff_inverted_fraction <= 0.5)
- n = 1371
- Pearson r = 0.9946
- Syn2b = 1.0007 * dnadiff + -0.0025
- R² = 0.9891

### Partial correlations
- raw r = 0.4507
- partial | anim_ani = 0.4450
- partial | anim_ani + observable_fraction = 0.4362

### Bland-Altman
- mean difference (Syn2b - dnadiff) = -0.0663
- median difference = -0.0121
- SD of difference = 0.1422

## held_out_50k
- n = 43312
- Pearson r = 0.1771
- Spearman rho = 0.0604
- Syn2b = 0.1379 * dnadiff + 0.3239
- R² = 0.0313

### Unsaturated subset (dnadiff_inverted_fraction <= 0.5)
- n = 22277
- Pearson r = 0.9281
- Syn2b = 0.9488 * dnadiff + 0.0081
- R² = 0.8613

### Partial correlations
- raw r = 0.1771
- partial | observable_fraction = 0.1744

### Bland-Altman
- mean difference (Syn2b - dnadiff) = -0.0995
- median difference = -0.0431
- SD of difference = 0.1687

## Interpretation
The two inverted-fraction metrics agree almost one-to-one when Syn2b is below its 0.5 saturation ceiling.
Above 0.5 dnadiff reports higher values because Syn2b flips its majority frame. This is the expected
behaviour, not a failure of the invariance argument. The count-based breakpoint comparison (SV_REANALYSIS.md)
shows a 290-unit intercept against dnadiff; the ratio comparison shows no intercept in the unsaturated range.

Controlling for fragmentation (observable_fraction) does not increase the ratio correlation, as predicted
for a length-weighted quantity. In contrast, controlling fragmentation increases the count-based breakpoint
correlation (SV_REANALYSIS.md).
