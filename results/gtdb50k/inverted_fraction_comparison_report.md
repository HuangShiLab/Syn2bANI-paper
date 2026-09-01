# Comparison of dnadiff and Syn2b inverted fractions

Pairs: high-ANI GTDB test set, n = 2376

## Overall relationship
- Pearson r = 0.4507
- Spearman rho = 0.2598
- Syn2b = 0.3609 * dnadiff + 0.2336
- R² = 0.2031

## Unsaturated subset (dnadiff_inverted_fraction <= 0.5)
- n = 1371
- Pearson r = 0.9946
- Syn2b = 1.0007 * dnadiff + -0.0025
- Within the unsaturated range the intercept is essentially zero, consistent with both metrics being
  length-weighted ratios invariant to fragmentation. The overall intercept is driven by Syn2b saturation
  at 0.5 (MATH_REVIEW.md §7).

## Partial correlations (controlling for confounders)
- raw r = 0.4507
- partial | anim_ani = 0.4450
- partial | anim_ani + syn2b_observable_fraction = 0.4362
- Correlation does not improve when fragmentation (observable_fraction) is controlled, as predicted for a ratio metric.

## Bland-Altman
- mean difference (Syn2b - dnadiff) = -0.0663
- median difference = -0.0121
- SD of difference = 0.1422

## Interpretation
The two inverted-fraction metrics agree almost one-to-one when Syn2b is below its 0.5 saturation ceiling.
Above 0.5 dnadiff reports higher values because Syn2b flips its majority frame. This is the expected
behaviour, not a failure of the invariance argument. The count-based breakpoint comparison (SV_REANALYSIS.md)
shows a 290-unit intercept against dnadiff; the ratio comparison shows no intercept in the unsaturated range.
