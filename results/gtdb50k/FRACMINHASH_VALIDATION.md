# FracMinHash scale sweep validation

## Summary

We compared Syn2b `inverted_fraction` estimates from FracMinHash landmark sources
at scales 250, 750, 1582, 2000, and 6000 against the four-enzyme panel baseline,
using dnadiff-derived inverted fraction as truth on the held-out GTDB-R207
43,334-pair benchmark. The raw (majority-corrected) inverted fraction is used
throughout because it is the validated estimator.

## Accuracy by scale

| scale | n | MAE | RMSE | Pearson r | mean shared_tags | median shared_tags |
|---|---|---:|---:|---:|---:|---:|
| four | 43312 | 0.0370 | 0.0555 | 0.935 | 561 | 314 |
| fmh250 | 43334 | 0.0323 | 0.0477 | 0.951 | 1,266 | 761 |
| fmh750 | 43329 | 0.0399 | 0.0579 | 0.930 | 423 | 254 |
| fmh1582 | 43251 | 0.0488 | 0.0704 | 0.902 | 201 | 120 |
| fmh2000 | 43177 | 0.0528 | 0.0759 | 0.889 | 159 | 95 |
| fmh6000 | 41319 | 0.0784 | 0.1106 | 0.798 | 55 | 34 |

## Error-variance model

Fitted `Var(residual) = c0 + c1 · p(1-p)/m`, where `p` is dnadiff inverted
fraction and `m` is the number of shared landmarks.

| scale | c0 | c1 | r² | n |
|---|---|---:|---:|---:|---:|
| four | 0.00055 | 1.484 | 0.276 | 43312 |
| fmh250 | 0.00032 | 3.834 | 0.177 | 43334 |
| fmh750 | 0.00031 | 1.977 | 0.237 | 43329 |
| fmh1582 | 0.00054 | 1.332 | 0.261 | 43251 |
| fmh2000 | 0.00059 | 1.239 | 0.264 | 43177 |
| fmh6000 | -0.00024 | 1.174 | 0.284 | 41319 |

The irreducible variance `c0` is near zero for all scales, as expected for a
ratio that is bounded and corrected. The sampling coefficient `c1` is close to
the theoretical 1.504 at the four-enzyme density and at very sparse scales
(1582–6000), but it is elevated at the two densest FracMinHash scales (250 and
750). This suggests that when FracMinHash landmarks are much denser than the
enzyme panel, the binomial independence assumption breaks down and landmarks
begin to carry correlated information.

## Interpretation

1. **Accuracy improves with landmark density up to a point.** FracMinHash-250
gives the lowest MAE (0.032) and highest correlation (r = 0.951), beating the
four-enzyme panel (MAE 0.037). However, the computational cost and storage also
increase with density.

2. **The four-enzyme panel and FracMinHash-750 are density-matched** (~400–600
shared landmarks on a typical genome). Their accuracy is similar (MAE 0.037 vs
0.040), confirming that FracMinHash can substitute for enzyme digestion when the
scale is chosen to match the desired landmark density.

3. **Very sparse FracMinHash (6000) degrades gracefully** but is still usable
for coarse screens (MAE 0.078, r = 0.798). This is the regime where the method
approaches a k-mer-sketch-like sparsity.

4. **The error-model constant is not universal.** The elevated `c1` at fmh250
and fmh750 means the published standard-error formula should be recalibrated if
a much denser FracMinHash scale is adopted. For a drop-in replacement of the
four-enzyme panel, scales around 1500–2000 give the closest error-model match.

## Recommendation

- For best accuracy and a stable error model: use the four-enzyme panel or
  FracMinHash scale ~1500–2000.
- For maximum accuracy regardless of error-model fit: FracMinHash-250.
- For very fast coarse screens: FracMinHash-6000, with the understanding that
  inverted-fraction precision is reduced.

## Files

- `results/gtdb50k/inverted_fraction_truth_fmh{250,750,1582,2000,6000}.tsv`
- `results/gtdb50k/syn2b_inverted_fraction_50k_fmh{250,2000,6000}.tsv` (HPC)
- `scripts/analyze_fmh_sweep.py`
- `scripts/gtdb50k/validate_inverted_fraction_truth.py`
