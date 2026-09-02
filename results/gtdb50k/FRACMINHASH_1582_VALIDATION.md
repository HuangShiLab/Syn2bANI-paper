# FracMinHash scale 1582 validation (GTDB 50k held-out)

**Date:** 2026-09-02  
**Reproduce:**
```bash
python3 scripts/gtdb50k/validate_inverted_fraction_truth.py \
  --truth results/gtdb50k/dnadiff_inverted_fraction.tsv \
  --syn2b results/gtdb50k/syn2b_inverted_fraction_50k_fmh1582.tsv \
  --label FracMinHash-1582 \
  --out results/gtdb50k/inverted_fraction_truth_fmh1582.tsv
```

All comparisons use the fixed-reference `raw_inverted_fraction`, which is
comparable to dnadiff across the full [0,1] range.

---

## 1. Scale 1582 matches BcgI's landmark count

On a 4.5 Mb genome, `--scale 1582` gives ~2,870 FracMinHash landmarks, against
~2,872 for BcgI. The two runs therefore isolate the effect of landmark
*distribution* while holding `m` fixed.

| panel | median shared tags | mean shared tags |
|---|---:|---:|
| BcgI | 142 | 295 |
| FracMinHash-1582 | 105 | 254 |

FracMinHash-1582 is slightly sparser on these pairs, but the same order of
magnitude as BcgI.

## 2. Accuracy vs dnadiff truth

| panel | n | Pearson r | slope | intercept | RMSE | MAE | bias |
|---|---:|---:|---:|---:|---:|---:|---:|
| four-enzyme | 43,312 | **0.9355** | 1.0039 | −0.0024 | **0.0555** | 0.0370 | −0.0005 |
| FracMinHash-1582 | 43,251 | **0.9020** | 1.0077 | −0.0041 | **0.0704** | 0.0488 | −0.0004 |
| BcgI | 41,485 | **0.8303** | 0.9989 | −0.0000 | **0.0961** | 0.0576 | −0.0006 |

FracMinHash-1582 sits between the four-enzyme panel and BcgI, despite having
roughly the same landmark count as BcgI.

## 3. Error-model calibration: the key test

The fitted four-enzyme error model

```
SE = sqrt(1.504 * p(1-p) / m + 0.0205^2)
```

was trained on four-enzyme landmarks. Applied to other panels:

| panel | mean predicted SE | observed RMSE | z SD |
|---|---:|---:|---:|
| four-enzyme | 0.0463 | 0.0555 | **1.08** |
| FracMinHash-1582 | 0.0634 | 0.0704 | **1.18** |
| BcgI | 0.0697 | 0.0961 | **2.88** |

FracMinHash-1582 is well-calibrated (z SD ≈ 1.2), while BcgI is substantially
over-dispersed (z SD ≈ 2.9). Both have similar `m`, so the difference is not
landmark count. The model's `1/m` term is correct; what breaks BcgI is the
*distribution* of its landmarks — restriction-site clustering and Hamming-1
near-duplicates — which FracMinHash avoids by using a Poisson-uniform k-mer
sketch.

## 4. Implications

1. **The four-enzyme error model is a property of the estimator, not of the
   enzymes.** It transfers cleanly to FracMinHash landmarks.
2. **BcgI's under-performance is not a small-m problem.** At matched density,
   FracMinHash reaches r = 0.90 where BcgI reaches 0.83, with a much better
   calibrated residual.
3. **The design rule is "use uniform landmarks", not "recalibrate per panel".**
   Restriction sites cluster on motifs and leave large gaps; FracMinHash fills
   those gaps uniformly.

## 5. Next step

Run `--scale 750` (four-enzyme-matched density, ~6,030 landmarks on a 4.5 Mb
genome) to check whether the same constants hold at higher density and to make
a direct four-enzyme vs FracMinHash comparison.
