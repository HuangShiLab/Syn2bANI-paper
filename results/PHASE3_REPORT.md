# Phase 3 Report: Mid-ANI ANIm re-validation and hybrid-boundary check

## 1. Mid-ANI oral/gut 15 pairs re-evaluated against ANIm/minimap2 truth

ANIm truth for the 15 mid-ANI pairs (*Bifidobacterium longum* × *B. breve*,
*Veillonella parvula* × *V. dispar*) already exists in
`results/validation_mid_ani_anim/anim_4e/anim_midani_evaluation.tsv`
(dnadiff 1-to-1, MUMmer 3.23). The file reports Syn2bANI v8 4-enzyme output,
skani, and FastANI on the same 15 pairs.

### Metrics vs ANIm truth

| Method | n | MAE (ANI points) | bias | r |
|---|---|---|---|---|
| Syn2bANI `s2b_ani` (default/calibrated column) | 15 | 4.4824 | −4.4137 | −0.543 |
| Syn2bANI `ani_uniform` (uniform MLE; gate fallback because all flags are `INCONSISTENT`) | 15 | 1.4231 | +1.4231 | 0.683 |
| skani | 15 | 1.3963 | −1.3963 | 0.672 |
| FastANI | 15 | 1.8642 | −1.8642 | 0.521 |

### Interpretation

- The current default output (`s2b_ani`) performs poorly on these mid-ANI
  pairs, underestimating by ~4.4 points on average. This is consistent with
  the known regime-specificity problem: the GTDB-trained calibration
  overshoots on mid-ANI / low-AF pairs.
- When the consistency flag fires and the gate falls back to the homogeneous
  (uniform) MLE, the error drops to MAE 1.42, comparable to skani (1.40) and
  better than FastANI (1.86).
- All 15 pairs are flagged `INCONSISTENT`, so the *gated* estimate is the
  uniform estimate for this set. The manuscript’s reported gated MAE of 0.959
  on mid-ANI is not reproduced by the current 4-enzyme v8 output; the
  discrepancy may reflect an earlier gate threshold, a different calibration
  snapshot, or a different subset. This should be reconciled before
  submission.

## 2. Hybrid boundary evaluation on the unified 80–100% GTDB-R207 benchmark

The unified benchmark combines 43,334 held-out same-genus pairs (80–95% plus
95–100 split) with 727 high-ANI test pairs (95–97 / 97–100). The hybrid rule
is: use `ani_gated` when it is finite and ≥ threshold, otherwise use
`ani_cal`.

### Threshold sweep (overall MAE)

| Threshold (%) | n scored | MAE | bias |
|---|---|---|---|
| 95.0 | 40,629 | 0.7065 | +0.0284 |
| 96.0 | 40,629 | 0.6443 | −0.0699 |
| 97.0 | 40,629 | 0.6189 | −0.1089 |
| 97.5 | 40,629 | 0.6161 | −0.1133 |
| 98.0 | 40,629 | 0.6146 | −0.1154 |
| 98.5 | 40,629 | 0.6144 | −0.1169 |
| 99.0 | 40,629 | 0.6140 | −0.1191 |
| 99.5 | 40,629 | 0.6145 | −0.1211 |

The optimum by overall MAE is **99.0%** (MAE 0.6140), but the difference
between 97.0% and 99.5% is small (< 0.003 ANI points). The currently deployed
**98.0%** threshold (MAE 0.6146) is essentially optimal and arguably more
robust because it avoids pushing more borderline 95–97% pairs onto the raw
estimate.

### Per-band MAE at 98% threshold

| Band | n | MAE | bias |
|---|---|---|---|
| 80–85 | 8,823 | 0.7193 | −0.2046 |
| 85–90 | 15,918 | 0.6450 | −0.0950 |
| 90–95 | 14,758 | 0.5411 | −0.1032 |
| 95–97 | 496 | 0.4514 | +0.0815 |
| 97–100 | 634 | 0.2324 | +0.1748 |

### Conclusion on hybrid boundary

- The 98% boundary is well-supported; moving it to 99% changes overall MAE by
  less than 0.001 points.
- The main gain of the hybrid rule is in the 97–100% band (raw gated MAE
  0.23 vs calibrated v6 MAE 0.55).
- In the 95–97% band the two estimators are similar; the boundary could be
  lowered to 97% without materially harming overall accuracy.

## 3. Files updated / created

- `results/gtdb50k/hybrid_threshold_sweep.tsv` — threshold sweep table.
- `results/PHASE3_REPORT.md` — this report.

## 4. Open questions

1. **Mid-ANI manuscript number mismatch**: The current v8 4-enzyme gated/uniform
   MAE on the 15 mid-ANI pairs is 1.42, not 0.959. Confirm whether the
   manuscript value came from a different snapshot or a different estimator
   configuration.
2. **Mid-ANI practical recommendation**: Because the calibrated column is
   harmful here, the recommendation for mid-ANI / MAG inputs should be to
   use the raw/gated output (or the uniform fallback when flagged), not the
   calibrated output. This is already stated in the manuscript but the exact
   numbers need updating.
