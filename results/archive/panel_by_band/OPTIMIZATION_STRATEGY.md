# Syn2bANI GTDB-R207 optimization strategy

## Current benchmark result (vs ANIm truth, 2074 pairs)

| ANI band   | Syn2bANI default MAE (%) | skani MAE (%) | Greedy panel MAE (%) | Calibrated* MAE (%) |
|------------|--------------------------|---------------|----------------------|---------------------|
| 0.80–0.85  | 2.21                     | 1.74          | 2.12 (5 enzymes)     | 0.94                |
| 0.85–0.90  | 1.69                     | 0.85          | 1.66 (6 enzymes)     | 0.73                |
| 0.90–0.95  | 0.99                     | 0.45          | 0.89 (4 enzymes)     | 0.78                |
| 0.95–0.99  | 0.74                     | 0.34          | 0.68 (4 enzymes)     | 1.08                |
| **All**    | **1.54**                 | **0.91**      | **1.48**             | **0.81**            |

\*Ridge regression band-holdout calibration using only Syn2bANI-internal
features (`s2b_ani`, `s2b_ani_uniform`, `s2b_af_q/r`, `s2b_std_err`,
`s2b_retention`, `s2b_n_anchors/chains/tags`). No skani feature.

## Key findings

1. **Calibration is the biggest immediate win.**
   - Raw default output: MAE 1.54%, bias +0.55%.
   - Linear calibration (band holdout): **MAE 0.81%**, bias −0.06%.
   - This beats skani overall (0.91%) without using any skani information.
   - If `skani_ani` is included as an optional feature, MAE drops to **0.68%**.

2. **The diagnostic flag is currently misleading.**

   | Flag              | n    | Syn2bANI MAE (%) | Syn2bANI bias (%) | skani MAE (%) |
   |-------------------|------|------------------|-------------------|---------------|
   | `ok`              | 832  | 1.78             | **+1.62**         | 0.85          |
   | `INCONSISTENT`    | 1137 | 1.19             | +0.01             | 0.88          |
   | `BELOW_DETECTION` | 105  | 3.52             | −2.19             | 1.55          |

   Pairs labeled `ok` are *more* biased than `INCONSISTENT` pairs. The
   ok/inconsistent gate is not capturing the dominant error source and
   should be revisited.

3. **Error is strongly correlated with signal-strength features.**
   - `s2b_af_q` and `s2b_af_r` correlate negatively with absolute error
     (r ≈ −0.33): low alignment fraction → larger error.
   - `s2b_std_err` correlates positively with absolute error (r ≈ +0.30).
   - `s2b_retention` correlates positively with signed error (r ≈ +0.33):
     high tag retention is associated with overestimation.

4. **Adding more enzymes has diminishing returns.**
   - Greedy panel on the full set stops improving after 7 enzymes
     (MAE 1.529%).
   - Per-band optimal panels use 4–6 enzymes.
   - A smaller panel would reduce runtime without large accuracy loss,
     especially after calibration.

## Recommended optimization path

### 1. Embed a calibration model in `syn2bani ani` (highest priority)
- Train a small linear or GBRT model on a GTDB-R207 sample with ANIm truth.
- Features: the Syn2bANI-internal quantities above.
- Use strict band-stratified CV (and later an independent oral/gut set) to
  avoid overfitting.
- Expected target: **MAE ≈ 0.8% overall**, i.e. skani-level accuracy.

### 2. Fix or replace the `ok`/`INCONSISTENT` flag
- The current LRT-based heterogeneity test is not the right gate for ANI
  accuracy.
- Consider using the calibration residual / `s2b_std_err` / AF to assign a
  reliability score (e.g. high-confidence / medium / below-detection).

### 3. Offer smaller default enzyme panels for speed
- For example, per-band greedy panels already approach the best accuracy.
- After calibration, a 4–5 enzyme panel is likely both fast and accurate.
- Benchmark runtime vs enzyme count to quantify the trade-off.

### 4. Handle `BELOW_DETECTION` explicitly
- These 105 pairs have very low shared fraction and large negative bias.
- Better to report "no reliable ANI" or a lower bound than a point estimate.

### 5. Validate on independent oral/gut genomes
- Before shipping the calibration, confirm it generalizes outside GTDB-R207.

## Files

- `results/panel_by_band/eval_pairs.tsv` – merged truth + predictions.
- `results/panel_by_band/evaluation_summary.json` – per-band/per-phylum MAE.
- `results/panel_by_band/error_driver_report.json` – correlations and
  calibration probe results.
