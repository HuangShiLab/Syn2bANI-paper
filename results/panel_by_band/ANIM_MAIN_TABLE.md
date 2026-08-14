# ANIm main table — panel rescore + 4-enzyme ridge calibration

Date: 2026-08-13. All work local on the Mac Studio; no git operations.
Binary: `/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani`
(syn2bani 0.1.0, main @ 69ce9f4). Truth: 2,074 GTDB-R207 pairs with dnadiff
ANIm (`results/panel_by_band/eval_pairs.tsv`, `anim_ani` in percent).

## 1. Panel rescoring from stored strata (no new compute)

Commands:

```
syn2bani panel --strata results/gtdb_r207_100k_strata.tsv \
    --truth results/panel_by_band/truth_all.tsv --panels 'BcgI,AlfI,AloI,FalI' \
  > results/panel_by_band/rescore_4e_current.txt
syn2bani panel --strata ... --panels \
    'AlfI,AloI,BcgI,BplI,BsaXI,Bsp24I,CjeI,CjePI,FalI,PpiI,PsrI' \
  > results/panel_by_band/rescore_11e_full.txt
syn2bani panel --strata ... --panels \
    'AlfI,AloI,BplI,BsaXI,CjeI,CjePI,FalI;CjePI;CjeI,CjePI;AlfI,BplI,CjeI,CjePI' \
  > results/panel_by_band/rescore_intermediate.txt
```

`panel` prints summary MAE/bias only (no per-pair output). Results:

| panel | MAE | bias | n |
|---|---|---|---|
| BcgI,AlfI,AloI,FalI (current default) | 1.8465 | +0.6624 | 2074 |
| CjePI | 1.6058 | +0.6436 | 2074 |
| CjeI,CjePI | 1.5687 | +0.5105 | 2074 |
| AlfI,BplI,CjeI,CjePI | 1.5426 | +0.5158 | 2074 |
| AlfI,AloI,BplI,BsaXI,CjeI,CjePI,FalI (greedy-best-7) | 1.5289 | +0.4862 | 2074 |
| all 11 enzymes | 1.5437 | +0.5453 | 2074 |

**Provenance confirmation: YES.** The 11-enzyme strata rescore gives
MAE 1.5437 / bias +0.5453, identical to the numbers computed from
`eval_pairs.tsv`'s `s2b_ani` column (1.5437 / +0.5453, this analysis). The old
`eval_pairs` s2b_ani values came from the same 11-enzyme run that wrote the
strata file.

**Caveat / surprise:** the strata rescore of the *current* 4-enzyme panel
(1.847, n = 2074) does NOT match the current v8-final matrix numbers
(gamma MAE 2.874, n = 1,969 finite). The strata file predates the v8 estimator
changes: it records per-enzyme statistics from the older run, so `panel`
rescores that older estimator's behavior. The v8 matrix run (current binary,
gamma estimator, BELOW_DETECTION masking) is noticeably worse uncalibrated
than the same 4 enzymes under the old estimator (2.874 vs 1.847 raw MAE).
This gap — not the enzyme panel itself — is what the ridge calibration
recovers (below).

## 2. Band-holdout ridge calibration on current 4-enzyme features

Script: `scripts/anim_main_table_4e.py` (run: `python3
scripts/anim_main_table_4e.py`). Methodology follows
`scripts/analyze_error_drivers.py::calibration_experiment` and
`results/panel_by_band/OPTIMIZATION_STRATEGY.md`:

- Features (Syn2bANI-internal only, from
  `results/matrix_gtdb_r207_100k_v8_final.tsv` joined to the 2,074 truth pairs
  on query/ref assembly accessions): `s2b_ani`, `s2b_ani_uniform`, `s2b_af_q`,
  `s2b_af_r`, `s2b_std_err`, `s2b_retention`, `s2b_n_anchors`, `s2b_n_chains`,
  `s2b_n_tags`.
- Band-holdout CV: train on 3 bands, predict the held-out band, rotate over
  all 4 bands. Median imputation + StandardScaler + RidgeCV
  (alphas = 1e-3..100, inner cv = 5), all fitted on the training bands only.
  (Difference vs the original probe: it did not standardize inside the CV
  loop; we do. Selected alphas: 0.001 / 1 / 10 / 10 per held-out band.)
- **BELOW_DETECTION handling:** the 105 pairs with NaN 4e `s2b_ani` (all
  flagged BELOW_DETECTION) are EXCLUDED from the 4e gamma/uniform/ridge rows
  (ridge needs `s2b_ani` as its primary feature; these pairs carry no point
  estimate). 4e n = 1,969. They remain in the 11e/skani rows (n = 2,074).
- The old 11-enzyme feature matrix is not available (the matrix file now
  holds 4e columns), so the ridge was fit for 4e only; the old 11e embedded
  calibration (`s2b_ani_cal`) is reported as-is.

Result: ridge-CV overall MAE **0.906** (bias −0.037, RMSE 1.307, r 0.913),
vs raw 4e gamma 2.874 and skani 0.906. Per band: 0.876 / 0.852 / 0.938 /
1.375 (80–85 / 85–90 / 90–95 / 95–99). Largest mean coefficients:
`s2b_retention` +2.85, `s2b_ani_uniform` −1.12, `s2b_af_r` +0.97 (details in
`ridge_cv_4e_report.json`). Calibration fully closes the raw 4e gap and
matches skani overall; the 95–99% band (n = 72) is the weakest (MAE 1.375,
bias −1.19, r 0.134) — the model is trained without that band's regime and
compresses its small dynamic range.

## 3. Outputs

- `results/panel_by_band/anim_main_table.tsv` — consolidated table,
  method × band (plus band = "all"): n, MAE, RMSE, bias, r. Methods:
  syn2bani_4e_gamma, syn2bani_4e_uniform, syn2bani_4e_ridge_cv, syn2bani_11e,
  syn2bani_11e_cal, skani, FastANI_subset (363 pairs only, n marked).
- `results/panel_by_band/ridge_cv_preds_4e.tsv` — per-pair out-of-fold ridge
  predictions (1,969 rows).
- `results/panel_by_band/ridge_cv_4e_report.json` — CV settings, fold alphas,
  mean coefficients, exclusion counts.
- `results/panel_by_band/rescore_*.txt` — `syn2bani panel` outputs.
- `figures/report/fig7_anim_by_band.{png,pdf}` + caption in
  `figures/report/CAPTIONS.md` (`fig7_anim_by_band()` in
  `analysis/plot_report_figures.py`).

All pre-established control numbers reproduced exactly in this re-analysis:
4e gamma 2.874 (bias +2.601, n = 1,969; per-band 3.689 / 3.529 / 1.817 /
1.159), 4e uniform 3.566, 11e 1.544, 11e_cal 1.231, skani 0.906, FastANI
1.056 (n = 363).
