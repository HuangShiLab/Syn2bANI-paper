# Calibration v2 — retrained on current-binary features, externally validated

Date: 2026-08-14. Binary whose features were used: syn2bani 0.1.0, Syn2bANI
main @ `69ce9f4` (feature matrix `results/anim_truth_2074_v8current.tsv`).
Script: `scripts/calibration_v2.py` (run: `python3 scripts/calibration_v2.py`).
Truth: ANIm (dnadiff) on 2,074 GTDB-R207 pairs
(`results/panel_by_band/eval_pairs.tsv`), joined via
`results/anim_2074_acc2seqid.tsv` (seqid → assembly accession).

Protocol (unchanged from the July analysis, `scripts/anim_main_table_4e.py`):
band-holdout CV — train on 3 ANI bands, test on the held-out band, rotate;
per-fold median imputation + StandardScaler + RidgeCV (alphas 1e-3..100,
inner cv = 5). **No flag filtering** (the ok/INCONSISTENT flag is inverted on
GTDB); only the 21 rows with non-finite `ani` were dropped → n = 2,053.
`het_shape = inf` is treated as missing and median-imputed, matching the
Rust `predict()` non-finite → imputer-median path.

## 1. Feature-set comparison (band-holdout CV, MAE/bias/r vs ANIm)

Set A (base 9): `ani, ani_uniform, af_query, af_reference, std_err,
retention, n_anchors, n_chains, n_tags`.
Set B (expanded 18): A + `synteny_score, breakpoint_count, enzyme_spread,
enzyme_chi2, het_shape, ani_from_loss, ani_from_hist, max_block_anchors,
mean_block_anchors`.
Set C: Set B features with GradientBoostingRegressor (sklearn defaults).

| experiment | 80–85 | 85–90 | 90–95 | 95–99 | all (MAE) | bias | r |
|---|---|---|---|---|---|---|---|
| raw gamma (no calibration) | 3.571 | 3.529 | 1.818 | 1.159 | 2.881 | +2.411 | 0.706 |
| A base-9 ridge | 0.960 | 0.871 | 1.106 | 1.388 | 0.988 | −0.074 | 0.906 |
| **B expanded ridge** | **0.951** | **0.846** | **1.078** | **1.286** | **0.963** | −0.066 | **0.910** |
| C expanded GBRT | 1.059 | 1.120 | 1.446 | 1.347 | 1.224 | −0.247 | 0.887 |

(Full table with RMSE: `results/panel_by_band/calibration_v2_cv.tsv`.)

Answers:

- The new synteny/dispersion features **do improve out-of-band prediction**,
  but modestly: B beats A in all four bands, overall MAE 0.963 vs 0.988
  (−2.5% relative). Largest Set-B coefficients: `retention` +2.81,
  `ani_uniform` −1.67, `ani_from_loss` +1.11, `af_reference` +0.99,
  `synteny_score` +0.49.
- **Nonlinearity buys nothing at n ≈ 2k**: GBRT is worse than ridge overall
  (1.224 vs 0.963) and in every band. Consistent with the learning-curve
  plateau (~1,500 pairs) found earlier — the residual error is not
  expressible structure, it is noise/truth-floor.
- Note vs the July run: the 80–85 band now has n = 463 (was 380) because the
  current binary emits finite `ani` for 84 of the 105 BELOW_DETECTION-flagged
  pairs; per instructions only non-finite `ani` rows were dropped. Overall
  raw gamma MAE is 2.881 on n = 2,053 (July: 2.874 on n = 1,969) — the
  picture is unchanged.

## 2. External validation

Final models trained on ALL 2,053 pairs (`train_final_ridge`, RidgeCV with
alphas extended to 1000). Set B chosen for production by the pre-set rule
(B beats A by >2% relative overall; it wins every band).
`results/panel_by_band/linear_cal_v2.json` = Set B (18 features);
`results/panel_by_band/linear_cal_v2_setA.json` = Set A fallback.

### (a) Oral/gut, 100 same-species pairs vs FastANI (`data/oral_gut_validation_merged_v8.tsv`)

The July merged file has only Set-A-era columns; the 1,225 pairs were
therefore re-run with the current binary on the HPC
(`results/oral_gut_1225_v8current.tsv`, seqid-keyed, mapping
`results/oral_gut_1225_acc2seqid.tsv`; values byte-identical to July). Only
122/1,225 pairs have finite `ani` (88 ok + 34 INCONSISTENT; the other 1,103
are BELOW_DETECTION NaN — divergent pairs below the tag detection floor,
expected), and FastANI reports only the 100 same-species (`high`-label)
pairs among them, so the evaluable set is the same 100 pairs as before.
**Both models are now validated here.**

| method | n | MAE | bias | r |
|---|---|---|---|---|
| raw gamma | 100 | 0.552 | +0.421 | 0.980 |
| raw uniform | 100 | 1.165 | +1.165 | 0.975 |
| **calibrated v2 (Set A)** | 100 | **0.460** | +0.152 | **0.995** |
| calibrated v2 (Set B) | 100 | 0.636 | −0.095 | 0.927 |

**Set B does not transfer to oral/gut**: its in-distribution CV advantage
reverses out-of-distribution — MAE 0.636 vs Set A's 0.460, and r drops from
0.995 to 0.927. The degradation is concentrated on the cleanest pairs
(ok-flag: MAE 0.577 for B vs 0.274 for A; INCONSISTENT: 0.750 vs 0.820) —
the near-clonal oral/gut pairs sit far outside the GTDB training feature
distribution (synteny_score ~0.94 vs ~0.33, retention ~0.96 vs ~0.3), and
the 18-feature linear model extrapolates worse there than the 9-feature one.
Caveat: the reference is FastANI, which reads ~1 point low vs ANIm at
divergent bands; on these same-species pairs (ANI ≈ 97–99.5) that offset is
small, and the calibrated biases are +0.15 (A) / −0.09 (B). Set A
generalizes outside its GTDB training distribution (oral/gut species,
different lab pipeline).

### (b) Mid-ANI, 15 pairs vs ANIm truth (`results/validation_mid_ani_anim/anim/`)

Current-binary features, all Set-B columns present. This is the only
external ANIm-truth check.

| method | n | MAE | bias | r |
|---|---|---|---|---|
| raw gamma | 15 | 4.482 | −4.414 | −0.543 |
| raw uniform | 15 | 1.423 | +1.423 | 0.683 |
| calibrated v2 (Set A) | 15 | 1.343 | −0.887 | −0.483 |
| **calibrated v2 (Set B)** | 15 | **1.229** | **−0.845** | −0.206 |

Per-pair detail: `results/panel_by_band/calibration_v2_midani_pairs.tsv`.
Calibration removes the gamma overshoot (−4.4 → −0.85) and beats raw
uniform. The near-zero/negative r is expected and not alarming: all 15
truths sit in a 2.5-point window (87.6–90.2), so ranking power is
near-meaningless at n = 15; what matters is the error magnitude. All 15
pairs are flagged INCONSISTENT, so a flag-filtering user never sees these
numbers either way.

## 3. Overfitting / leakage assessment

- In-sample MAE 0.745 vs band-holdout CV MAE 0.963 — a normal train/test gap
  for ridge at n = 2,053; no runaway fit (selected alpha = 1.0, mid-grid).
- No leakage: imputer/scaler/ridge are fitted on the training bands only;
  held-out bands are entirely unseen. All features are tool-internal
  quantities, not functions of the ANIm truth. `ani_from_loss/hist` are the
  tool's own alternative estimates — legitimate internal features.
- GBRT (more capacity) doing *worse* out-of-band is itself evidence the
  ridge is not underfit; the learning curve says more of the same data would
  not help either.

## 4. Recommendation

**Ship Set A (`linear_cal_v2_setA.json`, 9 features).** The full evidence
(updated 2026-08-14 with the oral/gut Set-B re-run):

| evidence | Set A | Set B | winner |
|---|---|---|---|
| GTDB band-holdout CV (n = 2,053) | 0.988 | 0.963 | B, narrowly (−2.5%) |
| mid-ANI 15 pairs vs ANIm | 1.343 | 1.229 | B (n = 15, 2.5-point window — weak) |
| oral/gut 100 same-species vs FastANI | **0.460 / r 0.995** | 0.636 / r 0.927 | **A, clearly** |

Set B's CV win is small and does not survive distribution shift: on the
near-clonal oral/gut pairs — the tool's main operating regime — B's extra
synteny/dispersion features extrapolate badly (ok-flag MAE 0.577 vs A's
0.274; r 0.927 vs 0.995). B's mid-ANI win rests on n = 15 inside a
2.5-point truth window. Set A is also the zero-Rust-change option (§5): the
existing 9-feature `predict_from_result` vector matches
`linear_cal_v2_setA.json` exactly, so shipping is a pure JSON swap with no
test changes. Keep the per-band caveat — the 95–99% GTDB band stays
underpowered (n = 72, MAE 1.39 for A / 1.29 for B).

**Do not ship Set B as-is.** It stays documented (`linear_cal_v2.json`) as a
candidate for divergent-pair accuracy; revisit only if a larger external
ANIm set (n ≫ 15) confirms its mid-ANI advantage without the oral/gut
regression. A band-dependent choice (Set B below ~95% ANI, Set A above) is a
possible follow-up but needs more than 15 external pairs to justify.

## 5. Production spec — Rust changes required to ship

Model schema (`src/core/calibration.rs::LinearCalModel`): `name`,
`feature_names`, `means`, `scales`, `coefficients`, `intercept`,
`imputer_medians`, `training_n`, `training_mae` — the v2 JSONs above conform
(feature order = the order Rust must supply them, see below).

TSV ↔ `ChainAniResult` mapping (from `src/cli/ani.rs:347-396`), confirmed:

| TSV column | Rust expression in ani.rs | note |
|---|---|---|
| `ani` | `res.ani_het * 100.0` | gamma estimate, percent |
| `ani_uniform` | `res.ani * 100.0` | uniform estimate, percent |
| `af_query` / `af_reference` | direct | fractions |
| `std_err` | `res.std_err * 100.0` | percent-scaled |
| `retention`, `het_shape` | direct | het_shape can be inf |
| `ani_from_loss` / `ani_from_hist` | field `* 100.0` | percent |
| `enzyme_spread` | `res.agreement.spread * 100.0` | percent |
| `enzyme_chi2` | `res.agreement.reduced_chi2` | direct |
| `n_anchors` / `n_chains` / `n_tags` | direct / direct / `res.n_tags_in_chains` | counts |
| `synteny_blocks` / `synteny_score` / `breakpoint_count` | direct | |
| `max_block_anchors` / `mean_block_anchors` | direct | |

**Recommended path (Set A, per §4): zero code changes.** Replace
`models/gtdb_r207_linear_cal.json` with the contents of
`results/panel_by_band/linear_cal_v2_setA.json` (loaded via
`include_str!("../../models/gtdb_r207_linear_cal.json")` in
`load_embedded_model`, re-exported as `load_embedded_cal_model`, used by
`src/cli/ani.rs:266-273` when `--calibrate` is passed without
`--calibrate-model`). The current 9-feature `predict_from_result` vector
already matches the Set-A JSON exactly (same order, same scaling), and the
`embedded_model_loads` test stays at 9 features.

Exact edits **if Set B is ever shipped** (kept for reference; not
recommended per §4):

1. `src/core/calibration.rs::predict_from_result` — extend the `features`
   array from 9 to 18 entries, in this exact order (matches
   `linear_cal_v2.json::feature_names`):
   `[res.ani_het*100.0, res.ani*100.0, res.af_query, res.af_reference,
   res.std_err*100.0, res.retention, res.n_anchors as f64, res.n_chains as
   f64, res.n_tags_in_chains as f64, res.synteny_score,
   res.breakpoint_count as f64, res.agreement.spread*100.0,
   res.agreement.reduced_chi2, res.het_shape, res.ani_from_loss*100.0,
   res.ani_from_hist*100.0, res.max_block_anchors as f64,
   res.mean_block_anchors]`.
   The `predict()` imputation path already handles `het_shape = inf`
   (non-finite → `imputer_medians[i]`), matching training. The early-NaN
   guard (`!res.ani.is_finite() || res.below_detection || n_chains == 0 ||
   n_anchors == 0`) stays as-is — note it keys on `res.ani` (uniform),
   which is finite whenever `ani_het` is, so no change needed.
2. Replace `models/gtdb_r207_linear_cal.json` with the contents of
   `results/panel_by_band/linear_cal_v2.json` (loaded via
   `include_str!("../../models/gtdb_r207_linear_cal.json")` in
   `load_embedded_model`, re-exported as `load_embedded_cal_model`, used by
   `src/cli/ani.rs:266-273` when `--calibrate` is passed without
   `--calibrate-model`).
3. Update the test `embedded_model_loads` in `calibration.rs`:
   `assert_eq!(model.feature_names.len(), 9)` → `18`. The other two tests
   (`predict_is_finite_for_typical_result`, `predict_is_nan_for_failed_estimate`)
   construct full `ChainAniResult` dummies and need no changes — all 18
   fields already exist on the struct.
4. Optional but recommended: bump `name` handling — none needed, `name` is
   informational only.

After either swap, sanity-check with:
`syn2bani ani --calibrate` on a few GTDB pairs and compare the `ani_cal`
column against applying the JSON in Python (already verified equivalent:
the Rust `predict()` is the same impute→standardize→dot-product this
script's `model_predict` performs).
