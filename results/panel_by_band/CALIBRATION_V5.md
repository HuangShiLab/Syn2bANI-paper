# Calibration v5 — v4 protocol retrained on post-rescue feature matrices

Script: `scripts/calibration_v5.py` (imports `calibration_v2`/`v3`/`v4`
helpers; identical band-holdout CV protocol as v4). Outputs:
`calibration_v5_cv.tsv`, `calibration_v5_external.tsv`, `linear_cal_v5.json`.

## Why v5

The short-contig rescue pass (Syn2bANI 0aabd0c) changed estimator features
(`af_query`, `retention`, `n_tags_in_chains`, …) on fragmented inputs, so the
deployed v4 model — trained on pre-rescue features — was inconsistent with
the shipped estimator on draft genomes. All feature matrices were re-run on
HPC with the 068119c binary (`syn2bani ani --verbose -t 1`, default 4-enzyme
panel) and v5 was trained on the post-rescue matrices. Note the v4→v5 binary
span (98177dc → 068119c) also includes the `breakpoint_count` formula fix
(22b853b) and the IUPAC-aware digest fix (f054dbb), so the feature diff is
broader than the rescue pass alone (see "Feature diff" below).

## Training set

- Old: `results/anim_truth_2074_v9rescue.tsv` joined to `eval_pairs.tsv`
  → 2,053 rows (21 dropped, non-finite `ani_gated`, gate = `none` — the same
  21 pairs as v3/v4).
- New: `results/anim_truth_hi95_v9rescue.tsv` ⋈ `results/anim_truth_hi95.tsv`
  (accession-keyed) → 467 rows, 0 dropped, 0 gate fallbacks, all in
  [95.0, 99.5) → band `0.95-0.99`.
- **Merged n = 2,520**; top-band n: 72 → 539 (same composition as v4).
- QA: the shipped `gate` column matches the verified Python gate rule on
  2,053/2,053 and 467/467 rows (asserted in the script).

## Protocol

Identical to v4: band-holdout CV (train 3 bands, test held-out, rotate),
per-fold median impute + StandardScaler + RidgeCV, training rows shuffled
before CV (RidgeCV inner KFold is order-sensitive), seeds 0–4, seed 0
primary. Variants v5a (9 gated features) and v5b (v5a + `gate_fallback`);
selection rule: v5b only if overall CV MAE < 0.995·v5a. v4a/v4b are rerun
on the old 2,053-row subset of the post-rescue matrix as the reference line.

Seed sensitivity (overall CV MAE, seeds 0–4):

| experiment | seed MAEs | mean | spread |
|---|---|---|---|
| v4a_ref (n=2,053, post-rescue) | 0.7716 / 0.7716 / 0.7729 / 0.7793 / 0.7716 | 0.7734 | 0.0077 |
| v4b_ref (n=2,053, post-rescue) | 0.7712 / 0.7712 / 0.7725 / 0.7806 / 0.7712 | 0.7733 | 0.0094 |
| v5a (n=2,520) | 0.7313 / 0.7313 / 0.7327 / 0.7313 / 0.7327 | 0.7319 | 0.0014 |
| v5b (n=2,520) | 0.7313 / 0.7313 / 0.7328 / 0.7313 / 0.7328 | 0.7319 | 0.0015 |

## Variant selection

v5a = 0.7313, v5b = 0.7313 (seed 0). Rule (v5b only if < 0.995·v5a =
0.7276): **not met on any of the 5 seeds → v5a selected** (9 features, no
gate indicator — same outcome as v4). Final fit: in-sample MAE 0.6414,
**alpha = 10.0**.

## Per-band CV MAE (seed 0)

| band | n | raw_gated | v4a (pre-rescue) | v4a_ref (post-rescue) | v5a |
|---|---|---|---|---|---|
| 0.8-0.85 | 463 | 2.9520 | 0.9457 | 0.8059 | 0.8218 |
| 0.85-0.9 | 829 | 2.8083 | 0.9117 | 0.7315 | 0.7657 |
| 0.9-0.95 | 689 | 1.4677 | 0.7103 | 0.7639 | 0.6441 |
| 0.95-0.99 | 539 | 0.6664 | 0.8621 | 1.0859 (n=72) | **0.7119** |
| all | 2,520 | 2.0100 | 0.8523 | 0.7716 (n=2,053) | **0.7313** |

Post-rescue features are substantially more informative: the same v4a
protocol on the post-rescue 2,053-row subset drops overall CV MAE
0.9371 → 0.7716 (v3a_ref pre-rescue vs v4a_ref post-rescue), and v5a on
the merged set improves 0.8523 → 0.7313 vs v4a. Top-band MAE
0.8621 → 0.7119; raw gated MAE 2.4387 → 2.0100 (merged) and
2.8194 → 2.3231 (2,074-pair set).

## External validation (post-rescue feature re-runs)

Both external sets were re-run with the 068119c binary:
`results/oral_gut_1225_v9rescue.tsv` (oral/gut, vs FastANI) and
`results/gating_flag/midani_15_v9rescue.tsv` (mid-ANI, vs ANIm). No
version skew: v5 external rows use post-rescue features throughout.

| dataset | method | n | MAE | bias | r |
|---|---|---|---|---|---|
| oral/gut vs FastANI | raw_gated | 100 | 0.5505 | +0.4200 | 0.9804 |
| oral/gut vs FastANI | calibrated_v5 | 100 | 0.4965 | +0.4444 | 0.9947 |
| oral/gut vs FastANI | (v4 was, pre-rescue feats) | 100 | 0.4634 | +0.3744 | 0.9949 |
| mid-ANI vs ANIm | raw_gated | 15 | 0.9588 | +0.8902 | 0.8734 |
| mid-ANI vs ANIm | calibrated_v5 | 15 | 0.8061 | −0.4262 | 0.6938 |
| mid-ANI vs ANIm | (v4 was) | 15 | 0.7942 | −0.2192 | 0.5467 |

The mid-ANI raw features are bit-identical pre/post rescue (the 15 pairs
involve genomes unaffected by the estimator changes); only the calibrated
value moved. The oral/gut raw MAE barely moved (0.5518 → 0.5505).

## Gates vs v4 — **FAIL (3 of 5)**

Per-band gates compare v5a against the v4a protocol rerun on the same
post-rescue 2,053-row subset; external gates compare against v4's external
MAEs (computed on pre-rescue features, the closest available reference).

| gate | reference | v5 | verdict |
|---|---|---|---|
| band 0.8-0.85 Δ ≤ +0.02 | 0.8059 (v4a_ref) | 0.8218 | PASS (+0.0159) |
| band 0.85-0.9 Δ ≤ +0.02 | 0.7315 (v4a_ref) | 0.7657 | **FAIL (+0.0342)** |
| band 0.9-0.95 Δ ≤ +0.02 | 0.7639 (v4a_ref) | 0.6441 | PASS (−0.1198) |
| oral/gut external ≤ 0.4634 | 0.4634 (v4) | 0.4965 | **FAIL (+0.0331)** |
| mid-ANI external ≤ 0.7942 | 0.7942 (v4) | 0.8061 | **FAIL (+0.0119)** |

Same pattern as v4 vs v3 (which failed 3 of 5): a large overall/top-band
improvement traded for a modest 0.85-0.9 regression and slightly worse
external MAEs. Numbers delivered as-is; no tuning was done to pass the
gates.

## Feature diff (v4-era 98177dc binary → 068119c)

| set | pairs | any training-feature change | \|Δani_gated\| > 0.01 | > 0.1 | > 1.0 | gate flips |
|---|---|---|---|---|---|---|
| 2074 GTDB | 2,074 | 1,607 | 1,497 | 1,136 | 403 | 279 |
| 467 hi95 | 467 | 315 | 248 | 141 | 18 | 17 |
| combined | 2,541 | 1,922 | 1,745 | 1,277 | 421 | 296 |

Changed-pair \|Δani_gated\| distribution: 2074 set median 0.342, p90 2.037,
max 9.675; hi95 set median 0.072, p90 0.682, max 2.285. No pair entered or
left the non-finite `ani_gated` set. The diff is much broader than the
rescue pass alone because the binary span also contains the IUPAC-aware
digest fix (f054dbb) and the `breakpoint_count` fix (22b853b).

## Final model (`linear_cal_v5.json`, `gtdb_r207_linear_cal_v5`)

- Variant v5a, 9 features (Rust order): `ani_gated, ani_uniform, af_query,
  af_reference, std_err, retention, n_anchors, n_chains, n_tags_in_chains`
- Coefficients: [0.5364, −0.0344, 0.5907, 0.6586, −0.4012, 2.6278, 0.1881,
  0.0853, −0.4595]; intercept 90.5256
- alpha = 10.0; training_n = 2,520; training (in-sample) MAE = 0.6414
- Final fit rows shuffled with seed 0 (RidgeCV inner-KFold order
  sensitivity).
- Not deployed: `models/gtdb_r207_linear_cal.json` in the code repo still
  holds v4; deployment is a separate reviewed step.
