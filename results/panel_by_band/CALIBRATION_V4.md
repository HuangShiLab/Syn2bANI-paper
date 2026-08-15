# Calibration v4 — v3 protocol + 467 new 95–99.5% ANIm pairs

Script: `scripts/calibration_v4.py` (imports `calibration_v2`/`calibration_v3`
helpers; same band-holdout CV protocol as v3). Outputs:
`calibration_v4_cv.tsv`, `calibration_v4_external.tsv`, `linear_cal_v4.json`.

## Training set

- Old: v3 gated matrix joined to `eval_pairs.tsv` → 2,053 rows (21 dropped,
  non-finite `ani_gated`; identical to v3).
- New: `results/anim_truth_hi95_gated.tsv` ⋈ `results/anim_truth_hi95.tsv`
  (accession-keyed, no acc2seqid needed) → 467 rows, 0 dropped, 0 gate
  fallbacks. All 467 have anim_ani in [95.0, 99.5) (range 95.00–99.15), so
  all land in the single matching `eval_pairs` band **`0.95-0.99`**.
- **Merged n = 2,520**; top-band n: **72 → 539**.

## Protocol deviation from v3 (deliberate)

RidgeCV's default KFold is order-sensitive, so training rows are shuffled
with a fixed seed before CV. Seeds 0–4 reported; seed 0 is the primary
number. v3a/v3b were rerun on the old 2,053-row subset with the same seeded
protocol as the reference line (matches the original v3 run within ~0.005
MAE, confirming the order sensitivity is real but small on this data).

Seed sensitivity (overall CV MAE, seeds 0–4):

| experiment | seed MAEs | mean | spread |
|---|---|---|---|
| v3a_ref (n=2,053) | 0.9371 / 0.9363 / 0.9324 / 0.9363 / 0.9376 | 0.9359 | 0.0052 |
| v3b_ref (n=2,053) | 0.9338 / 0.9328 / 0.9287 / 0.9328 / 0.9353 | 0.9327 | 0.0066 |
| v4a (n=2,520) | 0.8523 / 0.8515 / 0.8515 / 0.8515 / 0.8515 | 0.8517 | 0.0008 |
| v4b (n=2,520) | 0.8501 / 0.8492 / 0.8492 / 0.8506 / 0.8492 | 0.8497 | 0.0014 |

## Variant selection

v4a = 0.8523, v4b = 0.8501 (seed 0). Rule (v4b only if < 0.995·v4a =
0.8480): **not met on any of the 5 seeds → v4a selected** (9 features, no
gate indicator). Final fit: in-sample MAE 0.7254, **alpha = 1.0**.

## Per-band CV MAE (seed 0), v3a reference vs v4a

| band | n (v3→v4) | raw_gated | v3a_ref | v4a | Δ |
|---|---|---|---|---|---|
| 0.8-0.85 | 463 | 3.3686 | 0.9705 | 0.9457 | −0.0248 |
| 0.85-0.9 | 829 | 3.4897 | 0.8686 | 0.9117 | **+0.0431** |
| 0.9-0.95 | 689 | 1.8175 | 0.9564 | 0.7103 | −0.2461 |
| 0.95-0.99 | 72→539 | 1.1589→0.8174 | 1.3260 | **0.8621** | −0.4639 |
| all | 2,053→2,520 | 2.4387 | 0.9371 | 0.8523 | −0.0848 |

The key question: top-band (95–99) MAE dropped **1.326 → 0.862** with n
72 → 539, and the bias shrank from −1.12 to −0.60 (r 0.15 → 0.42). The old
72-row band MAE is identical across seeds (held-out model never retrains on
it); the merged-set numbers are stable to ±0.001.

## External validation (identical to v3)

| dataset | method | n | MAE | bias | r |
|---|---|---|---|---|---|
| oral/gut vs FastANI | raw_gated | 100 | 0.5518 | +0.4213 | 0.9804 |
| oral/gut vs FastANI | calibrated_v4 | 100 | 0.4634 | +0.3744 | 0.9949 |
| oral/gut vs FastANI | (v3 was) | 100 | 0.4243 | +0.1508 | 0.9948 |
| mid-ANI vs ANIm | raw_gated | 15 | 0.9588 | +0.8902 | 0.8734 |
| mid-ANI vs ANIm | calibrated_v4 | 15 | 0.7942 | −0.2192 | 0.5467 |
| mid-ANI vs ANIm | (v3 was) | 15 | 0.7391 | −0.0435 | 0.5668 |

## Gates vs v3 — **FAIL (3 of 5)**

| gate | v3 | v4 | verdict |
|---|---|---|---|
| band 0.8-0.85 Δ ≤ +0.02 | 0.9705 | 0.9457 | PASS (−0.0248) |
| band 0.85-0.9 Δ ≤ +0.02 | 0.8686 | 0.9117 | **FAIL (+0.0431)** |
| band 0.9-0.95 Δ ≤ +0.02 | 0.9564 | 0.7103 | PASS (−0.2461) |
| oral/gut external ≤ 0.4243 | 0.4243 | 0.4634 | **FAIL (+0.0391)** |
| mid-ANI external ≤ 0.7391 | 0.7391 | 0.7942 | **FAIL (+0.0551)** |

v4 trades a large top-band/0.9-0.95 improvement and better overall CV MAE
for a modest 0.85-0.9 regression and slightly worse external MAEs. Numbers
delivered as-is; no tuning was done to pass the gates.

## Final model (`linear_cal_v4.json`, `gtdb_r207_linear_cal_v4`)

- Variant v4a, 9 features (Rust order): `ani_gated, ani_uniform, af_query,
  af_reference, std_err, retention, n_anchors, n_chains, n_tags_in_chains`
- Coefficients: [0.8729, −1.1472, 0.6772, 0.8535, −0.4475, 3.3199, 0.2257,
  0.1351, −0.5221]; intercept 90.5256
- alpha = 1.0; training_n = 2,520; training (in-sample) MAE = 0.7254
- Final fit rows shuffled with seed 0 (RidgeCV inner-KFold order sensitivity).
