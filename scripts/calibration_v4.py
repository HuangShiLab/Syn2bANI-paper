#!/usr/bin/env python3
"""Calibration v4: v3 protocol + 467 new 95-99.5% ANIm band pairs.

Training set = v3 gated matrix (results/anim_truth_2074_gated.tsv joined to
eval_pairs.tsv; n = 2,053 after dropping 21 non-finite ani_gated) plus the
new hi95 pairs (results/anim_truth_hi95_gated.tsv joined to
results/anim_truth_hi95.tsv; accession-keyed, n = 467). All 467 new pairs
have anim_ani in [95, 99.5) and are assigned the matching eval_pairs band
label "0.95-0.99" (band cutoffs verified: that band holds the >=95 pairs).
Merged n = 2,520; top-band n goes 72 -> 539.

Protocol identical to v3 (scripts/calibration_v3.py): band-holdout CV
(train 3 bands, test held-out, rotate), per-fold median impute +
StandardScaler + RidgeCV. One fix: RidgeCV's default KFold is
order-sensitive, so training rows are shuffled with a fixed seed before CV
and the final overall number is reported across 5 seeds (0..4); seed 0 is
the primary/reporting seed.

Variants:
  v4a: V3A features (ani_gated, ani_uniform, af_query, af_reference,
       std_err, retention, n_anchors, n_chains, n_tags)
  v4b: v4a + gate_fallback (0/1 indicator)
Selection rule (same as v3): v4b only if overall CV MAE < 0.995 * v4a.
v3a/v3b are rerun on the old 2,053-row subset with the same seeded
protocol as the reference line.

External validation identical to v3 (oral/gut 100 same-species pairs vs
FastANI with apply_gate_rule recomputation; mid-ANI 15 pairs vs ANIm).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import calibration_v2 as c2  # protocol helpers (band_holdout_cv, metric_rows)
import calibration_v3 as c3  # v3 loaders / external validation / feature sets

ROOT = Path(__file__).resolve().parent.parent
HI95_GATED = ROOT / "results/anim_truth_hi95_gated.tsv"
HI95_TRUTH = ROOT / "results/anim_truth_hi95.tsv"
OUT_CV = ROOT / "results/panel_by_band/calibration_v4_cv.tsv"
OUT_EXT = ROOT / "results/panel_by_band/calibration_v4_external.tsv"
OUT_JSON = ROOT / "results/panel_by_band/linear_cal_v4.json"

SEEDS = [0, 1, 2, 3, 4]
PRIMARY_SEED = 0
TOP_BAND = "0.95-0.99"
HI95_LO, HI95_HI = 95.0, 99.5  # ANIm band cutoffs for the new pairs


def load_hi95():
    """New 95-99.5% band pairs, accession-keyed (no acc2seqid needed)."""
    mat = pd.read_csv(HI95_GATED, sep="\t")
    truth = pd.read_csv(HI95_TRUTH, sep="\t")
    df = mat.merge(truth[["query", "reference", "anim_ani"]],
                   on=["query", "reference"], how="inner")
    assert len(df) == 467, len(df)
    in_band = df["anim_ani"].between(HI95_LO, HI95_HI, inclusive="left")
    print(f"hi95 pairs in [{HI95_LO}, {HI95_HI}): {int(in_band.sum())}/467 "
          f"(anim_ani range {df['anim_ani'].min():.2f}-"
          f"{df['anim_ani'].max():.2f})")
    assert in_band.all(), "hi95 pairs outside the 95-99.5 band"
    df["band"] = TOP_BAND
    df = df.replace([np.inf, -np.inf], np.nan)
    df["gate_fallback"] = (df["gate"] == "uniform_fallback").astype(float)
    n_drop = int(df["ani_gated"].isna().sum())
    df = df[df["ani_gated"].notna()].reset_index(drop=True)
    print(f"hi95 training rows: {len(df)} (dropped {n_drop} non-finite; "
          f"fallbacks: {int(df['gate_fallback'].sum())})")
    return df


def shuffled(df, seed):
    return df.sample(frac=1.0, random_state=seed).reset_index(drop=True)


def cv_experiment(df, name, features, seeds):
    """Band-holdout CV per seed; returns metric rows with a seed column."""
    rows = []
    for seed in seeds:
        d = shuffled(df, seed)
        preds = c2.band_holdout_cv(d, features, "ridge")
        for r in c2.metric_rows(name, preds, d["anim_ani"].values,
                                d["band"].values):
            r["seed"] = seed
            rows.append(r)
    return rows


def add_mean_rows(cv):
    """Append seed='mean' rows (mean of per-seed metrics) per exp/band."""
    mean = (cv.groupby(["experiment", "band"], as_index=False)
              [["n", "MAE", "RMSE", "bias", "r"]].mean())
    mean["seed"] = "mean"
    return pd.concat([cv, mean[cv.columns]], ignore_index=True)


def overall(cv, name, seed=PRIMARY_SEED):
    m = cv[(cv.experiment == name) & (cv.band == "all") & (cv.seed == seed)]
    return m["MAE"].iloc[0]


def train_final(df, features, rust_names, seed):
    """Final fit on all rows, shuffled with a fixed seed (RidgeCV's inner
    KFold is order-sensitive)."""
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler
    d = shuffled(df, seed)
    X = d[features].values
    y = d["anim_ani"].values
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    Xs = scaler.fit_transform(imputer.fit_transform(X))
    model = RidgeCV(alphas=c2.ALPHAS + [1000], cv=5)
    model.fit(Xs, y)
    mae = float(np.abs(model.predict(Xs) - y).mean())
    cal = {
        "name": "gtdb_r207_linear_cal_v4",
        "feature_names": rust_names,
        "means": scaler.mean_.tolist(),
        "scales": scaler.scale_.tolist(),
        "coefficients": model.coef_.tolist(),
        "intercept": float(model.intercept_),
        "imputer_medians": imputer.statistics_.tolist(),
        "training_n": int(len(d)),
        "training_mae": mae,
    }
    return cal, ("ridge", imputer, scaler, model), mae, float(model.alpha_)


def main():
    old = c3.load_gated()          # 2,053 rows, bands from eval_pairs
    new = load_hi95()              # 467 rows, band = 0.95-0.99
    df = pd.concat([old, new], ignore_index=True)
    print(f"merged training rows: {len(df)} "
          f"(top band n: {int((old.band == TOP_BAND).sum())} -> "
          f"{int((df.band == TOP_BAND).sum())})")

    rows = []
    # raw reference on the merged set
    d0 = shuffled(df, PRIMARY_SEED)
    for r in c2.metric_rows("raw_gated", d0["ani_gated"].values,
                            d0["anim_ani"].values, d0["band"].values):
        r["seed"] = PRIMARY_SEED
        rows.append(r)
    # v3 reference rerun on the old 2,053-row subset, same seeded protocol
    rows += cv_experiment(old, "v3a_ref_gated9", c3.V3A, SEEDS)
    rows += cv_experiment(old, "v3b_ref_gated9+gate", c3.V3B, SEEDS)
    # v4 variants on the merged 2,520-row set
    rows += cv_experiment(df, "v4a_gated9", c3.V3A, SEEDS)
    rows += cv_experiment(df, "v4b_gated9+gate", c3.V3B, SEEDS)

    cv = pd.DataFrame(rows)
    cv = add_mean_rows(cv).round(4)
    cv = cv[["experiment", "seed", "band", "n", "MAE", "RMSE", "bias", "r"]]
    cv.to_csv(OUT_CV, sep="\t", index=False)
    pd.set_option("display.width", 160)
    print(cv[cv.band == "all"].to_string(index=False))

    print("\n=== seed sensitivity (overall CV MAE) ===")
    for name in ["v3a_ref_gated9", "v3b_ref_gated9+gate", "v4a_gated9",
                 "v4b_gated9+gate"]:
        per_seed = [overall(cv, name, s) for s in SEEDS]
        print(f"{name:22s} " +
              " ".join(f"{v:.4f}" for v in per_seed) +
              f"  mean={np.mean(per_seed):.4f} spread="
              f"{max(per_seed) - min(per_seed):.4f}")

    mae_a, mae_b = overall(cv, "v4a_gated9"), overall(cv, "v4b_gated9+gate")
    use_b = mae_b < 0.995 * mae_a  # gate indicator must earn its parameter
    # check the decision is not an artifact of the primary seed
    decisions = {s: overall(cv, "v4b_gated9+gate", s)
                 < 0.995 * overall(cv, "v4a_gated9", s) for s in SEEDS}
    feats = c3.V3B if use_b else c3.V3A
    names = c3.RUST_NAMES_V3B if use_b else c3.RUST_NAMES_V3B[:9]
    tag = "v4b" if use_b else "v4a"
    print(f"\nv4a={mae_a:.4f} v4b={mae_b:.4f} (seed {PRIMARY_SEED}) "
          f"-> final: {tag}; use_b per seed: {decisions}")

    cal, fitted, mae_in, alpha = train_final(df, feats, names, PRIMARY_SEED)
    with open(OUT_JSON, "w") as fh:
        json.dump(cal, fh, indent=2)
    print(f"final v4 in-sample MAE {mae_in:.4f}, alpha {alpha}; "
          f"wrote {OUT_JSON}")

    ext, og = c3.external_oralgut(fitted, feats, "calibrated_v4")
    mrows, md = c3.external_midani(fitted, feats, "calibrated_v4")
    ext += mrows
    ext = pd.DataFrame(ext).round(4)
    ext.to_csv(OUT_EXT, sep="\t", index=False)
    print("\n=== external validation (v4) ===")
    print(ext.to_string(index=False))
    print(f"\noral/gut fallbacks recomputed: "
          f"{int(og['gate_fallback'].sum())}/100; mid-ANI fallbacks: "
          f"{int(md['gate_fallback'].sum())}/15")

    # gates vs v3 (primary seed, per-band v4-selected vs v3a reference)
    print("\n=== gates vs v3 ===")
    v3_ref = "v3a_ref_gated9"
    v4_sel = "v4b_gated9+gate" if use_b else "v4a_gated9"
    ok = True
    for band in c2.BAND_ORDER:
        if band == TOP_BAND:
            continue
        m3 = cv[(cv.experiment == v3_ref) & (cv.band == band) &
                (cv.seed == PRIMARY_SEED)]["MAE"].iloc[0]
        m4 = cv[(cv.experiment == v4_sel) & (cv.band == band) &
                (cv.seed == PRIMARY_SEED)]["MAE"].iloc[0]
        passed = m4 <= m3 + 0.02
        ok &= passed
        print(f"band {band}: v3 {m3:.4f} -> v4 {m4:.4f} "
              f"(delta {m4 - m3:+.4f}) {'PASS' if passed else 'FAIL'}")
    og4 = ext[(ext.dataset == "oralgut_same_species_vs_fastani_gated") &
              (ext.method == "calibrated_v4")]["MAE"].iloc[0]
    mid4 = ext[(ext.dataset == "midani_15_vs_anim_gated") &
               (ext.method == "calibrated_v4")]["MAE"].iloc[0]
    og_ok = og4 <= 0.4243
    mid_ok = mid4 <= 0.7391
    print(f"oral/gut external: v3 0.4243 -> v4 {og4:.4f} "
          f"{'PASS' if og_ok else 'FAIL'}")
    print(f"mid-ANI external:  v3 0.7391 -> v4 {mid4:.4f} "
          f"{'PASS' if mid_ok else 'FAIL'}")
    ok &= og_ok and mid_ok
    print(f"ALL GATES {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
