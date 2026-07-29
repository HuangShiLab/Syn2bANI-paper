#!/usr/bin/env python3
"""Analyze why Syn2bANI deviates from ANIm and skani on the GTDB sample."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error


def load_data(eval_pairs, matrix):
    ep = pd.read_csv(eval_pairs, sep="\t")
    mat = pd.read_csv(matrix, sep="\t")
    mat = mat.rename(columns={
        "query": "query_asm",
        "reference": "ref_asm",
    })
    df = ep.merge(mat, on=["query_asm", "ref_asm"], how="left", suffixes=("", "_mat"))
    df["err_s2b"] = df["s2b_ani"] - df["anim_ani"]
    df["err_skani"] = df["skani_ani"] - df["anim_ani"]
    df["abs_err_s2b"] = df["err_s2b"].abs()
    df["abs_err_skani"] = df["err_skani"].abs()
    return df


def summarize_by_band(df):
    rows = []
    for band, sub in df.groupby("band"):
        rows.append({
            "band": band,
            "n": len(sub),
            "s2b_mae": sub["abs_err_s2b"].mean(),
            "s2b_bias": sub["err_s2b"].mean(),
            "skani_mae": sub["abs_err_skani"].mean(),
            "skani_bias": sub["err_skani"].mean(),
        })
    return pd.DataFrame(rows)


def summarize_by_flag(df):
    rows = []
    for flag, sub in df.groupby("s2b_flag"):
        if len(sub) < 10:
            continue
        rows.append({
            "flag": flag,
            "n": len(sub),
            "s2b_mae": sub["abs_err_s2b"].mean(),
            "s2b_bias": sub["err_s2b"].mean(),
            "skani_mae": sub["abs_err_skani"].mean(),
            "skani_bias": sub["err_skani"].mean(),
            "mean_anim": sub["anim_ani"].mean(),
        })
    return pd.DataFrame(rows)


def bin_summary(df, col, n_bins=5):
    df = df.copy()
    df["bin"] = pd.qcut(df[col], q=n_bins, duplicates="drop")
    rows = []
    for b, sub in df.groupby("bin"):
        rows.append({
            "bin": str(b),
            "n": len(sub),
            "s2b_mae": sub["abs_err_s2b"].mean(),
            "s2b_bias": sub["err_s2b"].mean(),
            "skani_mae": sub["abs_err_skani"].mean(),
            "skani_bias": sub["err_skani"].mean(),
        })
    return pd.DataFrame(rows)


def correlation_table(df, cols):
    out = []
    for c in cols:
        valid = df[[c, "err_s2b", "abs_err_s2b"]].dropna()
        if len(valid) < 30:
            continue
        out.append({
            "feature": c,
            "n": len(valid),
            "corr_err": valid[c].corr(valid["err_s2b"]),
            "corr_abs_err": valid[c].corr(valid["abs_err_s2b"]),
        })
    return pd.DataFrame(out)


def calibration_experiment(df, include_skani=True):
    """Try simple Ridge regression to predict ANIm from available features."""
    features = ["s2b_ani", "s2b_ani_uniform", "s2b_af_q", "s2b_af_r",
                "s2b_std_err", "s2b_retention", "s2b_n_anchors",
                "s2b_n_chains", "s2b_n_tags"]
    if include_skani:
        features.append("skani_ani")
    Xdf = df[features].copy()
    if include_skani:
        Xdf["skani_ani"] = Xdf["skani_ani"].fillna(Xdf["s2b_ani"])
    imputer = SimpleImputer(strategy="median")
    X_full = imputer.fit_transform(Xdf)
    y = df["anim_ani"].values
    bands = df["band"].values
    unique_bands = df["band"].unique()
    preds = np.empty(len(df))
    preds[:] = np.nan
    coef_sum = np.zeros(len(features))
    n_fits = 0
    for band in unique_bands:
        test = bands == band
        train = ~test
        if train.sum() < 50 or test.sum() < 10:
            continue
        model = RidgeCV(alphas=[1e-3, 1e-2, 0.1, 1, 10, 100], cv=5)
        model.fit(X_full[train], y[train])
        preds[test] = model.predict(X_full[test])
        coef_sum += model.coef_
        n_fits += 1
    valid = ~np.isnan(preds)
    if valid.sum() == 0:
        return None
    calib_err = preds[valid] - y[valid]
    out = {
        "n": int(valid.sum()),
        "mae": float(np.mean(np.abs(calib_err))),
        "bias": float(np.mean(calib_err)),
        "rmse": float(np.sqrt(np.mean(calib_err ** 2))),
        "features": features,
        "mean_coef": (coef_sum / n_fits).tolist() if n_fits else None,
    }
    raw_err = df.loc[valid, "s2b_ani"].values - y[valid]
    out["raw_mae"] = float(np.mean(np.abs(raw_err)))
    out["raw_bias"] = float(np.mean(raw_err))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eval-pairs", default="results/panel_by_band/eval_pairs.tsv")
    ap.add_argument("--matrix", default="results/matrix_gtdb_r207_100k_v8_final.tsv")
    ap.add_argument("--out", default="results/panel_by_band/error_driver_report.json")
    args = ap.parse_args()

    df = load_data(args.eval_pairs, args.matrix)
    print(f"loaded {len(df)} pairs")

    report = {}
    report["by_band"] = summarize_by_band(df).to_dict(orient="records")
    report["by_flag"] = summarize_by_flag(df).to_dict(orient="records")

    numeric_features = ["s2b_af_q", "s2b_af_r", "s2b_std_err", "s2b_retention",
                        "s2b_n_anchors", "s2b_n_chains", "s2b_n_tags",
                        "skani_align_frac_ref", "skani_align_frac_query"]
    report["correlations"] = correlation_table(df, numeric_features).to_dict(orient="records")

    for col in ["s2b_af_q", "s2b_retention", "s2b_n_anchors", "s2b_n_tags"]:
        report[f"by_{col}"] = bin_summary(df, col).to_dict(orient="records")

    report["calibration_ridge_band_holdout"] = calibration_experiment(df, include_skani=True)
    report["calibration_ridge_no_skani"] = calibration_experiment(df, include_skani=False)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"wrote {args.out}")

    print("\n=== by flag ===")
    print(pd.DataFrame(report["by_flag"]).to_string(index=False))
    print("\n=== correlations with s2b error ===")
    print(pd.DataFrame(report["correlations"]).to_string(index=False))
    print("\n=== calibration experiment (band holdout, with skani) ===")
    print(json.dumps(report["calibration_ridge_band_holdout"], indent=2))
    print("\n=== calibration experiment (band holdout, NO skani) ===")
    print(json.dumps(report["calibration_ridge_no_skani"], indent=2))


if __name__ == "__main__":
    sys.exit(main())
