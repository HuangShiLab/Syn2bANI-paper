#!/usr/bin/env python3
"""Validate Syn2b inverted-fraction estimates against dnadiff truth.

Reports MAE, RMSE, Pearson correlation, and checks the error model
SE = sqrt(1.504 * p * (1-p) / m + 0.0205^2) for both four-enzyme and BcgI.
"""
import argparse
import numpy as np
import pandas as pd
from scipy import stats


def rmse(a, b):
    return np.sqrt(np.mean((a - b) ** 2))


def mae(a, b):
    return np.mean(np.abs(a - b))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--truth", required=True, help="dnadiff inverted-fraction TSV")
    p.add_argument("--syn2b", required=True, help="Syn2b inverted-fraction TSV")
    p.add_argument("--label", required=True, help="label for output (four/BcgI)")
    p.add_argument("--out", required=True, help="output TSV of per-pair residuals")
    args = p.parse_args()

    truth = pd.read_csv(args.truth, sep="\t")
    syn = pd.read_csv(args.syn2b, sep="\t")
    syn = syn[syn["status"] == "ok"].copy()

    df = pd.merge(
        truth[["pairid", "dnadiff_inverted_fraction"]],
        syn[[
            "pairid",
            "syn2b_inverted_fraction",
            "syn2b_raw_inverted_fraction",
            "syn2b_shared_tags",
            "syn2b_observable_fraction",
            "syn2b_breakpoints",
        ]],
        on="pairid",
        how="inner",
    ).dropna()

    t = df["dnadiff_inverted_fraction"]
    p_est = df["syn2b_inverted_fraction"]
    p_raw = df["syn2b_raw_inverted_fraction"]
    m = df["syn2b_shared_tags"]

    c = 1.504
    sigma0 = 0.0205
    df["se_pred_min"] = np.sqrt(c * p_est * (1 - p_est) / m + sigma0**2)
    df["se_pred_raw"] = np.sqrt(c * p_raw * (1 - p_raw) / m + sigma0**2)
    df["se_pred"] = df["se_pred_raw"]
    df["resid"] = p_raw - t
    df["z"] = df["resid"] / df["se_pred"]

    print(f"\n=== {args.label} vs dnadiff truth (raw_inverted_fraction) ===")
    print(f"n = {len(df)}")
    print(f"MAE  (raw) = {mae(p_raw, t):.4f}")
    print(f"RMSE (raw) = {rmse(p_raw, t):.4f}")
    print(f"MAE  (min) = {mae(p_est, t):.4f}")
    print(f"RMSE (min) = {rmse(p_est, t):.4f}")
    print(f"mean resid = {df['resid'].mean():.4f}")
    print(f"sd resid   = {df['resid'].std():.4f}")
    r, rp = stats.pearsonr(t, p_raw)
    print(f"Pearson r (raw) = {r:.4f} (p={rp:.2e})")
    slope, intercept, rv, pv, se = stats.linregress(t, p_raw)
    print(f"p_raw = {slope:.4f} * truth + {intercept:.4f}; R2={rv**2:.4f}")
    r_min, _ = stats.pearsonr(t, p_est)
    print(f"Pearson r (min) = {r_min:.4f}")

    print("\n--- error model calibration ---")
    print(f"mean predicted SE: {df['se_pred'].mean():.4f}")
    print(f"observed RMSE:     {rmse(p_raw, t):.4f}")
    print(f"z mean: {df['z'].mean():.3f}; sd: {df['z'].std():.3f}")
    print(f"|z| > 2: {(np.abs(df['z']) > 2).mean():.3f}")
    print(f"|z| > 3: {(np.abs(df['z']) > 3).mean():.3f}")

    # tag-count dependence: bin by m and compare observed RMSE to predicted
    df["m_bin"] = pd.qcut(m, q=10, duplicates="drop")
    print("\n--- by shared-tag decile ---")
    print("m_low\tm_high\tn\tmean_m\ttruth_mean\tp_raw_mean\tobs_rmse\tpred_se")
    for interval, g in df.groupby("m_bin", observed=False):
        if len(g) < 30:
            continue
        print(
            f"{interval.left:.0f}\t{interval.right:.0f}\t{len(g)}\t"
            f"{g['syn2b_shared_tags'].mean():.0f}\t{g['dnadiff_inverted_fraction'].mean():.3f}\t"
            f"{g['syn2b_raw_inverted_fraction'].mean():.3f}\t"
            f"{rmse(g['syn2b_raw_inverted_fraction'], g['dnadiff_inverted_fraction']):.4f}\t"
            f"{g['se_pred'].mean():.4f}"
        )

    # ANI band dependence
    if "band" in truth.columns:
        df = df.merge(truth[["pairid", "band"]], on="pairid", how="left")
        print("\n--- by ANI band ---")
        for band, g in sorted(df.groupby("band"), key=lambda x: x[0]):
            if len(g) < 30:
                continue
            print(
                f"{band}: n={len(g)} truth={g['dnadiff_inverted_fraction'].mean():.3f} "
                f"p_raw={g['syn2b_raw_inverted_fraction'].mean():.3f} "
                f"rmse={rmse(g['syn2b_raw_inverted_fraction'], g['dnadiff_inverted_fraction']):.4f} "
                f"pred={g['se_pred'].mean():.4f}"
            )

    df.to_csv(args.out, sep="\t", index=False)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
