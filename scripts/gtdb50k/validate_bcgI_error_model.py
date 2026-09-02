#!/usr/bin/env python3
"""Validate the single-enzyme error model against the four-enzyme panel.

The model predicts:
    SE = sqrt(1.504 * p * (1 - p) / m + 0.0205^2)
where p is the inverted fraction and m is the number of shared tags.
Dropping from four enzymes to one should reduce m by ~4x and raise the SE
by the amount the formula predicts.
"""
import argparse
import numpy as np
import pandas as pd
from scipy import stats


def load(path):
    df = pd.read_csv(path, sep="\t")
    df = df[df["status"] == "ok"].copy()
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--four", required=True, help="four-enzyme inverted-fraction TSV")
    p.add_argument("--bcgi", required=True, help="BcgI inverted-fraction TSV")
    p.add_argument("--out", required=True, help="output validation TSV")
    p.add_argument("--plot", help="optional output PNG for residual diagnostic")
    args = p.parse_args()

    df4 = load(args.four)
    df1 = load(args.bcgi)

    # keep columns needed for the model
    keep = [
        "pairid",
        "syn2b_inverted_fraction",
        "syn2b_raw_inverted_fraction",
        "syn2b_shared_tags",
        "syn2b_observable_fraction",
    ]
    df4 = df4[keep].rename(
        columns={
            "syn2b_inverted_fraction": "p4",
            "syn2b_raw_inverted_fraction": "p4_raw",
            "syn2b_shared_tags": "m4",
            "syn2b_observable_fraction": "obs4",
        }
    )
    df1 = df1[keep].rename(
        columns={
            "syn2b_inverted_fraction": "p1",
            "syn2b_raw_inverted_fraction": "p1_raw",
            "syn2b_shared_tags": "m1",
            "syn2b_observable_fraction": "obs1",
        }
    )

    df = pd.merge(df4, df1, on="pairid", how="inner")
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    n = len(df)
    print(f"pairs with both four-enzyme and BcgI results: {n}")

    c = 1.504
    sigma0 = 0.0205

    # predicted SE for each measurement
    df["se4"] = np.sqrt(c * df["p4_raw"] * (1 - df["p4_raw"]) / df["m4"] + sigma0**2)
    df["se1"] = np.sqrt(c * df["p1_raw"] * (1 - df["p1_raw"]) / df["m1"] + sigma0**2)

    # expected SE of the difference if measurements are independent
    df["se_diff_pred"] = np.sqrt(df["se4"] ** 2 + df["se1"] ** 2)
    df["diff"] = df["p1_raw"] - df["p4_raw"]
    df["diff_min"] = df["p1"] - df["p4"]
    df["z"] = df["diff"] / df["se_diff_pred"]

    # tag count ratio
    df["m_ratio"] = df["m1"] / df["m4"]

    print("\n=== tag counts ===")
    print(f"m4 median: {df['m4'].median():.1f}; mean: {df['m4'].mean():.1f}")
    print(f"m1 median: {df['m1'].median():.1f}; mean: {df['m1'].mean():.1f}")
    print(f"m1/m4 median: {df['m_ratio'].median():.3f}; mean: {df['m_ratio'].mean():.3f}")

    print("\n=== inverted fractions (raw) ===")
    print(f"p4_raw mean: {df['p4_raw'].mean():.4f}; sd: {df['p4_raw'].std():.4f}")
    print(f"p1_raw mean: {df['p1_raw'].mean():.4f}; sd: {df['p1_raw'].std():.4f}")
    print(f"diff (raw) mean: {df['diff'].mean():.4f}; sd: {df['diff'].std():.4f}")
    print(f"diff (min) mean: {df['diff_min'].mean():.4f}; sd: {df['diff_min'].std():.4f}")

    print("\n=== predicted vs observed SE of difference ===")
    print(f"mean predicted SE(diff): {df['se_diff_pred'].mean():.4f}")
    print(f"observed SD(diff):       {df['diff'].std():.4f}")

    # robust scale estimate (MAD * 1.4826)
    mad = np.median(np.abs(df["diff"] - np.median(df["diff"]))) * 1.4826
    print(f"observed MAD*1.4826(diff): {mad:.4f}")

    # z-score distribution under the model
    print("\n=== z = diff / predicted_SE ===")
    print(f"mean z: {df['z'].mean():.3f}; sd: {df['z'].std():.3f}")
    print(f"|z| > 2 fraction: {(np.abs(df['z']) > 2).mean():.3f}")
    print(f"|z| > 3 fraction: {(np.abs(df['z']) > 3).mean():.3f}")

    # binned comparison by m4 deciles or p4 bands
    df["p4_band"] = pd.cut(df["p4_raw"], bins=[0, 0.05, 0.2, 0.5, 1.0], include_lowest=True)
    print("\n=== by inverted-fraction band (raw) ===")
    for band, g in df.groupby("p4_band", observed=False):
        if len(g) < 50:
            continue
        print(
            f"{band}: n={len(g):5d} "
            f"p4_raw={g['p4_raw'].mean():.3f} p1_raw={g['p1_raw'].mean():.3f} "
            f"sd(diff)={g['diff'].std():.4f} mean_pred_se={g['se_diff_pred'].mean():.4f} "
            f"m1/m4={g['m_ratio'].mean():.3f}"
        )

    # correlation of the two estimates
    r, rp = stats.pearsonr(df["p4_raw"], df["p1_raw"])
    print(f"\nPearson r(p4_raw, p1_raw) = {r:.4f} (p={rp:.2e})")

    # regression of p1 on p4
    slope, intercept, rv, pv, se = stats.linregress(df["p4_raw"], df["p1_raw"])
    print(f"p1_raw = {slope:.4f} * p4_raw + {intercept:.4f}; R2={rv**2:.4f}")

    # save detailed table
    df.to_csv(args.out, sep="\t", index=False)
    print(f"\nwrote {args.out}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(14, 4))

        ax = axes[0]
        ax.hexbin(df["p4_raw"], df["p1_raw"], gridsize=80, mincnt=1, cmap="YlOrRd")
        ax.plot([0, 1], [0, 1], "k--", lw=1)
        ax.set_xlabel("four-enzyme raw inverted fraction")
        ax.set_ylabel("BcgI raw inverted fraction")
        ax.set_title(f"r={r:.3f}, n={n}")

        ax = axes[1]
        ax.hist(df["z"], bins=100, range=(-5, 5), density=True, alpha=0.7)
        xs = np.linspace(-5, 5, 200)
        ax.plot(xs, stats.norm.pdf(xs, 0, 1), "k--", lw=1.5)
        ax.set_xlabel("z = (p1 - p4) / predicted SE")
        ax.set_ylabel("density")
        ax.set_title(f"mean z={df['z'].mean():.2f}, sd={df['z'].std():.2f}")

        ax = axes[2]
        ax.scatter(df["se_diff_pred"], np.abs(df["diff"]), alpha=0.1, s=5)
        ax.plot([0, 0.5], [0, 0.5], "k--", lw=1)
        ax.set_xlabel("predicted SE of difference")
        ax.set_ylabel("|observed difference|")
        ax.set_title("predicted vs observed spread")

        fig.tight_layout()
        fig.savefig(args.plot, dpi=300)
        print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
