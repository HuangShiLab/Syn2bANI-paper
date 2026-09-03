#!/usr/bin/env python3
"""Analyze FracMinHash scale sweep against dnadiff inverted-fraction truth.

Expected inputs (in results/gtdb50k/):
    inverted_fraction_truth_four.tsv   -- four-enzyme baseline
    inverted_fraction_truth_fmh250.tsv
    inverted_fraction_truth_fmh750.tsv
    inverted_fraction_truth_fmh1582.tsv
    inverted_fraction_truth_fmh2000.tsv
    inverted_fraction_truth_fmh6000.tsv

Output:
    results/gtdb50k/FRACMINHASH_VALIDATION.md
    results/gtdb50k/fracminhash_scale_comparison.tsv

Model:
    Var(err) = c1 * p*(1-p)/m + c0
where p = dnadiff_inverted_fraction, m = syn2b_shared_tags.
"""
import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


SCALES = {
    "four": None,
    "fmh250": 250,
    "fmh750": 750,
    "fmh1582": 1582,
    "fmh2000": 2000,
    "fmh6000": 6000,
}


def fit_model(df):
    """Fit Var(resid) = c1 * p*(1-p)/m + c0 by weighted least squares."""
    df = df.dropna(subset=["dnadiff_inverted_fraction", "syn2b_inverted_fraction", "syn2b_shared_tags"])
    p = df["dnadiff_inverted_fraction"].values
    resid = (df["syn2b_inverted_fraction"] - p).values
    m = df["syn2b_shared_tags"].values
    ok = (m > 0) & (p >= 0) & (p <= 1)
    p, resid, m = p[ok], resid[ok], m[ok]
    x = p * (1 - p) / m
    y = resid ** 2
    # simple OLS on y ~ x
    X = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    c0, c1 = beta[0], beta[1]
    pred = c0 + c1 * x
    r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - np.mean(y)) ** 2)
    return c0, c1, r2, len(df)


def summarize(df, label):
    df = df.dropna(subset=["dnadiff_inverted_fraction", "syn2b_inverted_fraction"])
    err = df["syn2b_inverted_fraction"] - df["dnadiff_inverted_fraction"]
    return {
        "scale": label,
        "n": len(df),
        "mae": np.mean(np.abs(err)),
        "rmse": np.sqrt(np.mean(err ** 2)),
        "pearson_r": np.corrcoef(df["dnadiff_inverted_fraction"], df["syn2b_inverted_fraction"])[0, 1],
        "mean_shared_tags": df["syn2b_shared_tags"].mean() if "syn2b_shared_tags" in df else np.nan,
        "median_shared_tags": df["syn2b_shared_tags"].median() if "syn2b_shared_tags" in df else np.nan,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", default="results/gtdb50k")
    args = p.parse_args()
    res = Path(args.results_dir)

    rows = []
    model_rows = []
    for label, scale in SCALES.items():
        path = res / f"inverted_fraction_truth_{label}.tsv"
        if not path.exists():
            print(f"[skip] {path} not found")
            continue
        df = pd.read_csv(path, sep="\t")
        rows.append(summarize(df, label))
        c0, c1, r2, n = fit_model(df)
        model_rows.append({
            "scale": label,
            "c0": c0,
            "c1": c1,
            "r2": r2,
            "n": n,
        })
        print(f"{label}: n={n}, MAE={rows[-1]['mae']:.4f}, c0={c0:.6f}, c1={c1:.3f}, r2={r2:.3f}")

    summary = pd.DataFrame(rows)
    models = pd.DataFrame(model_rows)
    out_tsv = res / "fracminhash_scale_comparison.tsv"
    summary.to_csv(out_tsv, sep="\t", index=False)
    print(f"wrote {out_tsv}")

    def df_to_md(df):
        lines = ["| " + " | ".join(df.columns) + " |",
                 "| " + " | ".join(["---"] * len(df.columns)) + " |"]
        for _, row in df.iterrows():
            lines.append("| " + " | ".join(str(v) for v in row) + " |")
        return "\n".join(lines)

    md = []
    md.append("# FracMinHash scale sweep validation\n")
    md.append("## Accuracy by scale\n")
    md.append(df_to_md(summary))
    md.append("\n## Error-variance model\n")
    md.append("Fitted `Var(residual) = c0 + c1 * p*(1-p)/m`, where `p` is dnadiff inverted fraction and `m` is shared_tags.\n")
    md.append(df_to_md(models))
    md.append("\n## Interpretation\n")
    md.append("- `c0` is the irreducible variance (measurement noise).\n")
    md.append("- `c1` scales the binomial sampling term; it should be stable across landmark sources if the model is correct.\n")
    md.append("- Four-enzyme and FracMinHash-750 are density-matched (~6,000 landmarks on a 4.5 Mb genome); their MAE and model constants should be most similar.\n")

    out_md = res / "FRACMINHASH_VALIDATION.md"
    out_md.write_text("\n".join(md))
    print(f"wrote {out_md}")


if __name__ == "__main__":
    main()
