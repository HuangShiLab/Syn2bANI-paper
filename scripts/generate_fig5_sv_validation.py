#!/usr/bin/env python3
"""Generate Figure 5: SV outputs agree with alignment-based truth.

Publication-quality 2x2 panel figure.
"""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import Ridge
import warnings
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS
import matplotlib.pyplot as plt

set_publication_style()
D = Path("results/gtdb50k")
OUT = Path("paper/figures/main/fig5_sv_validation")


def rank(a):
    """Return average-rank transform, handling NaN."""
    from scipy.stats import rankdata
    return pd.Series(a).rank(method='average')


def _regress_out(Z, y):
    """Stable OLS via tiny ridge to avoid singular-design warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return y - Ridge(alpha=1e-6, fit_intercept=True).fit(Z, y).predict(Z)


def robust_partial_residuals(x, y, ctrls):
    """Return (residuals_x, residuals_y, n, spearman_partial_r).

    Spearman partial correlation is computed on rank-transformed variables
    after regressing out rank-transformed controls. For visualization, the
    function returns the original-scale residuals (x and y after OLS on
    standardized controls) so that axes are interpretable.
    """
    df = pd.DataFrame({'x': x, 'y': y})
    for i, c in enumerate(ctrls):
        df[f'c{i}'] = c
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    if len(df) < 10:
        return None, None, 0, np.nan
    # Spearman partial correlation (ranks)
    x_r, y_r = rank(df['x']), rank(df['y'])
    Z_r = pd.DataFrame({f'c{i}': rank(df[f'c{i}']) for i in range(len(ctrls))}).values
    Z_r = (Z_r - Z_r.mean(axis=0)) / (Z_r.std(axis=0) + 1e-9)
    rx_r = _regress_out(Z_r, x_r)
    ry_r = _regress_out(Z_r, y_r)
    spr, _ = stats.spearmanr(rx_r, ry_r)
    # Visualization residuals (original scale, standardized controls)
    Z = pd.DataFrame({f'c{i}': df[f'c{i}'] for i in range(len(ctrls))}).values
    Z = (Z - Z.mean(axis=0)) / (Z.std(axis=0) + 1e-9)
    rx = _regress_out(Z, df['x'].values)
    ry = _regress_out(Z, df['y'].values)
    return rx, ry, len(df), spr


def main():
    # Inverted fraction data
    s2b = pd.read_csv(D / "syn2b_inverted_fraction_50k.tsv", sep="\t")
    dna = pd.read_csv(D / "dnadiff_inverted_fraction.tsv", sep="\t")
    inv = s2b.merge(dna[["pairid", "dnadiff_inverted_fraction"]], on="pairid", how="inner")
    truth = pd.read_csv(D / "truth_50k.tsv", sep="\t")
    inv = inv.merge(truth[["pairid", "anim_ani"]], on="pairid", how="inner")

    inv_hi = pd.read_csv(D / "syn2b_inverted_fraction_high_ani_all.tsv", sep="\t").merge(
        pd.read_csv(D / "dnadiff_inverted_fraction_high_ani_all.tsv", sep="\t")[["pairid", "dnadiff_inverted_fraction"]],
        on="pairid", how="inner")
    # Need ANIm truth to select strain-range subset
    hi_truth = pd.read_csv(D / "high_ani_truth.tsv", sep="\t")
    inv_hi = inv_hi.merge(hi_truth[["pairid", "anim_ani"]], on="pairid", how="inner")

    # SV comparison + truth
    sv = pd.read_csv(D / "sv_comparison_merged.tsv", sep="\t")
    sv = sv.merge(truth[["pairid", "anim_ani", "anim_af_qry"]], on="pairid", how="inner")
    sv["n_ctg"] = (sv["synteny_blocks"] - sv["breakpoint_count"]).clip(lower=1)
    sv = sv.dropna(subset=["anim_ani", "breakpoint_count", "dnadiff_breakpoints", "af_query"])

    fig, axes = plt.subplots(2, 2, figsize=figure_size(17.8, aspect=0.85))

    # (a) Inverted fraction all pairs
    ax = axes[0, 0]
    x, y = inv["syn2b_raw_inverted_fraction"], inv["dnadiff_inverted_fraction"]
    m = ~(x.isna() | y.isna())
    r, _ = stats.pearsonr(x[m], y[m])
    hb = ax.hexbin(x[m], y[m], gridsize=60, cmap="Blues", mincnt=1, linewidths=0)
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("Syn2bANI raw inverted fraction")
    ax.set_ylabel("dnadiff inverted fraction")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    label_panel(ax, "a")
    ax.set_title(f"All GTDB-R207 pairs (r = {r:.3f}, n = {m.sum()})", fontsize=8)

    # (b) High-ANI zoom (ANIm ≥ 97%)
    ax = axes[0, 1]
    sub = inv_hi[inv_hi["anim_ani"] >= 97]
    x2, y2 = sub["syn2b_raw_inverted_fraction"], sub["dnadiff_inverted_fraction"]
    mm = ~(x2.isna() | y2.isna())
    r2, _ = stats.pearsonr(x2[mm], y2[mm])
    slope, intercept, _, _, _ = stats.linregress(x2[mm], y2[mm])
    ax.scatter(x2[mm], y2[mm], alpha=0.25, s=6, c=COLORS['blue'], edgecolors='none')
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.plot([0, 1], [intercept, intercept + slope], "-", lw=1.0, c=COLORS['vermillion'])
    ax.set_xlabel("Syn2bANI raw inverted fraction")
    ax.set_ylabel("dnadiff inverted fraction")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    label_panel(ax, "b")
    ax.set_title(f"ANIm ≥ 97% (r = {r2:.3f}, slope = {slope:.3f}, n = {mm.sum()})", fontsize=8)

    # (c) Breakpoint count partial correlation
    ax = axes[1, 0]
    n_ctg_safe = sv["n_ctg"].replace([np.inf, -np.inf], np.nan)
    rx, ry, n, rp = robust_partial_residuals(
        sv["breakpoint_count"], sv["dnadiff_breakpoints"],
        [sv["anim_ani"], n_ctg_safe]
    )
    ax.scatter(rx, ry, alpha=0.25, s=5, c=COLORS['vermillion'], edgecolors='none')
    ax.axhline(0, color="k", ls="--", lw=0.8)
    ax.axvline(0, color="k", ls="--", lw=0.8)
    # symmetric axis limits clipped for visibility
    lim_x = np.percentile(np.abs(rx), 99)
    lim_y = np.percentile(np.abs(ry), 99)
    ax.set_xlim(-lim_x, lim_x)
    ax.set_ylim(-lim_y, lim_y)
    ax.set_xlabel("Breakpoint count residual")
    ax.set_ylabel("dnadiff breakpoints residual")
    label_panel(ax, "c")
    ax.set_title(f"Partial correlation (ρ = {rp:.3f}, n = {n})", fontsize=8)

    # (d) af_query vs dnadiff aligned fraction
    ax = axes[1, 1]
    x4, y4 = sv["af_query"], sv["anim_af_qry"] / 100.0
    m4 = ~(x4.isna() | y4.isna())
    r4, _ = stats.pearsonr(x4[m4], y4[m4])
    ax.scatter(x4[m4], y4[m4], alpha=0.15, s=4, c=COLORS['bluish_green'], edgecolors='none')
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("Syn2bANI af_query")
    ax.set_ylabel("dnadiff aligned fraction")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    label_panel(ax, "d")
    ax.set_title(f"Aligned fraction (r = {r4:.3f}, n = {m4.sum()})", fontsize=8)

    plt.tight_layout()
    save_figure(fig, OUT)


if __name__ == "__main__":
    main()
