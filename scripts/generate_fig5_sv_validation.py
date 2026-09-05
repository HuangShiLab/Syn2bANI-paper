#!/usr/bin/env python3
"""Generate Figure 5: SV outputs agree with alignment-based truth.

Publication-quality 3x2 panel figure.  Top two rows show accuracy
metrics against dnadiff; bottom row shows that this accuracy is
achieved with orders-of-magnitude lower run time.
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
E = Path("results/efficiency_v8")
OUT = Path("paper/figures/main/fig5_sv_validation")


def rank(a):
    """Return average-rank transform, handling NaN."""
    from scipy.stats import rankdata
    return pd.Series(a).rank(method="average")


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
    df = pd.DataFrame({"x": x, "y": y})
    for i, c in enumerate(ctrls):
        df[f"c{i}"] = c
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    if len(df) < 10:
        return None, None, 0, np.nan
    # Spearman partial correlation (ranks)
    x_r, y_r = rank(df["x"]), rank(df["y"])
    Z_r = pd.DataFrame({f"c{i}": rank(df[f"c{i}"]) for i in range(len(ctrls))}).values
    Z_r = (Z_r - Z_r.mean(axis=0)) / (Z_r.std(axis=0) + 1e-9)
    rx_r = _regress_out(Z_r, x_r)
    ry_r = _regress_out(Z_r, y_r)
    spr, _ = stats.spearmanr(rx_r, ry_r)
    # Visualization residuals (original scale, standardized controls)
    Z = pd.DataFrame({f"c{i}": df[f"c{i}"] for i in range(len(ctrls))}).values
    Z = (Z - Z.mean(axis=0)) / (Z.std(axis=0) + 1e-9)
    rx = _regress_out(Z, df["x"].values)
    ry = _regress_out(Z, df["y"].values)
    return rx, ry, len(df), spr


def load_sv_efficiency():
    """Return aggregated SV timing data for panels e-f."""
    sv = pd.read_csv(E / "sv_benchmark.tsv", sep="\t")
    sv_summary = (
        sv.groupby(["mode", "n_genomes", "n_pairs"])
        .agg(skani_wall_s=("skani_wall_s", "median"),
             dnadiff_wall_s=("dnadiff_wall_s", "median"))
        .reset_index()
    )
    sv_summary["total_wall_s"] = sv_summary["skani_wall_s"] + sv_summary["dnadiff_wall_s"]

    st = pd.read_csv(E / "syn2b_struct_benchmark.tsv", sep="\t")
    st_summary = (
        st.groupby(["mode", "n_genomes", "n_pairs"])
        .agg(struct_wall_s=("struct_wall_s", "median"))
        .reset_index()
    )

    rt = pd.read_csv(E / "runtime_scaling.tsv", sep="\t")
    sk = pd.read_csv(E / "sketch_benchmark.tsv", sep="\t")

    # syn2bANI sketch-reuse workflow: sketch + ani from sketches + struct
    ani_sk = rt[(rt["tool"] == "syn2bani") & (rt["mode"] == "ani_sketches")][["n_genomes", "n_pairs", "wall_s"]]
    sk_syn = sk[sk["tool"] == "syn2bani"][["n_genomes", "wall_s"]].rename(columns={"wall_s": "sketch_s"})
    s2b_sk = (st_summary.merge(ani_sk, on=["n_genomes", "n_pairs"], how="inner")
              .merge(sk_syn, on="n_genomes", how="inner").copy())
    s2b_sk["total_wall_s"] = s2b_sk["struct_wall_s"] + s2b_sk["wall_s"] + s2b_sk["sketch_s"]

    return sv_summary, s2b_sk


def panel_a(ax):
    # Inverted fraction data
    s2b = pd.read_csv(D / "syn2b_inverted_fraction_50k.tsv", sep="\t")
    dna = pd.read_csv(D / "dnadiff_inverted_fraction.tsv", sep="\t")
    inv = s2b.merge(dna[["pairid", "dnadiff_inverted_fraction"]], on="pairid", how="inner")
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


def panel_b(ax):
    inv_hi = pd.read_csv(D / "syn2b_inverted_fraction_high_ani_all.tsv", sep="\t").merge(
        pd.read_csv(D / "dnadiff_inverted_fraction_high_ani_all.tsv", sep="\t")[["pairid", "dnadiff_inverted_fraction"]],
        on="pairid", how="inner")
    hi_truth = pd.read_csv(D / "high_ani_truth.tsv", sep="\t")
    inv_hi = inv_hi.merge(hi_truth[["pairid", "anim_ani"]], on="pairid", how="inner")
    sub = inv_hi[inv_hi["anim_ani"] >= 97]
    x2, y2 = sub["syn2b_raw_inverted_fraction"], sub["dnadiff_inverted_fraction"]
    mm = ~(x2.isna() | y2.isna())
    r2, _ = stats.pearsonr(x2[mm], y2[mm])
    slope, intercept, _, _, _ = stats.linregress(x2[mm], y2[mm])
    ax.scatter(x2[mm], y2[mm], alpha=0.25, s=6, c=COLORS["blue"], edgecolors="none")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.plot([0, 1], [intercept, intercept + slope], "-", lw=1.0, c=COLORS["vermillion"])
    ax.set_xlabel("Syn2bANI raw inverted fraction")
    ax.set_ylabel("dnadiff inverted fraction")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    label_panel(ax, "b")
    ax.set_title(f"ANIm \u2265 97% (r = {r2:.3f}, slope = {slope:.3f}, n = {mm.sum()})", fontsize=8)


def panel_c(ax):
    sv = pd.read_csv(D / "sv_comparison_merged.tsv", sep="\t")
    truth = pd.read_csv(D / "truth_50k.tsv", sep="\t")
    sv = sv.merge(truth[["pairid", "anim_ani", "anim_af_qry"]], on="pairid", how="inner")
    sv["n_ctg"] = (sv["synteny_blocks"] - sv["breakpoint_count"]).clip(lower=1)
    sv = sv.dropna(subset=["anim_ani", "breakpoint_count", "dnadiff_breakpoints", "af_query"])
    n_ctg_safe = sv["n_ctg"].replace([np.inf, -np.inf], np.nan)
    rx, ry, n, rp = robust_partial_residuals(
        sv["breakpoint_count"], sv["dnadiff_breakpoints"],
        [sv["anim_ani"], n_ctg_safe]
    )
    ax.scatter(rx, ry, alpha=0.25, s=5, c=COLORS["vermillion"], edgecolors="none")
    ax.axhline(0, color="k", ls="--", lw=0.8)
    ax.axvline(0, color="k", ls="--", lw=0.8)
    lim_x = np.percentile(np.abs(rx), 99)
    lim_y = np.percentile(np.abs(ry), 99)
    ax.set_xlim(-lim_x, lim_x)
    ax.set_ylim(-lim_y, lim_y)
    ax.set_xlabel("Breakpoint count residual")
    ax.set_ylabel("dnadiff breakpoints residual")
    label_panel(ax, "c")
    ax.set_title(f"Partial correlation (\u03c1 = {rp:.3f}, n = {n})", fontsize=8)


def panel_d(ax):
    sv = pd.read_csv(D / "sv_comparison_merged.tsv", sep="\t")
    truth = pd.read_csv(D / "truth_50k.tsv", sep="\t")
    sv = sv.merge(truth[["pairid", "anim_ani", "anim_af_qry"]], on="pairid", how="inner")
    x4, y4 = sv["af_query"], sv["anim_af_qry"] / 100.0
    m4 = ~(x4.isna() | y4.isna())
    r4, _ = stats.pearsonr(x4[m4], y4[m4])
    ax.scatter(x4[m4], y4[m4], alpha=0.15, s=4, c=COLORS["bluish_green"], edgecolors="none")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("Syn2bANI af_query")
    ax.set_ylabel("dnadiff aligned fraction")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    label_panel(ax, "d")
    ax.set_title(f"Aligned fraction (r = {r4:.3f}, n = {m4.sum()})", fontsize=8)


def panel_e(ax, sv_summary, s2b_sk):
    """Bar chart of wall time at the largest subset (n=22)."""
    n = 22
    dnadiff = sv_summary[(sv_summary["mode"] == "dnadiff") & (sv_summary["n_genomes"] == n)]["dnadiff_wall_s"].values[0]
    skani_dna = sv_summary[(sv_summary["mode"] == "skani_dnadiff") & (sv_summary["n_genomes"] == n)]["total_wall_s"].values[0]
    s2b = s2b_sk[s2b_sk["n_genomes"] == n]["total_wall_s"].values[0]

    labels = ["dnadiff", "skani +\ndnadiff", "syn2bANI\n(ANI+SV)"]
    values = [dnadiff, skani_dna, s2b]
    colors = [COLORS["vermillion"], COLORS["reddish_purple"], COLORS["blue"]]
    bars = ax.bar(labels, values, color=colors, edgecolor="white", linewidth=0.4)
    ax.set_ylabel("Wall time (s)")
    ax.set_yscale("log")
    # Annotate bars.
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, val * 1.2, f"{val:.1f}s",
                ha="center", va="bottom", fontsize=7)
    label_panel(ax, "e")
    ax.set_title(f"SV workflow time (n={n}, {n*n} pairs)", fontsize=8)


def panel_f(ax, sv_summary, s2b_sk):
    """Speedup of syn2bANI (ANI+SV) over dnadiff across subset sizes."""
    dnadiff = sv_summary[sv_summary["mode"] == "dnadiff"][["n_pairs", "dnadiff_wall_s"]].sort_values("n_pairs")
    s2b = s2b_sk[["n_pairs", "total_wall_s"]].sort_values("n_pairs")
    merged = dnadiff.merge(s2b, on="n_pairs", suffixes=("_dna", "_s2b"))
    merged["speedup"] = merged["dnadiff_wall_s"] / merged["total_wall_s"]

    ax.plot(merged["n_pairs"], merged["speedup"], marker="o", color=COLORS["blue"],
            lw=1.2, ms=5, markeredgecolor="white", markeredgewidth=0.4)
    ax.axhline(1, color="k", ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of pairs (all-vs-all, n\u00b2)")
    ax.set_ylabel("Speedup vs. dnadiff")
    ax.set_xlim(3, 700)
    label_panel(ax, "f")
    ax.set_title("syn2bANI (ANI+SV) speedup", fontsize=8)


def main():
    sv_summary, s2b_sk = load_sv_efficiency()

    fig, axes = plt.subplots(3, 2, figsize=figure_size(17.8, aspect=1.28))

    panel_a(axes[0, 0])
    panel_b(axes[0, 1])
    panel_c(axes[1, 0])
    panel_d(axes[1, 1])
    panel_e(axes[2, 0], sv_summary, s2b_sk)
    panel_f(axes[2, 1], sv_summary, s2b_sk)

    plt.tight_layout()
    save_figure(fig, OUT)


if __name__ == "__main__":
    main()
