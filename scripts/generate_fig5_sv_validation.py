#!/usr/bin/env python3
"""Generate Figure 5: SV outputs agree with alignment-based truth."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

D = "results/gtdb50k"


def partial_residuals(x, y, ctrls):
    """Return residuals of x and y after OLS on standardized controls."""
    m = ~(x.isna() | y.isna())
    for c in ctrls:
        m &= ~c.isna()
    Z = np.column_stack([np.ones(m.sum())] + [((c[m].values - c[m].mean()) / c[m].std()) for c in ctrls])
    rx = x[m].values - Z @ np.linalg.lstsq(Z, x[m].values, rcond=None)[0]
    ry = y[m].values - Z @ np.linalg.lstsq(Z, y[m].values, rcond=None)[0]
    return rx, ry, int(m.sum())


def main():
    # Load inverted fraction data
    s2b = pd.read_csv(f"{D}/syn2b_inverted_fraction_50k.tsv", sep="\t")
    dna = pd.read_csv(f"{D}/dnadiff_inverted_fraction.tsv", sep="\t")
    inv = s2b.merge(dna[["pairid", "dnadiff_inverted_fraction"]], on="pairid", how="inner")

    # Load SV comparison + truth for breakpoint / AF
    sv = pd.read_csv(f"{D}/sv_comparison_merged.tsv", sep="\t")
    truth = pd.read_csv(f"{D}/truth_50k.tsv", sep="\t")
    sv = sv.merge(truth[["pairid", "anim_ani", "anim_af_qry"]], on="pairid", how="inner")
    sv["n_ctg"] = sv["synteny_blocks"] - sv["breakpoint_count"]
    sv = sv.dropna(subset=["anim_ani", "breakpoint_count", "dnadiff_breakpoints", "af_query"])

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # (a) Inverted fraction all pairs
    ax = axes[0, 0]
    x, y = inv["syn2b_raw_inverted_fraction"], inv["dnadiff_inverted_fraction"]
    m = ~(x.isna() | y.isna())
    r, p = stats.pearsonr(x[m], y[m])
    ax.hexbin(x[m], y[m], gridsize=50, cmap="Blues", mincnt=1)
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("Syn2bANI raw inverted fraction")
    ax.set_ylabel("dnadiff inverted fraction")
    ax.set_title(f"(a) All GTDB-R207 pairs (r = {r:.3f}, n = {m.sum()})")

    # (b) High-ANI zoom
    ax = axes[0, 1]
    inv_hi = pd.read_csv(f"{D}/syn2b_inverted_fraction_high_ani.tsv", sep="\t").merge(
        pd.read_csv(f"{D}/dnadiff_inverted_fraction_high_ani.tsv", sep="\t")[["pairid", "dnadiff_inverted_fraction"]],
        on="pairid", how="inner")
    x2, y2 = inv_hi["syn2b_inverted_fraction"], inv_hi["dnadiff_inverted_fraction"]
    mm = ~(x2.isna() | y2.isna())
    r2, _ = stats.pearsonr(x2[mm], y2[mm])
    ax.scatter(x2[mm], y2[mm], alpha=0.3, s=10)
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("Syn2bANI inverted fraction")
    ax.set_ylabel("dnadiff inverted fraction")
    ax.set_title(f"(b) ANIm ≥ 95% high-ANI set (r = {r2:.3f}, n = {mm.sum()})")

    # (c) Breakpoint count partial correlation
    ax = axes[1, 0]
    rx, ry, n = partial_residuals(sv["breakpoint_count"], sv["dnadiff_breakpoints"],
                                  [sv["anim_ani"], sv["n_ctg"]])
    rp, _ = stats.pearsonr(rx, ry)
    # Clip extreme residuals for visualization only
    rx_clip = np.clip(rx, -500, 2000)
    ry_clip = np.clip(ry, -2000, 5000)
    ax.scatter(rx_clip, ry_clip, alpha=0.3, s=10)
    ax.axhline(0, color="k", ls="--", lw=1)
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.set_xlim(-600, 2200)
    ax.set_ylim(-2500, 5500)
    ax.set_xlabel("Breakpoint count residual (after ANIm + contigs)")
    ax.set_ylabel("dnadiff breakpoints residual")
    ax.set_title(f"(c) Partial correlation (r = {rp:.3f}, n = {n})")

    # (d) af_query vs dnadiff aligned fraction
    ax = axes[1, 1]
    x4, y4 = sv["af_query"], sv["anim_af_qry"] / 100.0
    m4 = ~(x4.isna() | y4.isna())
    r4, _ = stats.pearsonr(x4[m4], y4[m4])
    ax.scatter(x4[m4], y4[m4], alpha=0.3, s=10)
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("Syn2bANI af_query")
    ax.set_ylabel("dnadiff aligned fraction")
    ax.set_title(f"(d) Aligned fraction (r = {r4:.3f}, n = {m4.sum()})")

    plt.tight_layout()
    out = "paper/figures/main/fig5_sv_validation.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
