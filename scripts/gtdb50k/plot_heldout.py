#!/usr/bin/env python3
"""Publication-grade figure for the GTDB 50k held-out benchmark.

Panels:
  a. syn2bani calibrated v5 vs ANIm (hexbin, common subset)
  b. skani vs ANIm (hexbin, common subset)
  c. per-band MAE bars (calibrated vs skani, common subset)
  d. error-vs-truth distributions (calibrated vs skani)

Output: figures/report/fig_gtdb50k_heldout.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
RES = os.path.join(ROOT, "results", "gtdb50k")
OUT = os.path.join(ROOT, "figures", "report", "fig_gtdb50k_heldout.png")

BAND_ORDER = ["80-85", "85-90", "90-95", "95-100"]

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
})

def main():
    truth = pd.read_csv(os.path.join(RES, "truth_50k.tsv"), sep="\t")
    s2b = pd.read_csv(os.path.join(RES, "s2b_50k.tsv"), sep="\t")
    pairs = pd.read_csv(os.path.join(RES, "pairs_50k.tsv"), sep="\t")
    pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]
    df = truth.merge(pairs[["pairid", "skani_ani", "band"]], on="pairid")
    df = df.merge(s2b[["pairid", "ani_cal"]], on="pairid")
    df = df[np.isfinite(df["ani_cal"])]

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4))
    lo, hi = 79.5, 100.5

    for ax, col, title, mae in (
        (axes[0, 0], "ani_cal", "a  Syn2bANI calibrated (v5)", 0.619),
        (axes[0, 1], "skani_ani", "b  skani", 0.957),
    ):
        hb = ax.hexbin(df["anim_ani"], df[col], gridsize=120, extent=(lo, hi, lo, hi),
                       bins="log", mincnt=1, cmap="viridis")
        ax.plot([lo, hi], [lo, hi], color="red", lw=0.8, ls="--")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel("ANIm (dnadiff) [%]")
        ax.set_ylabel("estimated ANI [%]")
        ax.set_title(f"{title}   MAE {mae}", loc="left", fontweight="bold")
        ax.set_aspect("equal")
        fig.colorbar(hb, ax=ax, label="pairs (log)", shrink=0.85)

    # c: per-band MAE
    ax = axes[1, 0]
    mae_cal = [np.mean(np.abs(df.loc[df.band == b, "ani_cal"] - df.loc[df.band == b, "anim_ani"])) for b in BAND_ORDER]
    mae_sk = [np.mean(np.abs(df.loc[df.band == b, "skani_ani"] - df.loc[df.band == b, "anim_ani"])) for b in BAND_ORDER]
    x = np.arange(len(BAND_ORDER)); w = 0.36
    ax.bar(x - w/2, mae_cal, w, label="Syn2bANI cal", color="#219ebc")
    ax.bar(x + w/2, mae_sk, w, label="skani", color="#fb8500")
    for xi, v in zip(x - w/2, mae_cal):
        ax.text(xi, v + 0.03, f"{v:.2f}", ha="center", fontsize=7.5)
    for xi, v in zip(x + w/2, mae_sk):
        ax.text(xi, v + 0.03, f"{v:.2f}", ha="center", fontsize=7.5)
    ax.set_xticks(x, BAND_ORDER)
    ax.set_xlabel("ANI band [%]")
    ax.set_ylabel("MAE vs ANIm [ANI points]")
    ax.set_title("c  Per-band MAE (held-out)", loc="left", fontweight="bold")
    ax.legend(frameon=False)
    ax.set_ylim(0, max(mae_sk) * 1.18)

    # d: error distributions
    ax = axes[1, 1]
    e_cal = (df["ani_cal"] - df["anim_ani"]).to_numpy()
    e_sk = (df["skani_ani"] - df["anim_ani"]).to_numpy()
    bins = np.linspace(-5, 5, 161)
    ax.hist(e_sk, bins=bins, density=True, histtype="step", lw=1.2, color="#fb8500",
            label=f"skani (bias {np.mean(e_sk):+.2f})")
    ax.hist(e_cal, bins=bins, density=True, histtype="step", lw=1.2, color="#219ebc",
            label=f"Syn2bANI cal (bias {np.mean(e_cal):+.2f})")
    ax.axvline(0, color="grey", lw=0.7, ls=":")
    ax.set_xlabel("estimate − ANIm [ANI points]")
    ax.set_ylabel("density")
    ax.set_title("d  Signed-error distribution", loc="left", fontweight="bold")
    ax.legend(frameon=False)

    fig.suptitle("Held-out GTDB-R207 benchmark: 43,334 pairs, training genomes excluded "
                 f"(scored subset n = {len(df):,})", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=300)
    print("wrote", OUT)

if __name__ == "__main__":
    main()
