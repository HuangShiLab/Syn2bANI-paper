#!/usr/bin/env python3
"""Generate Figure 6: ANI and SV are decoupled in GTDB-R207."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

D = "results/gtdb50k"


def main():
    inv = pd.read_csv(f"{D}/syn2b_inverted_fraction_50k.tsv", sep="\t")
    dna = pd.read_csv(f"{D}/dnadiff_inverted_fraction.tsv", sep="\t")
    inv = inv.merge(dna[["pairid", "dnadiff_inverted_fraction"]], on="pairid", how="inner")
    sv = pd.read_csv(f"{D}/sv_comparison_merged.tsv", sep="\t")
    truth = pd.read_csv(f"{D}/truth_50k.tsv", sep="\t")
    df = inv.merge(truth[["pairid", "anim_ani"]], on="pairid", how="inner").merge(
        sv[["pairid", "breakpoint_count"]], on="pairid", how="inner")
    df = df.dropna()

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # (a) ANIm vs inverted fraction
    ax = axes[0, 0]
    ax.scatter(df["anim_ani"], df["dnadiff_inverted_fraction"], alpha=0.2, s=5)
    ax.set_xlabel("ANIm (%)")
    ax.set_ylabel("dnadiff inverted fraction")
    ax.set_title("(a) ANIm vs inverted fraction (GTDB-R207)")

    # (b) ANIm vs breakpoint count
    ax = axes[0, 1]
    ax.scatter(df["anim_ani"], df["breakpoint_count"], alpha=0.2, s=5)
    ax.set_xlabel("ANIm (%)")
    ax.set_ylabel("Syn2bANI breakpoint count")
    ax.set_title("(b) ANIm vs breakpoint count")

    # (c) inverted fraction by ANI band
    ax = axes[1, 0]
    bands = [(80, 85), (85, 90), (90, 95), (95, 97), (97, 99), (99, 101)]
    vals = [df[(df["anim_ani"] >= lo) & (df["anim_ani"] < hi)]["dnadiff_inverted_fraction"].values
            for lo, hi in bands]
    labels = [f"{lo}-{hi}" for lo, hi in bands]
    bp = ax.boxplot(vals, labels=labels, showfliers=False, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("lightblue")
    ax.set_xlabel("ANIm band (%)")
    ax.set_ylabel("dnadiff inverted fraction")
    ax.set_title("(c) Inverted fraction distribution by ANI band")

    # (d) high-ANI high-SV examples
    ax = axes[1, 1]
    disc = df[(df["anim_ani"] >= 99) & (df["dnadiff_inverted_fraction"] >= 0.3)]
    ax.scatter(disc["anim_ani"], disc["dnadiff_inverted_fraction"], color="red", s=20)
    ax.set_xlabel("ANIm (%)")
    ax.set_ylabel("dnadiff inverted fraction")
    ax.set_title(f"(d) Discordant pairs: ANI ≥ 99% & inv-frac ≥ 0.3 (n = {len(disc)})")

    plt.tight_layout()
    out = "paper/figures/main/fig6_ani_sv_discordance.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
