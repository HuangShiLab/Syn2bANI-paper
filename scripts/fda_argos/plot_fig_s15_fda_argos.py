#!/usr/bin/env python3
"""Generate Supplementary Figure S15 for FDA-ARGOS E. coli comparisons."""

from __future__ import annotations

from pathlib import Path

import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / "case_studies/fda_argos"
OUT = ROOT / "paper/figures/supplementary/fig_s15_fda_argos_ecoli"


def read_paf_chromosome(path: Path, min_length: int = 100_000):
    rows = []
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            block_len = int(f[10])
            if block_len < min_length:
                continue
            rows.append({
                "q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]),
                "strand": f[4], "r": f[5], "rlen": int(f[6]), "rs": int(f[7]),
                "re": int(f[8]), "nmatch": block_len,
            })
    return pd.DataFrame(rows).sort_values("qs")


def main() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
        "legend.fontsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "pdf.fonttype": 42, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False,
    })

    pairs = pd.read_csv(CASE / "pairs99.9_ecoli.tsv", sep="\t")
    summary = pd.read_csv(CASE / "validate/validation_summary.tsv", sep="\t")
    discordant = pd.read_csv(CASE / "amr_discordance.tsv", sep="\t")
    dkeys = {
        tuple(sorted((row.acc1, row.acc2)))
        for row in discordant.itertuples(index=False)
    }
    pairs["amr_diff"] = [
        tuple(sorted((row.acc1, row.acc2))) in dkeys
        for row in pairs.itertuples(index=False)
    ]

    fig = plt.figure(figsize=(180 / 25.4, 135 / 25.4))
    gs = fig.add_gridspec(
        2, 2, left=0.075, right=0.975, top=0.93, bottom=0.09,
        wspace=0.28, hspace=0.46,
    )

    ax = fig.add_subplot(gs[0, 0])
    ax.hist(pairs["ani"], bins=np.arange(99.895, 100.005, 0.005),
            color="#4C78A8", edgecolor="white", linewidth=0.3)
    ax.set_xlabel("skani ANI (%)")
    ax.set_ylabel("Reference-genome pairs")
    ax.set_title("a  Very-high-ANI pairs", loc="left", weight="bold")
    ax.text(0.03, 0.90, f"n = {len(pairs)}", transform=ax.transAxes)

    ax = fig.add_subplot(gs[0, 1])
    order = ["global_flip_or_rotation", "partial_inversion", "collinear_local"]
    labels = ["Global flip /\nrotation", "Partial\ninversion", "Collinear /\nlocal"]
    counts = summary["class"].value_counts().reindex(order).fillna(0)
    bars = ax.bar(labels, counts, color=["#BBBBBB", "#B2182B", "#4C78A8"],
                  edgecolor="white", linewidth=0.3)
    for bar, n in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.7,
                f"{int(n)}", ha="center", va="bottom")
    ax.set_ylabel("Pairs")
    ax.set_ylim(0, 54)
    ax.set_title("b  Minimap2 structural validation", loc="left", weight="bold")

    ax = fig.add_subplot(gs[1, 0])
    thresholds = [99.9, 99.99, 99.999]
    totals = [int((pairs["ani"] >= x).sum()) for x in thresholds]
    diffs = [int(((pairs["ani"] >= x) & pairs["amr_diff"]).sum())
             for x in thresholds]
    xpos = np.arange(len(thresholds))
    ax.bar(xpos - 0.18, totals, width=0.34, color="#D9D9D9",
           edgecolor="white", linewidth=0.3, label="Pairs with AMR data")
    ax.bar(xpos + 0.18, diffs, width=0.34, color="#B2182B",
           edgecolor="white", linewidth=0.3, label="Core AMR-genotype discordant")
    for x, total, diff in zip(xpos, totals, diffs):
        ax.text(x - 0.18, total + 1.2, str(total), ha="center", va="bottom")
        ax.text(x + 0.18, diff + 1.2, str(diff), ha="center", va="bottom")
    ax.set_xticks(xpos, ["ANI ≥ 99.9", "ANI ≥ 99.99", "ANI ≥ 99.999"])
    ax.set_ylabel("Pairs")
    ax.set_ylim(0, 118)
    ax.legend(loc="upper right")
    ax.set_title("c  AMRFinderPlus core-genotype discordance", loc="left",
                 weight="bold")

    ax = fig.add_subplot(gs[1, 1])
    paf = CASE / "validate/pafs/GCF_002209105.2__GCF_006364695.1.paf"
    aln = read_paf_chromosome(paf)
    qmax = int(aln["qlen"].iloc[0])
    rmax = int(aln["rlen"].iloc[0])
    for row in aln.itertuples(index=False):
        color = "#4C78A8" if row.strand == "+" else "#B2182B"
        ax.plot([row.qs / 1e6, row.qe / 1e6],
                [row.rs / 1e6, row.re / 1e6], lw=3.2, solid_capstyle="butt",
                color=color, alpha=0.88)
    ax.set_xlim(0, qmax / 1e6)
    ax.set_ylim(0, rmax / 1e6)
    ax.set_xlabel("FDAARGOS_348 coordinate (Mb)")
    ax.set_ylabel("FDAARGOS_772 coordinate (Mb)")
    ax.set_title("d  ANI = 100.0 pair: alignment blocks", loc="left",
                 weight="bold")
    ax.text(0.97, 0.05, "plus-strand blocks (blue); ≥100 kb",
            transform=ax.transAxes, color="0.25", ha="right", va="bottom")

    fig.suptitle(
        "FDA-ARGOS E. coli: structure and resistance genotypes at very high ANI",
        x=0.075, ha="left", y=0.985, fontsize=8, weight="bold",
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".png"), dpi=600, bbox_inches="tight")

    qa = {
        "pairs": int(len(pairs)),
        "class_counts": {key: int(counts[key]) for key in order},
        "amr_bands": {
            str(threshold): {"total": total, "discordant": diff}
            for threshold, total, diff in zip(thresholds, totals, diffs)
        },
        "panel_d_pair": "GCF_002209105.2__GCF_006364695.1",
        "panel_d_chromosome_blocks_ge100kb": int(len(aln)),
    }
    Path(str(OUT) + ".qa.json").write_text(json.dumps(qa, indent=2) + "\n")
    print(f"wrote {OUT}.pdf/png; QA counts: {qa}")


if __name__ == "__main__":
    main()
