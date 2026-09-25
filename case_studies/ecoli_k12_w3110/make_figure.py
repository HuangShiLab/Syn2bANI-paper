#!/usr/bin/env python3
"""K-12 MG1655 vs W3110 flagship vignette figure.

Two panels: (a) dotplot of the minimap2 asm20 alignment with the 782-kb
inversion and Syn2b junction positions; (b) what each tool reports for the
same pair (skani / ANI-only vs Syn2b + Syn2bANI struct).
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"

GENOME_LEN = 4_641_652  # MG1655 reference length
MB = 1e6


def parse_paf(path):
    blocks = []
    for line in open(path):
        f = line.split("\t")
        if len(f) < 12:
            continue
        # rstart,rend on query (MG1655 is query here), qstart,qend on ref (W3110)
        rstart, rend = int(f[2]), int(f[3])
        strand = f[4]
        qstart, qend = int(f[7]), int(f[8])
        blocks.append((rstart, rend, qstart, qend, strand))
    return blocks


def parse_junctions(path):
    pos = []
    for line in open(path):
        if line.startswith("#") or line.startswith("genome_A"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and f[2].isdigit():
            pos.append(int(f[2]))
    return pos


def parse_struct(path):
    svs = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        svs.append((f[3], int(f[1]), int(f[2])))  # name,start,end (BED)
    return svs


def main():
    blocks = parse_paf(OUT / "minimap2.paf")
    junctions = parse_junctions(OUT / "syn2b_synteny.junctions.tsv")
    svs = parse_struct(OUT / "struct_bed.tsv")
    inv = [s for s in svs if s[0].startswith("INV")]
    indels = [s for s in svs if not s[0].startswith("INV")]

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10})
    fig = plt.figure(figsize=(11.5, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.18)

    # --- panel (a): dotplot -------------------------------------------------
    ax = fig.add_subplot(gs[0])
    for rstart, rend, qstart, qend, strand in blocks:
        col = "#d62728" if strand == "-" else "#1f77b4"
        lw = 2.5 if strand == "-" else 1.6
        ax.plot([rstart / MB, rend / MB], [qstart / MB, qend / MB],
                color=col, lw=lw, solid_capstyle="round", zorder=3)
    for p in junctions:
        ax.axvline(p / MB, color="darkorange", ls="--", lw=1.1, zorder=2)
    for name, s, e in indels:
        ax.axvline(s / MB, color="gray", ls=":", lw=0.8, alpha=0.7, zorder=1)
    if inv:
        _, s, e = inv[0]
        ax.annotate("782-kb inversion\n(~700 genes, incl. rrnD–rrnE)",
                    xy=(0.5 * (s + e) / MB, 3.9), xytext=(2.6, 0.7),
                    fontsize=8.5, color="#d62728",
                    arrowprops=dict(arrowstyle="->", color="#d62728", lw=0.9))
    ax.set_xlabel("E. coli K-12 MG1655 (Mb)")
    ax.set_ylabel("E. coli K-12 W3110 (Mb)")
    ax.set_title("(a) Whole-genome alignment, MG1655 × W3110", loc="left")
    ax.set_xlim(0, 4.7)
    ax.set_ylim(0, 4.9)
    ax.set_aspect("equal")
    handles = [
        Line2D([], [], color="#1f77b4", lw=1.6, label="collinear block"),
        Line2D([], [], color="#d62728", lw=2.5, label="inverted block"),
        Line2D([], [], color="darkorange", ls="--", lw=1.1,
               label="Syn2b junction"),
        Line2D([], [], color="gray", ls=":", lw=0.9, label="indel (struct)"),
    ]
    ax.legend(handles=handles, fontsize=7.5, frameon=False, loc="upper left")

    # --- panel (b): what each tool reports ----------------------------------
    ax = fig.add_subplot(gs[1])
    ax.axis("off")

    def card(y, title, lines, face):
        h = 0.16 + 0.075 * len(lines)
        box = FancyBboxPatch((0.02, y - h), 0.96, h,
                             boxstyle="round,pad=0.012", fc=face, ec="none",
                             transform=ax.transAxes, zorder=1)
        ax.add_patch(box)
        ax.text(0.05, y - 0.035, title, transform=ax.transAxes,
                fontweight="bold", fontsize=9.5, va="top")
        for i, t in enumerate(lines):
            ax.text(0.07, y - 0.10 - 0.075 * i, t, transform=ax.transAxes,
                    fontsize=9, va="top", family="DejaVu Sans Mono")

    card(1.00, "skani / FastANI (ANI-only)",
         ["ANI 99.99   AF 100% / 100%",
          "no structural readout"],
         "#eef3fb")
    card(0.62, "Syn2b (structural metrics)",
         ["ANI ~99.99        breakpoints = 2",
          "inverted fraction = 0.160",
          "junctions: 3,423,157 / 4,207,508"],
         "#fdf3e7")
    card(0.28, "Syn2bANI struct (SV calls)",
         ["INV  3,428,880–4,211,231  (782 kb)",
          "9 short indels (IS-element losses)"],
         "#e9f6ec")
    ax.set_title("(b) Same pair, three readouts", loc="left", fontsize=10)

    figdir = HERE.parents[1] / "figures"
    figdir.mkdir(exist_ok=True)
    fig.savefig(figdir / "fig_k12_w3110_vignette.png", dpi=300,
                bbox_inches="tight")
    fig.savefig(figdir / "fig_k12_w3110_vignette.pdf", bbox_inches="tight")
    print("figure written")


if __name__ == "__main__":
    main()
