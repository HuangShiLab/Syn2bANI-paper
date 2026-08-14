#!/usr/bin/env python3
"""Generate publication-grade report figures for the syn2bani methods/validation report.

Figures
-------
F1  Simulation-ladder accuracy (syn2bani only; cross-tool overlay gated, see
    ADD_CROSS_TOOL_LADDER / add_cross_tool_ladder).
F2  Robustness: indel sweep, fragmentation, GC sweep, accessory fraction.
F3  Enzyme-panel optimization: per-enzyme bias, panel-size trade-off, tag
    composition vs GC.
F4  Mid-ANI validation against ANIm truth (15 pairs, all INCONSISTENT).
F5  GTDB R207 large-scale benchmark vs FastANI (672 pairs where FastANI reports).
F6  Computational efficiency.
F7  ANIm-truth benchmark by ANI band (2,074 GTDB R207 pairs; ridge-band-holdout
    calibration of the current 4-enzyme panel).

Pure pandas/matplotlib plotting from existing files; no heavy compute.
Run:  python3 analysis/plot_report_figures.py
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Paths (the only two roots used)
# --------------------------------------------------------------------------
PAPER = Path("/Users/macstudio/Downloads/Syn2bANI-paper")
S2B = Path("/Users/macstudio/Downloads/Syn2bANI")
FIGDIR = PAPER / "figures" / "report"

F1_RESULTS = S2B / "prototype" / "simindel_results_4e.tsv"
F1_MANIFEST = S2B / "prototype" / "simindel" / "manifest.tsv"
F2A_RESULTS = S2B / "prototype" / "simindel_sweep_results_4e.tsv"
F2A_MANIFEST = S2B / "prototype" / "simindel_sweep" / "manifest.tsv"
F2B_RESULTS = S2B / "prototype" / "simfrag_results_4e.tsv"
F2B_MANIFEST = S2B / "prototype" / "simfrag" / "manifest.tsv"
F4_EVAL = PAPER / "results" / "validation_mid_ani_anim" / "anim_4e" / "anim_midani_evaluation.tsv"
F4_METRICS = PAPER / "results" / "validation_mid_ani_anim" / "anim_4e" / "anim_midani_metrics.tsv"
F5_MATRIX = PAPER / "results" / "matrix_gtdb_r207_100k_v8_final.tsv"
F5_SUMMARY_OUT = FIGDIR / "gtdb_metrics_summary.tsv"
F6_RUNTIME = PAPER / "results" / "efficiency_v8" / "runtime_scaling.tsv"
F6_SKETCH = PAPER / "results" / "efficiency_v8" / "sketch_benchmark.tsv"
F7_TABLE = PAPER / "results" / "panel_by_band" / "anim_main_table.tsv"
F7_PREDS = PAPER / "results" / "panel_by_band" / "ridge_cv_preds_4e.tsv"

# Future cross-tool ladder TSV (columns: name, skani_ani, fastani_ani) on the same
# simulated genomes as F1. Does not exist yet; the code path is written and gated.
CROSS_TOOL_LADDER_TSV = S2B / "prototype" / "simindel_cross_tool_4e.tsv"
ADD_CROSS_TOOL_LADDER = True  # safe: add_cross_tool_ladder() no-ops if the file is absent

# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.labelsize": 8,
    "legend.fontsize": 6.5,
    "legend.frameon": False,
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
})

# Okabe-Ito tool colors
C_UNIFORM = "#0072B2"   # syn2bani uniform (blue)
C_GAMMA = "#D55E00"     # syn2bani gamma (vermillion)
C_SKANI = "#009E73"     # skani (green)
C_FASTANI = "#CC79A7"   # FastANI (reddish purple)
C_TRUTH = "#000000"     # ANIm / identity references
C_GREY = "#7F7F7F"

L_UNIFORM = "syn2bani (uniform)"
L_GAMMA = "syn2bani (gamma)"
L_SKANI = "skani"
L_FASTANI = "FastANI"

# --------------------------------------------------------------------------
# Hardcoded tables transcribed from the validation docs (sources in comments)
# --------------------------------------------------------------------------

# F2(c): ALGORITHM_MLE.md section 4.8, table "MAE (current 4) / MAE (balanced 5)"
GC_SWEEP = {
    # genome: (GC%, MAE current-4 panel, MAE balanced-5 panel)
    "F. nucleatum": (27.2, 0.162, 0.125),
    "S. mutans": (36.8, 0.135, 0.063),
    "E. coli K-12": (50.8, 0.074, 0.066),
    "B. longum": (60.1, 0.356, 0.312),
    "S. coelicolor": (72.1, 0.166, 0.200),
}

# F2(d): V8_MLE_VALIDATION.md section 3.2 — accessory-fraction混淆 experiment,
# true ANI fixed at 95.000. AF tracks 1 - F (doc: error <= 0.004).
ACCESSORY = {
    # accessory fraction (%): (estimate, error, AF)
    0: (95.044, 0.044, 1.000),
    10: (95.106, 0.106, 0.896),
    20: (95.074, 0.074, 0.796),
    30: (95.071, 0.071, 0.698),
    40: (95.250, 0.250, 0.597),
    50: (95.140, 0.140, 0.496),
}

# F3(a): V8_MLE_VALIDATION.md section 3.12 — per-enzyme divergence bias,
# truth 95.000, 10 replicates, uniform per-site divergence.
ENZYME_BIAS = {
    # enzyme: (bias, SD, significance in sigma)
    "AlfI": (0.273, 0.096, 9.0),
    "FalI": (-0.586, 0.211, 8.8),
    "BcgI": (0.095, 0.104, 2.9),
    "AloI": (0.011, 0.142, 0.2),
}

# F3(b): ALGORITHM_MLE.md section 4.7.1 — enzyme-panel size trade-off.
PANEL_SIZE = {
    # n_enzymes label: (tags in chains, MAE vs skani on real genomes, MAE simulated)
    "4": (4874, 0.094, 0.074),
    "6": (6655, 0.167, 0.062),
    "8": (7314, 0.141, 0.071),
    "8 (CjeI)": (22141, 0.290, 0.038),
    "11": (120717, 0.670, 0.133),
}

# F3(c): ALGORITHM_MLE.md section 4.8 — tag density / per-enzyme counts per Mb.
TAG_COMPOSITION = {
    # genome: (GC%, tags/Mb, BcgI, AlfI, AloI, FalI)
    "F. nucleatum": (27.2, 610, 22, 31, 79, 478),
    "S. mutans": (36.8, 1046, 174, 173, 75, 624),
    "E. coli K-12": (50.8, 1339, 632, 436, 113, 158),
    "B. longum": (60.1, 2101, 1173, 511, 222, 195),
    "S. coelicolor": (72.1, 1708, 979, 388, 223, 118),
}

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def panel_label(ax, letter):
    ax.text(-0.18, 1.10, f"({letter})", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="left")


def save(fig, name):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    png = FIGDIR / f"{name}.png"
    pdf = FIGDIR / f"{name}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png.name} / {pdf.name}")


def identity_line(ax, lo, hi):
    ax.plot([lo, hi], [lo, hi], color=C_TRUTH, lw=0.8, ls="--", zorder=1)


# --------------------------------------------------------------------------
# F1 — simulation ladder accuracy
# --------------------------------------------------------------------------

def load_ladder():
    res = pd.read_csv(F1_RESULTS, sep="\t")
    man = pd.read_csv(F1_MANIFEST, sep="\t")
    df = res.merge(man[["name", "true_ani"]], left_on="query", right_on="name", how="left")
    df["truth"] = df["true_ani"] * 100.0
    return df.sort_values("truth")


def add_cross_tool_ladder(ax_scatter, ax_error, tsv_path=CROSS_TOOL_LADDER_TSV):
    """Overlay skani/FastANI points on the F1 ladder panels.

    Expects a TSV with columns: name, skani_ani, fastani_ani, measured on the
    same simulated genomes as F1_MANIFEST (true ANI taken from the manifest).
    Returns the merged DataFrame (with a `truth` column in percent), or None if
    the file does not exist yet.
    """
    tsv_path = Path(tsv_path)
    if not tsv_path.exists():
        print(f"  [F1] cross-tool ladder TSV not present, skipping overlay: {tsv_path}")
        return None
    cross = pd.read_csv(tsv_path, sep="\t")
    man = pd.read_csv(F1_MANIFEST, sep="\t")
    df = cross.merge(man[["name", "true_ani"]], on="name", how="left")
    df["truth"] = df["true_ani"] * 100.0
    ax_scatter.scatter(df["truth"], df["skani_ani"], color=C_SKANI, marker="s",
                       label=L_SKANI, zorder=3)
    ax_scatter.scatter(df["truth"], df["fastani_ani"], color=C_FASTANI, marker="D",
                       label=L_FASTANI, zorder=3)
    ax_error.scatter(df["truth"], df["skani_ani"] - df["truth"], color=C_SKANI,
                     marker="s", label=L_SKANI, zorder=3)
    ax_error.scatter(df["truth"], df["fastani_ani"] - df["truth"], color=C_FASTANI,
                     marker="D", label=L_FASTANI, zorder=3)
    return df


def fig1():
    print("F1: simulation ladder")
    df = load_ladder()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3.1))

    lo, hi = 84, 101
    identity_line(ax1, lo, hi)
    ax1.scatter(df["truth"], df["ani"], color=C_GAMMA, marker="o", label=L_GAMMA, zorder=3)
    ax1.scatter(df["truth"], df["ani_uniform"], color=C_UNIFORM, marker="^",
                label=L_UNIFORM, zorder=3, facecolors="none", edgecolors=C_UNIFORM)
    ax1.set_xlim(lo, hi)
    ax1.set_ylim(lo, hi)
    ax1.set_xlabel("True ANI (%)")
    ax1.set_ylabel("Estimated ANI (%)")
    ax1.grid(True)
    panel_label(ax1, "a")

    err_g = df["ani"] - df["truth"]
    err_u = df["ani_uniform"] - df["truth"]
    ax2.axhline(0, color=C_TRUTH, lw=0.8, ls="--", zorder=1)
    ax2.scatter(df["truth"], err_g, color=C_GAMMA, marker="o", label=L_GAMMA, zorder=3)
    ax2.scatter(df["truth"], err_u, color=C_UNIFORM, marker="^", zorder=3,
                facecolors="none", edgecolors=C_UNIFORM, label=L_UNIFORM)
    ax2.set_xlabel("True ANI (%)")
    ax2.set_ylabel("Error (est − truth, ANI points)")
    ax2.grid(True)

    # Cross-tool overlay (skani/FastANI on the same simulated genomes). The
    # benchmark writes CROSS_TOOL_LADDER_TSV; until it exists this is a no-op.
    mae_lines = [f"MAE gamma = {err_g.abs().mean():.3f}",
                 f"MAE uniform = {err_u.abs().mean():.3f}"]
    if ADD_CROSS_TOOL_LADDER:
        cross = add_cross_tool_ladder(ax1, ax2)
        if cross is not None:
            for col, lab in (("skani_ani", "skani"), ("fastani_ani", "FastANI")):
                err = cross[col] - cross["truth"]
                mae_lines.append(f"MAE {lab} = {err.abs().mean():.3f}")
    # Panel (b) has no free corner once all four tools are plotted; the MAE
    # block goes in panel (a)'s empty upper-left triangle (above identity).
    ax1.text(0.03, 0.97, "\n".join(mae_lines),
             transform=ax1.transAxes, va="top", ha="left", fontsize=6.5)
    panel_label(ax2, "b")

    ax1.legend(loc="lower right")
    ax2.legend(loc="lower right")
    fig.tight_layout()
    save(fig, "fig1_simulation_ladder")


# --------------------------------------------------------------------------
# F2 — robustness (4 panels)
# --------------------------------------------------------------------------

def fig2():
    print("F2: robustness")
    fig, axes = plt.subplots(2, 2, figsize=(7, 5.6))

    # (a) indel sweep -------------------------------------------------------
    ax = axes[0, 0]
    res = pd.read_csv(F2A_RESULTS, sep="\t")
    man = pd.read_csv(F2A_MANIFEST, sep="\t")
    df = res.merge(man[["name", "indel_rate", "true_ani"]],
                   left_on="query", right_on="name").sort_values("indel_rate")
    truth = df["true_ani"] * 100.0
    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--")
    ax.plot(df["indel_rate"], df["ani"] - truth, color=C_GAMMA, marker="o", label=L_GAMMA)
    ax.plot(df["indel_rate"], df["ani_uniform"] - truth, color=C_UNIFORM, marker="^",
            label=L_UNIFORM)
    ax.set_xlabel("Indel rate (per 100 kb)")
    ax.set_ylabel("Error (est − 95.000, ANI points)")
    ax.grid(True)
    ax.legend(loc="best")
    panel_label(ax, "a")

    # (b) fragmentation -----------------------------------------------------
    # NOTE: ~half of the contigs in each draft are reverse-complemented and the
    # contig order is shuffled (manifest `flipped=True`); true ANI 95.000.
    ax = axes[0, 1]
    frag = pd.read_csv(F2B_RESULTS, sep="\t")
    frag["n_contigs"] = frag["query"].str.extract(r"q95_c(\d+)_")[0].astype(int)
    frag = frag.sort_values("n_contigs")
    truth = 95.000
    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--")
    ax.plot(frag["n_contigs"], frag["ani"] - truth, color=C_GAMMA, marker="o", label=L_GAMMA)
    ax.plot(frag["n_contigs"], frag["ani_uniform"] - truth, color=C_UNIFORM, marker="^",
            label=L_UNIFORM)
    ax.set_xscale("log")
    ax.set_xticks([20, 50, 100, 200])
    ax.set_xticklabels([20, 50, 100, 200])
    ax.minorticks_off()
    ax.set_xlabel("Number of contigs")
    ax.set_ylabel("Error (est − 95.000, ANI points)")
    ax.grid(True)
    ax.legend(loc="best")
    panel_label(ax, "b")

    # (c) GC sweep (hardcoded from ALGORITHM_MLE.md section 4.8) ------------
    ax = axes[1, 0]
    gc = np.array([v[0] for v in GC_SWEEP.values()])
    mae4 = np.array([v[1] for v in GC_SWEEP.values()])
    mae5 = np.array([v[2] for v in GC_SWEEP.values()])
    ax.plot(gc, mae4, color=C_GAMMA, marker="o", label="current 4-enzyme panel")
    ax.plot(gc, mae5, color=C_UNIFORM, marker="^", label="balanced 5-enzyme panel")
    for x, y, name in zip(gc, mae4, GC_SWEEP):
        ax.annotate(name, (x, y), textcoords="offset points", xytext=(0, 6),
                    ha="center", fontsize=5.5, color=C_GREY)
    ax.set_xlabel("Genome GC content (%)")
    ax.set_ylabel("MAE vs known truth (ANI points)")
    ax.set_ylim(0, None)
    ax.grid(True)
    ax.legend(loc="upper left")
    panel_label(ax, "c")

    # (d) accessory fraction (hardcoded from V8_MLE_VALIDATION.md §3.2) -----
    ax = axes[1, 1]
    acc = np.array(sorted(ACCESSORY))
    err = np.array([ACCESSORY[a][1] for a in acc])
    af = np.array([ACCESSORY[a][2] for a in acc])
    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--")
    ax.plot(acc, err, color=C_UNIFORM, marker="o", label="ANI error (est − 95.000)")
    ax.set_xlabel("Accessory fraction (%)")
    ax.set_ylabel("Error (ANI points)", color=C_UNIFORM)
    ax.tick_params(axis="y", labelcolor=C_UNIFORM)
    ax.set_ylim(-0.3, 0.35)
    ax.grid(True)
    ax2 = ax.twinx()
    ax2.spines["right"].set_visible(True)
    ax2.plot(acc, af, color=C_GREY, marker="s", ls=":", label="af_query")
    ax2.plot(acc, 1 - acc / 100.0, color=C_TRUTH, lw=0.8, ls="--", label="1 − F (truth)")
    ax2.set_ylabel("Aligned fraction (af_query)", color=C_GREY)
    ax2.tick_params(axis="y", labelcolor=C_GREY)
    ax2.set_ylim(0.4, 1.05)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="lower left", fontsize=6)
    panel_label(ax, "d")

    fig.tight_layout()
    save(fig, "fig2_robustness")


# --------------------------------------------------------------------------
# F3 — enzyme panel optimization (3 panels)
# --------------------------------------------------------------------------

def fig3():
    print("F3: enzyme panel")
    fig, axes = plt.subplots(1, 3, figsize=(7, 2.9))

    # (a) per-enzyme bias (V8_MLE_VALIDATION.md §3.12) ----------------------
    ax = axes[0]
    names = list(ENZYME_BIAS)
    bias = [ENZYME_BIAS[e][0] for e in names]
    sd = [ENZYME_BIAS[e][1] for e in names]
    sig = [ENZYME_BIAS[e][2] for e in names]
    colors = [C_GAMMA if abs(s) >= 3 else C_UNIFORM for s in sig]
    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--")
    ax.bar(names, bias, yerr=sd, capsize=3, color=colors, edgecolor="black",
           linewidth=0.5, error_kw=dict(lw=0.8))
    for i, (b, s) in enumerate(zip(bias, sig)):
        ax.annotate(f"{s:.1f}σ", (i, b), textcoords="offset points",
                    xytext=(0, 20 if b >= 0 else -30), ha="center", fontsize=6)
    ax.set_ylabel("Per-enzyme bias (ANI points)")
    ax.set_xlabel("Enzyme (truth 95.000, n = 10 replicates)")
    ax.set_ylim(-1.1, 0.65)
    ax.grid(True, axis="y")
    panel_label(ax, "a")

    # (b) panel-size trade-off (ALGORITHM_MLE.md §4.7.1) --------------------
    ax = axes[1]
    labels = list(PANEL_SIZE)
    tags = np.array([PANEL_SIZE[k][0] for k in labels], dtype=float)
    mae_real = np.array([PANEL_SIZE[k][1] for k in labels])
    mae_sim = np.array([PANEL_SIZE[k][2] for k in labels])
    ax.plot(tags, mae_real, color=C_GAMMA, marker="o", label="real genomes (vs skani)")
    ax.plot(tags, mae_sim, color=C_UNIFORM, marker="^", label="simulated (known truth)")
    offsets = {"4": ((-6, -12), "right"), "6": ((-4, 9), "center"),
               "8": ((8, -12), "left"), "8 (CjeI)": ((0, 9), "center"),
               "11": ((0, 9), "center")}
    for x, y, lab in zip(tags, mae_real, labels):
        off, ha = offsets[lab]
        ax.annotate(f"{lab} enz.", (x, y), textcoords="offset points",
                    xytext=off, ha=ha, fontsize=6, color=C_GREY)
    ax.set_xscale("log")
    ax.set_xlabel("Tags in chains (log scale)")
    ax.set_ylabel("MAE (ANI points)")
    ax.set_ylim(0, 0.8)
    ax.grid(True)
    ax.legend(loc="upper left", fontsize=6)
    panel_label(ax, "b")

    # (c) tag composition vs GC (ALGORITHM_MLE.md §4.8) ---------------------
    ax = axes[2]
    genomes = list(TAG_COMPOSITION)
    gc = [TAG_COMPOSITION[g][0] for g in genomes]
    enzymes = ["BcgI", "AlfI", "AloI", "FalI"]
    ecolors = {"BcgI": "#0072B2", "AlfI": "#E69F00", "AloI": "#009E73", "FalI": "#CC79A7"}
    counts = {e: np.array([TAG_COMPOSITION[g][2 + i] for g in genomes], dtype=float)
              for i, e in enumerate(enzymes)}
    totals = sum(counts.values())
    x = np.arange(len(genomes))
    bottom = np.zeros(len(genomes))
    for e in enzymes:
        share = 100.0 * counts[e] / totals
        ax.bar(x, share, bottom=bottom, color=ecolors[e], edgecolor="white",
               linewidth=0.4, label=e)
        bottom += share
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g}\n{c}%" for g, c in zip(genomes, gc)], fontsize=5.5,
                       rotation=30, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Share of tags (%)")
    ax.set_ylim(0, 100)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), fontsize=6, ncol=4,
              columnspacing=0.8, handlelength=1.2)
    panel_label(ax, "c")

    fig.tight_layout()
    save(fig, "fig3_enzyme_panel")


# --------------------------------------------------------------------------
# F4 — mid-ANI validation against ANIm truth
# --------------------------------------------------------------------------

F4_METHODS = [  # (column, color, marker, label)
    ("s2b_ani", C_GAMMA, "o", L_GAMMA),
    ("ani_uniform", C_UNIFORM, "^", L_UNIFORM),
    ("skani_ani", C_SKANI, "s", L_SKANI),
    ("fastani_ani", C_FASTANI, "D", L_FASTANI),
]


def fig4():
    print("F4: mid-ANI validation vs ANIm")
    ev = pd.read_csv(F4_EVAL, sep="\t")
    met = pd.read_csv(F4_METRICS, sep="\t").set_index("method")
    n_flagged = (ev["flag"] == "INCONSISTENT").sum()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3.2))

    # (a) scatter vs ANIm truth --------------------------------------------
    lo, hi = 80.5, 92
    identity_line(ax1, lo, hi)
    for col, color, marker, label in F4_METHODS:
        # All 15 pairs are flagged INCONSISTENT: syn2bani points get a black
        # edge to mark that downstream users would not see the raw numbers.
        edge = "black" if col in ("s2b_ani", "ani_uniform") else "none"
        ax1.scatter(ev["anim_ani"], ev[col], color=color, marker=marker,
                    edgecolors=edge, linewidths=0.6, label=label, zorder=3, s=16)
    ax1.set_xlim(lo, hi)
    ax1.set_ylim(lo, hi)
    ax1.set_xlabel("ANIm ANI (%, truth)")
    ax1.set_ylabel("Estimated ANI (%)")
    ax1.grid(True)
    ax1.text(0.03, 0.97, f"all {n_flagged}/{len(ev)} pairs flagged INCONSISTENT",
             transform=ax1.transAxes, va="top", fontsize=6.5, style="italic")
    ax1.legend(loc="upper left", bbox_to_anchor=(0.02, 0.88))
    panel_label(ax1, "a")

    # (b) per-pair error dumbbell, sorted by truth ---------------------------
    ev = ev.sort_values("anim_ani").reset_index(drop=True)
    y = np.arange(len(ev))
    ax2.axvline(0, color=C_TRUTH, lw=0.8, ls="--", zorder=1)
    for i in y:
        errs = [ev.loc[i, col] - ev.loc[i, "anim_ani"] for col, *_ in F4_METHODS]
        ax2.plot([min(errs), max(errs)], [i, i], color=C_GREY, lw=0.5, zorder=2)
    for col, color, marker, label in F4_METHODS:
        ax2.scatter(ev[col] - ev["anim_ani"], y, color=color, marker=marker,
                    label=label, zorder=3, s=14)
    ax2.set_yticks(y)
    ax2.set_yticklabels([f"{a:.1f}" for a in ev["anim_ani"]], fontsize=5.5)
    ax2.set_ylabel("Pair (ANIm truth, %) →")
    ax2.set_xlabel("Error (est − ANIm truth, ANI points)")
    ax2.grid(True, axis="x")
    ax2.legend(loc="lower left", fontsize=6)
    # Metrics summary goes in panel (a)'s empty lower-left triangle (below the
    # identity line); panel (b) has no free vertical band for it.
    lines = []
    for col, _, _, label in F4_METHODS:
        m = met.loc[col]
        lines.append(f"{label}: MAE {m['MAE']:.2f}, bias {m['bias']:+.2f}, r {m['r']:.2f}")
    ax1.text(0.03, 0.03, "\n".join(lines), transform=ax1.transAxes, va="bottom",
             ha="left", fontsize=5.8,
             bbox=dict(fc="white", ec="none", alpha=0.8))
    panel_label(ax2, "b")

    fig.tight_layout()
    save(fig, "fig4_midani_anim_validation")


# --------------------------------------------------------------------------
# F5 — GTDB R207 large-scale benchmark vs FastANI
# --------------------------------------------------------------------------

def fig5():
    print("F5: GTDB R207 benchmark")
    g = pd.read_csv(F5_MATRIX, sep="\t")
    sub = g[g["fastani_ani"].notna()].copy()
    # skani/FastANI are stored as fractions in this matrix; syn2bani as percent.
    sub["fastani_pct"] = sub["fastani_ani"] * 100.0
    sub["skani_pct"] = sub["skani_ani"] * 100.0
    n = len(sub)
    print(f"  pairs with FastANI report: {n} (of {len(g)}); flags: "
          + ", ".join(f"{k}={v}" for k, v in sub["s2b_flag"].value_counts().items()))

    # ---- overall metrics (printed + written to TSV) ----
    methods = [("s2b_ani_uniform", L_UNIFORM), ("s2b_ani", L_GAMMA), ("skani_pct", L_SKANI)]
    rows = []
    for col, label in methods:
        d = sub.dropna(subset=[col])
        err = d[col] - d["fastani_pct"]
        rows.append({
            "method": label, "n": len(d),
            "MAE": err.abs().mean(), "RMSE": float(np.sqrt((err ** 2).mean())),
            "bias": err.mean(), "r": d[col].corr(d["fastani_pct"]),
        })
    summary = pd.DataFrame(rows)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(F5_SUMMARY_OUT, sep="\t", index=False, float_format="%.4f")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"  wrote {F5_SUMMARY_OUT.name}")

    fig, axes = plt.subplots(1, 3, figsize=(7, 2.9))

    # (a) density + flag-colored overlay ------------------------------------
    ax = axes[0]
    lo, hi = 80, 100
    hb = ax.hexbin(sub["fastani_pct"], sub["s2b_ani_uniform"], gridsize=30,
                   extent=(lo, hi, lo, hi), bins="log", cmap="Greys", mincnt=1,
                   zorder=2, linewidths=0.1)
    fig.colorbar(hb, ax=ax, label="pairs (log)", pad=0.02).ax.tick_params(labelsize=6)
    identity_line(ax, lo, hi)
    flag_style = {"ok": (C_UNIFORM, "o"), "INCONSISTENT": (C_GAMMA, "s"),
                  "BELOW_DETECTION": (C_GREY, "x")}
    for flag, (color, marker) in flag_style.items():
        d = sub[sub["s2b_flag"] == flag]
        ax.scatter(d["fastani_pct"], d["s2b_ani_uniform"], color=color, marker=marker,
                   s=6, alpha=0.7, label=f"{flag} ({len(d)})", zorder=3, linewidths=0.5)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("FastANI ANI (%)")
    ax.set_ylabel("syn2bani uniform ANI (%)")
    ax.legend(loc="lower right", fontsize=5.5, markerscale=1.5)
    panel_label(ax, "a")

    # (b, c) MAE and bias binned by FastANI ANI ------------------------------
    bins = [(76, 80), (80, 85), (85, 90), (90, 95), (95, 100)]
    labels = [f"{a}–{b}" for a, b in bins]
    centers = np.arange(len(bins))
    width = 0.26
    mcolors = [C_UNIFORM, C_GAMMA, C_SKANI]
    for ax, metric, letter in ((axes[1], "MAE", "b"), (axes[2], "bias", "c")):
        for j, (col, label) in enumerate(methods):
            vals = []
            for a, b in bins:
                d = sub[(sub["fastani_pct"] >= a) & (sub["fastani_pct"] < b)]
                err = d[col] - d["fastani_pct"]
                vals.append(err.abs().mean() if metric == "MAE" else err.mean())
            ax.bar(centers + (j - 1) * width, vals, width, color=mcolors[j],
                   label=label, edgecolor="none")
        if metric == "bias":
            ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--")
        ax.set_xticks(centers)
        ax.set_xticklabels(labels, fontsize=6)
        ax.set_xlabel("FastANI ANI bin (%)")
        ax.set_ylabel(f"{metric} vs FastANI (ANI points)")
        ax.grid(True, axis="y")
        panel_label(ax, letter)
    axes[1].legend(loc="upper left", fontsize=6)
    axes[2].legend(loc="lower left", fontsize=6)

    fig.tight_layout()
    save(fig, "fig5_gtdb_r207_benchmark")


# --------------------------------------------------------------------------
# F6 — computational efficiency
# --------------------------------------------------------------------------

# (tool, mode) -> (color, marker, label). syn2bani in two blues, comparators in
# their usual Okabe-Ito colors.
F6_SERIES = [
    ("syn2bani", "ani_fasta", C_UNIFORM, "o", "syn2bani (FASTA input)"),
    ("syn2bani", "ani_sketches", "#56B4E9", "^", "syn2bani (sketch reuse)"),
    ("skani", "dist", C_SKANI, "s", "skani dist"),
    ("fastani", "all_vs_all", C_FASTANI, "D", "FastANI"),
]


def fig6():
    print("F6: efficiency")
    rt = pd.read_csv(F6_RUNTIME, sep="\t")
    sk = pd.read_csv(F6_SKETCH, sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(7, 2.9))

    def agg(df, val):
        g = df.groupby(["n_pairs"])[val]
        med = g.median()
        return (med.index.values, med.values,
                (med - g.min()).values, (g.max() - med).values)

    # (a) wall time vs n_pairs, log-log; (b) peak RSS vs n_pairs, log y.
    # Median over 3 reps with min–max range as error bars.
    for ax, val, ylabel, letter in (
            (axes[0], "wall_s", "Wall time (s)", "a"),
            (axes[1], "peak_rss_mb", "Peak RSS (MB)", "b")):
        for tool, mode, color, marker, label in F6_SERIES:
            d = rt[(rt["tool"] == tool) & (rt["mode"] == mode)]
            x, med, lo, hi = agg(d, val)
            ax.errorbar(x, med, yerr=[lo, hi], color=color, marker=marker,
                        label=label, capsize=2, elinewidth=0.6, ms=3.5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Number of pairs (all-vs-all, n²)")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both")
        panel_label(ax, letter)
    # Shared legend for panels (a)/(b) below the figure — no free corner
    # exists inside either log-log panel.
    handles, lbls = axes[0].get_legend_handles_labels()
    fig.legend(handles, lbls, loc="lower center", bbox_to_anchor=(0.38, -0.04),
               ncol=4, fontsize=6, columnspacing=1.0, handlelength=1.6)

    # (c) sketch total size vs n_genomes, log y -------------------------------
    ax = axes[2]
    for tool, color, marker, label in (
            ("syn2bani", C_UNIFORM, "o", "syn2bani (.s2ba)"),
            ("skani", C_SKANI, "s", "skani db")):
        d = sk[sk["tool"] == tool]
        g = d.groupby("n_genomes")["total_size_kb"]
        med = g.median()
        ax.errorbar(med.index.values, med.values,
                    yerr=[(med - g.min()).values, (g.max() - med).values],
                    color=color, marker=marker, label=label, capsize=2,
                    elinewidth=0.6, ms=3.5)
    s2b22 = sk[(sk["tool"] == "syn2bani") & (sk["n_genomes"] == 22)]["total_size_kb"].median()
    sk22 = sk[(sk["tool"] == "skani") & (sk["n_genomes"] == 22)]["total_size_kb"].median()
    ax.annotate(f"n = 22: {s2b22:,.0f} vs {sk22:,.0f} KB\n(~{sk22 / s2b22:.0f}× smaller)",
                xy=(22, s2b22), xytext=(0.97, 0.04), textcoords="axes fraction",
                fontsize=6, ha="right", va="bottom",
                arrowprops=dict(arrowstyle="-", lw=0.6, color=C_GREY))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([2, 5, 10, 15, 22])
    ax.set_xticklabels([2, 5, 10, 15, 22])
    ax.minorticks_off()
    ax.set_xlabel("Number of genomes")
    ax.set_ylabel("Sketch total size (KB)")
    ax.grid(True, which="both")
    ax.legend(loc="upper left", fontsize=6)
    panel_label(ax, "c")

    fig.tight_layout()
    save(fig, "fig6_efficiency")


# --------------------------------------------------------------------------
# F7 — ANIm-truth benchmark by ANI band (2,074 GTDB R207 pairs)
# --------------------------------------------------------------------------

# method -> (color, marker, label). FastANI is plotted faint/dashed because it
# covers only a 363-pair subset of the 2,074 ANIm-truth pairs.
F7_METHODS = [
    ("syn2bani_4e_gamma", C_GAMMA, "o", "syn2bani 4e (gamma)"),
    ("syn2bani_4e_ridge_cv", C_UNIFORM, "^", "syn2bani 4e (ridge CV)"),
    ("syn2bani_11e", "#E69F00", "s", "syn2bani 11e (old panel)"),
    ("skani", C_SKANI, "D", L_SKANI),
]
F7_BANDS = ["0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-0.99"]
F7_BAND_LABELS = ["80–85", "85–90", "90–95", "95–99"]
F7_BAND_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]


def fig7_anim_by_band():
    print("F7: ANIm benchmark by band")
    tab = pd.read_csv(F7_TABLE, sep="\t")
    tab = tab[tab["band"] != "all"]
    x = np.arange(len(F7_BANDS))

    fig, axes = plt.subplots(1, 3, figsize=(7, 2.9))

    def series(method, col):
        d = tab[tab["method"] == method].set_index("band")
        return np.array([d.loc[b, col] if b in d.index else np.nan
                         for b in F7_BANDS])

    # (a) MAE vs band --------------------------------------------------------
    ax = axes[0]
    for method, color, marker, label in F7_METHODS:
        ax.plot(x, series(method, "MAE"), color=color, marker=marker, label=label)
    fa = series("FastANI_subset", "MAE")
    ax.plot(x, fa, color=C_FASTANI, marker="D", ls="--", alpha=0.45,
            label="FastANI (subset)")
    ax.annotate("FastANI: 363-pair\nsubset only", xy=(0.97, 0.62),
                xycoords="axes fraction", ha="right", va="top",
                fontsize=5.5, color=C_FASTANI, style="italic")
    ax.set_ylim(0, 4.6)
    ax.set_xticks(x)
    ax.set_xticklabels(F7_BAND_LABELS, fontsize=6.5)
    ax.set_xlabel("ANIm band (%)")
    ax.set_ylabel("MAE vs ANIm (ANI points)")
    ax.grid(True, axis="y")
    ax.legend(loc="upper left", fontsize=5.8)
    panel_label(ax, "a")

    # (b) signed bias vs band -------------------------------------------------
    ax = axes[1]
    ax.axhline(0, color=C_TRUTH, lw=0.8, ls="--")
    for method, color, marker, label in F7_METHODS:
        ax.plot(x, series(method, "bias"), color=color, marker=marker, label=label)
    ax.plot(x, series("FastANI_subset", "bias"), color=C_FASTANI, marker="D",
            ls="--", alpha=0.45, label="FastANI (subset)")
    ax.set_xticks(x)
    ax.set_xticklabels(F7_BAND_LABELS, fontsize=6.5)
    ax.set_xlabel("ANIm band (%)")
    ax.set_ylabel("Bias (est − ANIm, ANI points)")
    ax.grid(True, axis="y")
    ax.legend(loc="upper right", fontsize=5.8)
    panel_label(ax, "b")

    # (c) scatter ridge-CV vs ANIm truth --------------------------------------
    ax = axes[2]
    d = pd.read_csv(F7_PREDS, sep="\t")
    lo, hi = 82, 99
    identity_line(ax, lo, hi)
    for band, color in zip(F7_BANDS, F7_BAND_COLORS):
        sub = d[d["band"] == band]
        ax.scatter(sub["anim_ani"], sub["ridge_pred"], color=color, s=4,
                   alpha=0.6, linewidths=0, label=f"{band} (n={len(sub)})",
                   zorder=3)
    mae = (d["ridge_pred"] - d["anim_ani"]).abs().mean()
    r = d["ridge_pred"].corr(d["anim_ani"])
    ax.text(0.03, 0.97, f"MAE = {mae:.3f}\nr = {r:.3f}\nn = {len(d)}",
            transform=ax.transAxes, va="top", fontsize=6.5)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("ANIm ANI (%, truth)")
    ax.set_ylabel("syn2bani 4e ridge-CV ANI (%)")
    ax.grid(True)
    ax.legend(loc="lower right", fontsize=5.8, markerscale=2)
    panel_label(ax, "c")

    fig.tight_layout()
    save(fig, "fig7_anim_by_band")


# --------------------------------------------------------------------------

def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    fig7_anim_by_band()
    print("done.")


if __name__ == "__main__":
    main()
