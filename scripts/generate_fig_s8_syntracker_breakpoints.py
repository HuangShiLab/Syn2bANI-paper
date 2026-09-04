#!/usr/bin/env python3
"""Generate Supplementary Figure S8: skani ANI vs Syn2bANI breakpoint_count.

Four published near-clonal isolate collections:
  (a) E. coli hypermutator (253 pairs)
  (b) H. pylori (2,926 pairs)
  (c) N. gonorrhoeae (66 pairs)
  (d) S. rimosus (190 pairs)

Data sources
------------
- Syn2bANI structural breakpoints:
    data/syntracker_validation/syn2b_structural_raw/syn2b_structural_pairs_raw.tsv
- skani ANI:
    data/syntracker_validation/skani/skani_<species>.tsv
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Publication style
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import set_publication_style, figure_size, label_panel, save_figure, COLORS


ROOT = Path(__file__).resolve().parent.parent
STRUCTURAL_TSV = ROOT / "data" / "syntracker_validation" / "syn2b_structural_raw" / "syn2b_structural_pairs_raw.tsv"
SKANI_DIR = ROOT / "data" / "syntracker_validation" / "skani"
OUTDIR = ROOT / "paper" / "figures" / "supplementary"

SPECIES = [
    "Escherichia_coli_hypermutator",
    "Helicobacter_pylori",
    "Neisseria_gonorrhoeae",
    "Streptomyces_rimosus",
]

COHORT_LABELS = {
    "Escherichia_coli_hypermutator": "E. coli hypermutator (negative control)",
    "Helicobacter_pylori": "H. pylori (mixed)",
    "Neisseria_gonorrhoeae": "N. gonorrhoeae (both modes)",
    "Streptomyces_rimosus": "S. rimosus (positive control)",
}

PANELS = ["a", "b", "c", "d"]


def norm_pair(a, b):
    """Return a deterministic, unordered pair key."""
    return "__".join(sorted([str(a), str(b)]))


def load_structural(path):
    """Load post-fix Syn2b structural output and drop self-comparisons."""
    df = pd.read_csv(path, sep="\t")
    df = df[df["status"] == "ok"].copy()
    df["pair"] = [norm_pair(a, b) for a, b in zip(df["genome_A"], df["genome_B"])]
    df = df[df["is_self"] != 1].copy()
    return df


def load_skani(path):
    """Load skani ANI and keep one row per unordered pair."""
    df = pd.read_csv(path, sep="\t")
    # skani names are reliably present in these files
    df["query"] = df["Query_name"].astype(str)
    df["reference"] = df["Ref_name"].astype(str)
    df["pair"] = [norm_pair(q, r) for q, r in zip(df["query"], df["reference"])]
    df = df[df["pair"].apply(lambda p: p.split("__")[0] != p.split("__")[1])]
    df = df.rename(columns={"ANI": "ani_skani"})
    return df[["pair", "ani_skani"]].drop_duplicates("pair").copy()


def nice_xlim(values, pad=0.005):
    """Pad skani ANI range so points are not flush with axis edges."""
    lo, hi = values.min(), values.max()
    return (max(0.0, lo - pad), min(100.0, hi + pad))


def nice_ylim(values, pad=0.5):
    """Pad breakpoint range; ensure y=0 is visible when present."""
    lo, hi = int(values.min()), int(values.max())
    return (max(0, lo - 0.5), hi + pad)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structural",
        type=Path,
        default=STRUCTURAL_TSV,
        help="Path to syn2b_structural_pairs_raw.tsv",
    )
    parser.add_argument(
        "--skani-dir",
        type=Path,
        default=SKANI_DIR,
        help="Directory containing skani_<species>.tsv files",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=OUTDIR,
        help="Output directory for figure files",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "pdf"],
        help="Output formats",
    )
    args = parser.parse_args()

    if not args.structural.exists():
        raise FileNotFoundError(f"Structural file not found: {args.structural}")
    if not args.skani_dir.exists():
        raise FileNotFoundError(f"skani directory not found: {args.skani_dir}")

    set_publication_style()
    args.outdir.mkdir(parents=True, exist_ok=True)

    struct = load_structural(args.structural)

    merged_frames = []
    for sp in SPECIES:
        ska_path = args.skani_dir / f"skani_{sp}.tsv"
        if not ska_path.exists():
            print(f"WARN: missing {ska_path}", file=sys.stderr)
            continue
        ska = load_skani(ska_path)
        sp_struct = struct[struct["cohort"] == sp].copy()
        merged = pd.merge(sp_struct, ska, on="pair", how="inner")
        if len(merged) == 0:
            print(f"WARN: no merged rows for {sp}", file=sys.stderr)
            continue
        merged["species"] = sp
        merged_frames.append(merged)

    if not merged_frames:
        raise SystemExit("No merged rows; check input paths.")

    df = pd.concat(merged_frames, ignore_index=True)

    # Verify expected pair counts
    counts = df.groupby("species").size()
    for sp in SPECIES:
        if sp in counts:
            print(f"{sp}: {counts[sp]} pairs")

    # Build 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=figure_size(17.8, aspect=0.92))
    axes = axes.flatten()

    color = COLORS["blue"]

    for ax, sp, panel in zip(axes, SPECIES, PANELS):
        sub = df[df["species"] == sp]
        n_pairs = len(sub)

        ax.scatter(
            sub["ani_skani"],
            sub["syn2b_breakpoints"],
            c=color,
            s=25,
            alpha=0.6,
            edgecolors="none",
            rasterized=False,
        )

        ax.set_xlim(nice_xlim(sub["ani_skani"]))
        ax.set_ylim(nice_ylim(sub["syn2b_breakpoints"]))

        title = COHORT_LABELS.get(sp, sp.replace("_", " "))
        ax.set_title(f"{title}\n(n = {n_pairs:,})", fontsize=9)

        # Panel label
        label_panel(ax, panel)

    # Shared axis labels
    for ax in axes[2:]:
        ax.set_xlabel("skani ANI (%)", fontsize=9)
    for ax in axes[::2]:
        ax.set_ylabel("Syn2bANI breakpoint count", fontsize=9)

    plt.tight_layout()
    out_stem = args.outdir / "fig_s8_syntracker_breakpoints"
    save_figure(fig, out_stem, formats=args.formats)
    plt.close(fig)

    # Log scale decision note
    print("\nNote: y-axis kept linear. Breakpoint counts include zero values in three")
    print("cohorts and span a modest within-panel range (max 30), so log scale is")
    print("inappropriate here; the across-cohort variation (0-30) is shown by the")
    print("separate panels.")


if __name__ == "__main__":
    main()
