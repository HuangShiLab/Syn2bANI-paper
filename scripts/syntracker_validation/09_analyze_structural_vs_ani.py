#!/usr/bin/env python3
"""Re-analyze SynTracker cohorts using the post-fix Syn2b structural channel.

Merges skani ANI with the `syn2b_structural_pairs_raw.tsv` output from
`08_syn2b_structural.py`, drops self-comparisons, and reproduces the four-cohort
control figures and correlation table.

This replaces `06_plot_results.py` for the structural-channel analysis; the old
script used the pre-fix `syn2bani/` outputs whose self-comparison floor was not
zero.
"""
import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


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


def norm_pair(a, b):
    return "__".join(sorted([str(a), str(b)]))


def load_skani(path):
    df = pd.read_csv(path, sep="\t")
    # prefer Ref_name/Query_name if present, else derive from file paths
    if "Ref_name" in df.columns and "Query_name" in df.columns:
        df["ref"] = df["Ref_name"]
        df["query"] = df["Query_name"]
    else:
        df["ref"] = df["Ref_file"].apply(lambda x: Path(x).stem)
        df["query"] = df["Query_file"].apply(lambda x: Path(x).stem)
    df["pair"] = [norm_pair(q, r) for q, r in zip(df["query"], df["ref"])]
    df = df.rename(columns={"ANI": "ani_skani"})
    df = df[["pair", "ani_skani"]].copy()
    df = df[df["pair"].apply(lambda p: p.split("__")[0] != p.split("__")[1])]
    # skani sometimes emits duplicate orientations
    df = df.drop_duplicates(subset="pair")
    return df


def load_structural(path):
    df = pd.read_csv(path, sep="\t")
    df = df[df["status"] == "ok"].copy()
    df["pair"] = [norm_pair(a, b) for a, b in zip(df["genome_A"], df["genome_B"])]
    df = df[df["is_self"] != 1].copy()
    return df


def add_host(df, meta_dir):
    meta = pd.read_csv(Path(meta_dir) / "samples_Helicobacter_pylori.tsv", sep="\t")
    iso_to_host = dict(zip(meta["isolate"].astype(str), meta["host"].astype(str)))
    df["host"] = df["pair"].apply(
        lambda p: iso_to_host.get(p.split("__")[0])
        or iso_to_host.get(p.split("__")[1], "unknown")
    )
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structural", required=True)
    p.add_argument("--skani-dir", required=True)
    p.add_argument("--metadata-dir", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    struct = load_structural(args.structural)

    frames = []
    for sp in SPECIES:
        ska_path = Path(args.skani_dir) / f"skani_{sp}.tsv"
        if not ska_path.exists():
            print(f"WARN: missing {ska_path}", file=sys.stderr)
            continue
        ska = load_skani(ska_path)
        sp_struct = struct[struct["cohort"] == sp].copy()
        merged = pd.merge(sp_struct, ska, on="pair", how="inner")
        merged["species"] = sp
        if sp == "Helicobacter_pylori":
            merged = add_host(merged, args.metadata_dir)
        frames.append(merged)

    df = pd.concat(frames, ignore_index=True)
    if len(df) == 0:
        raise SystemExit("No merged rows; check input paths.")

    # Save merged table
    df.to_csv(outdir / "merged_ani_synteny.tsv", sep="\t", index=False)

    # Correlation summary
    def _corr(g, x, y, method="spearman"):
        if g[x].notna().sum() < 3 or g[y].notna().sum() < 3:
            return np.nan
        return g[x].corr(g[y], method=method)

    corr = (
        df.groupby("species")
        .apply(
            lambda g: pd.Series({
                "n_pairs": len(g),
                "rho_ani_breakpoints": _corr(g, "ani_skani", "syn2b_breakpoints"),
                "rho_ani_inverted": _corr(g, "ani_skani", "syn2b_raw_inverted_fraction"),
                "rho_ani_scj": _corr(g, "ani_skani", "syn2b_scj_distance"),
                "mean_ani": g["ani_skani"].mean(),
                "median_breakpoints": g["syn2b_breakpoints"].median(),
                "median_inverted": g["syn2b_raw_inverted_fraction"].median(),
                "median_scj": g["syn2b_scj_distance"].median(),
            })
        )
        .reset_index()
    )
    corr["expectation"] = corr["species"].map({
        "Escherichia_coli_hypermutator": "SNP-driven: low structural signal",
        "Streptomyces_rimosus": "structural: high breakpoints",
        "Neisseria_gonorrhoeae": "both modes move together",
        "Helicobacter_pylori": "mixed, participant-dependent",
    })
    corr.to_csv(outdir / "correlation_summary.tsv", sep="\t", index=False)
    print(corr.to_string(index=False))

    # Plots
    sns.set_theme(style="whitegrid")
    df["species_label"] = df["species"].map(COHORT_LABELS)

    # Figure 1: ANI vs breakpoints
    g = sns.FacetGrid(df, col="species_label", col_wrap=2, sharex=False, sharey=False,
                      height=4, aspect=1.1)
    g.map_dataframe(sns.scatterplot, x="ani_skani", y="syn2b_breakpoints",
                    alpha=0.6, edgecolor=None, s=40)
    g.set_axis_labels("skani ANI (%)", "Syn2b breakpoints")
    g.set_titles(col_template="{col_name}")
    plt.tight_layout()
    out = outdir / "syntracker_ani_vs_breakpoints.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved {out}")

    # Figure 2: ANI vs raw inverted fraction
    g = sns.FacetGrid(df, col="species_label", col_wrap=2, sharex=False, sharey=False,
                      height=4, aspect=1.1)
    g.map_dataframe(sns.scatterplot, x="ani_skani", y="syn2b_raw_inverted_fraction",
                    alpha=0.6, edgecolor=None, s=40)
    g.set_axis_labels("skani ANI (%)", "Syn2b raw inverted fraction")
    g.set_titles(col_template="{col_name}")
    plt.tight_layout()
    out = outdir / "syntracker_ani_vs_inverted_fraction.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved {out}")

    # H. pylori: color by host
    if "host" in df.columns and df["host"].notna().any():
        hdf = df[df["species"] == "Helicobacter_pylori"].copy()
        plt.figure(figsize=(6, 5))
        sns.scatterplot(data=hdf, x="ani_skani", y="syn2b_breakpoints",
                        hue="host", alpha=0.7, edgecolor=None)
        plt.xlabel("skani ANI (%)")
        plt.ylabel("Syn2b breakpoints")
        plt.title("H. pylori: ANI vs breakpoints by host")
        plt.tight_layout()
        out = outdir / "syntracker_h_pylori_ani_vs_breakpoints_by_host.png"
        plt.savefig(out, dpi=300)
        plt.close()
        print(f"Saved {out}")


if __name__ == "__main__":
    main()
