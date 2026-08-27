#!/usr/bin/env python3
"""Publication-quality analysis of Syntracker validation cases.

Focus: pairs with near-clonal ANI (>99%) but low anchor adjacency, showing that
high ANI does not imply conserved genome architecture. Outputs figures to
figures/syntracker_validation/ and a case table to results/syntracker_validation/.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.patheffects as pe


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "syntracker_validation"
OUT_FIG = ROOT / "figures" / "syntracker_validation"
OUT_RES = ROOT / "results" / "syntracker_validation"


def norm_pair(a, b):
    return tuple(sorted([str(a), str(b)]))


def load_syn2bani(species):
    path = DATA / "syn2bani" / f"syn2bani_{species}.tsv"
    df = pd.read_csv(path, sep="\t")
    df["query"] = df["query"].astype(str)
    df["reference"] = df["reference"].astype(str)
    df = df[df["query"] != df["reference"]].copy()
    df["pair"] = [norm_pair(q, r) for q, r in zip(df["query"], df["reference"])]
    # Keep one row per unordered pair; use the row with the lower query name for determinism.
    df = df.sort_values(["pair", "query"]).drop_duplicates(subset="pair", keep="first")
    return df[["pair", "ani", "anchor_adjacency", "breakpoint_count", "af_query", "af_reference"]].copy()


def load_skani(species):
    path = DATA / "skani" / f"skani_{species}.tsv"
    df = pd.read_csv(path, sep="\t")
    # normalize column names
    colmap = {}
    for c in df.columns:
        low = c.lower().replace(" ", "_")
        if "query" in low:
            if "query_path" not in colmap.values() or "file" in low:
                colmap[c] = "query_path"
        elif ("reference" in low or "ref" in low) and "ani" not in low and "align" not in low:
            if "ref_path" not in colmap.values() or "file" in low:
                colmap[c] = "ref_path"
        elif low == "ani":
            colmap[c] = "ani_skani"
    df = df.rename(columns=colmap)
    df["query"] = df["query_path"].apply(lambda x: Path(x).stem)
    df["reference"] = df["ref_path"].apply(lambda x: Path(x).stem)
    df["pair"] = [norm_pair(q, r) for q, r in zip(df["query"], df["reference"])]
    return df[["pair", "ani_skani"]].drop_duplicates(subset="pair").copy()


def load_metadata(species):
    meta_path = DATA / "samples" / f"samples_{species}.tsv"
    if not meta_path.exists():
        return None
    meta = pd.read_csv(meta_path, sep="\t")
    meta["isolate"] = meta["isolate"].astype(str)
    return meta


def annotate_ecoli(row, info):
    q, r = row["pair"]
    q_info = info.get(q, "?")
    r_info = info.get(r, "?")
    return f"{q_info} vs {r_info}"


def annotate_hp(row, host_map):
    q, r = row["pair"]
    q_host = host_map.get(q, "?")
    r_host = host_map.get(r, "?")
    return f"host {q_host} vs {r_host}"


def main():
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    OUT_RES.mkdir(parents=True, exist_ok=True)

    species_list = ["Escherichia_coli_hypermutator", "Helicobacter_pylori",
                    "Neisseria_gonorrhoeae", "Streptomyces_rimosus"]

    all_cases = []
    species_frames = {}

    for species in species_list:
        syn = load_syn2bani(species)
        ska = load_skani(species)
        merged = syn.merge(ska, on="pair", how="outer")
        merged["species"] = species

        meta = load_metadata(species)
        if species == "Escherichia_coli_hypermutator":
            info = dict(zip(meta["isolate"].astype(str),
                            meta["mouse"].astype(str) + "_" + meta["day"].astype(str)))
            merged["q_info"] = merged["pair"].apply(lambda p: info.get(p[0], "?"))
            merged["r_info"] = merged["pair"].apply(lambda p: info.get(p[1], "?"))
            merged["same_mouse"] = merged["q_info"].str.split("_").str[0] == merged["r_info"].str.split("_").str[0]
            merged["time_delta"] = (merged["q_info"].str.split("_").str[1].astype(str).str.replace("d", "", regex=False).astype(int) -
                                    merged["r_info"].str.split("_").str[1].astype(str).str.replace("d", "", regex=False).astype(int)).abs()
        elif species == "Helicobacter_pylori":
            host_map = dict(zip(meta["isolate"].astype(str), meta["host"].astype(str)))
            merged["q_host"] = merged["pair"].apply(lambda p: host_map.get(p[0], "?"))
            merged["r_host"] = merged["pair"].apply(lambda p: host_map.get(p[1], "?"))
            merged["same_host"] = merged["q_host"] == merged["r_host"]

        # Focus on high-ANI pairs
        high = merged[merged["ani"] > 99.0].copy()
        if len(high) == 0:
            continue
        high["ani_rank"] = high["ani"].rank(pct=True)
        high["syn_rank"] = high["anchor_adjacency"].rank(pct=True)
        high["discordance"] = high["ani_rank"] - high["syn_rank"]

        top = high.nlargest(20, "discordance").copy()
        top["rank"] = range(1, len(top) + 1)
        all_cases.append(top)
        species_frames[species] = high

    cases_df = pd.concat(all_cases, ignore_index=True)
    cases_df.to_csv(OUT_RES / "top_discordant_cases.tsv", sep="\t", index=False, float_format="%.4f")

    # --- Figure 1: main 2x2 panel (Syn2bANI and skani ANI vs synteny) ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    def plot_one_panel(ax, df, x_col, title, ylim, color_by, annotate_n=2):
        if color_by == "time":
            sc = ax.scatter(df[x_col], df["anchor_adjacency"], c=df["time_delta"].astype(float),
                            cmap="viridis_r", s=35, alpha=0.75, edgecolors="none")
            plt.colorbar(sc, ax=ax, label="Days between isolates")
        elif color_by == "host":
            colors = df["same_host"].map({True: "#1f77b4", False: "#ff7f0e"})
            ax.scatter(df[x_col], df["anchor_adjacency"], c=colors, s=35, alpha=0.75, edgecolors="none")
            ax.scatter([], [], c="#1f77b4", s=30, label="Same host")
            ax.scatter([], [], c="#ff7f0e", s=30, label="Different host")
            ax.legend(loc="lower left", fontsize=8)
        xlabel = "Syn2bANI ANI (%)" if x_col == "ani" else "skani ANI (%)"
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Anchor adjacency")
        ax.set_title(title)
        ax.axhline(0.955, color="red", ls="--", lw=0.8, alpha=0.7, label="SynTracker cutoff")
        ax.set_xlim(99.0, 100.01)
        ax.set_ylim(ylim)
        # annotate top 1-2 discordant cases (defined by Syn2bANI discordance)
        top = cases_df[cases_df["species"] == df.name].head(annotate_n)
        offsets = [(25, 15), (-25, -15)]
        for i, (_, r) in enumerate(top.iterrows()):
            ox, oy = offsets[i]
            if "q_host" in r and pd.notna(r["q_host"]):
                label = f"{r['pair'][0]}({r['q_host']})\nvs {r['pair'][1]}({r['r_host']})"
            else:
                label = f"{r['pair'][0]}\nvs {r['pair'][1]}"
            ax.annotate(label, xy=(r[x_col], r["anchor_adjacency"]),
                        xytext=(ox, oy), textcoords="offset points",
                        fontsize=6, ha="center", va="center",
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.4,
                                       connectionstyle="arc3,rad=0.1"),
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="gray", alpha=0.85))

    eco = species_frames["Escherichia_coli_hypermutator"]
    eco.name = "Escherichia_coli_hypermutator"
    plot_one_panel(axes[0, 0], eco, "ani", "E. coli hypermutator — Syn2bANI", (0.70, 1.0), "time")
    plot_one_panel(axes[0, 1], eco, "ani_skani", "E. coli hypermutator — skani", (0.70, 1.0), "time")

    hp = species_frames["Helicobacter_pylori"]
    hp.name = "Helicobacter_pylori"
    plot_one_panel(axes[1, 0], hp, "ani", "H. pylori — Syn2bANI", (0.88, 0.985), "host")
    plot_one_panel(axes[1, 1], hp, "ani_skani", "H. pylori — skani", (0.88, 0.985), "host")

    plt.tight_layout()
    fig.savefig(OUT_FIG / "syntracker_high_ani_low_synteny.png", dpi=300)
    fig.savefig(OUT_FIG / "syntracker_high_ani_low_synteny.pdf")
    plt.close(fig)

    # --- Figure 2: breakpoints vs ANI, with marginal histograms ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, species, title in zip(axes, species_list,
                                   ["E. coli hypermutator", "H. pylori"]):
        sub = species_frames[species]
        ax.scatter(sub["ani"], sub["breakpoint_count"], s=30, alpha=0.6, edgecolors="none")
        # annotate top 3 with alternating offsets
        top = cases_df[cases_df["species"] == species].head(3)
        offsets = [(-40, 10), (40, -10), (-40, -10)]
        for i, (_, r) in enumerate(top.iterrows()):
            ox, oy = offsets[i]
            label = f"{r['pair'][0]} vs {r['pair'][1]}"
            ax.annotate(label, xy=(r["ani"], r["breakpoint_count"]),
                        xytext=(ox, oy), textcoords="offset points",
                        fontsize=7, ha="center", va="center",
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray", alpha=0.85))
        ax.set_xlabel("Syn2bANI ANI (%)")
        ax.set_ylabel("Breakpoint count")
        ax.set_title(title)
        ax.set_xlim(99.0, 100.01)
    plt.tight_layout()
    fig.savefig(OUT_FIG / "syntracker_breakpoints_vs_ani.png", dpi=300)
    fig.savefig(OUT_FIG / "syntracker_breakpoints_vs_ani.pdf")
    plt.close(fig)

    # --- Summary statistics ---
    summary_rows = []
    for species, df in species_frames.items():
        high = df[df["ani"] > 99.0]
        summary_rows.append({
            "species": species,
            "n_pairs": len(df),
            "n_high_ani": len(high),
            "mean_ani_high": high["ani"].mean(),
            "mean_synteny_high": high["anchor_adjacency"].mean(),
            "min_synteny_high": high["anchor_adjacency"].min(),
            "mean_breakpoints_high": high["breakpoint_count"].mean(),
            "max_breakpoints_high": high["breakpoint_count"].max(),
            "rho_ani_synteny": df["ani"].corr(df["anchor_adjacency"], method="spearman"),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_RES / "syntracker_summary.tsv", sep="\t", index=False, float_format="%.4f")

    # --- Supplementary figure: N. gonorrhoeae and S. rimosus ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    extra_species = ["Neisseria_gonorrhoeae", "Streptomyces_rimosus"]
    ylims = {"Neisseria_gonorrhoeae": (0.87, 0.91), "Streptomyces_rimosus": (0.90, 0.965)}
    for row, species in enumerate(extra_species):
        df = species_frames[species]
        for col, x_col in enumerate(["ani", "ani_skani"]):
            ax = axes[row, col]
            ax.scatter(df[x_col], df["anchor_adjacency"], s=35, alpha=0.75,
                       edgecolors="none", c="#2ca02c")
            xlabel = "Syn2bANI ANI (%)" if x_col == "ani" else "skani ANI (%)"
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Anchor adjacency")
            title = f"{species.replace('_', ' ')} — {'Syn2bANI' if x_col == 'ani' else 'skani'}"
            ax.set_title(title)
            ax.axhline(0.955, color="red", ls="--", lw=0.8, alpha=0.7)
            ax.set_xlim(99.0, 100.01)
            ax.set_ylim(ylims[species])
    plt.tight_layout()
    fig.savefig(OUT_FIG / "syntracker_supp_ngonorrhoeae_srimosus.png", dpi=300)
    fig.savefig(OUT_FIG / "syntracker_supp_ngonorrhoeae_srimosus.pdf")
    plt.close(fig)

    print(f"Wrote {OUT_FIG / 'syntracker_high_ani_low_synteny.png'}")
    print(f"Wrote {OUT_FIG / 'syntracker_breakpoints_vs_ani.png'}")
    print(f"Wrote {OUT_FIG / 'syntracker_supp_ngonorrhoeae_srimosus.png'}")
    print(f"Wrote {OUT_RES / 'top_discordant_cases.tsv'}")
    print(f"Wrote {OUT_RES / 'syntracker_summary.tsv'}")
    print("\nTop 5 discordant cases overall:")
    print(cases_df.head(5)[["species", "pair", "ani", "anchor_adjacency", "breakpoint_count"]].to_string(index=False))


if __name__ == "__main__":
    main()
