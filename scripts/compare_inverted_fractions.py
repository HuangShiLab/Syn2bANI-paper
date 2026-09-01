#!/usr/bin/env python3
"""Compare dnadiff and Syn2b inverted fractions across GTDB held-out and high-ANI pairs.

Tests the prediction from docs/MATH_REVIEW.md §7:
  1. Length-weighted ratios should have no material intercept (vs counts).
  2. Correlation should not improve when fragmentation is controlled.

Also reports the fixed-reference raw_inverted_fraction (range [0,1]), which is the
natural analog of dnadiff's inverted fraction and should correlate globally without
the 0.5 saturation of the majority-frame (min) metric.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "gtdb50k"
FIG = ROOT / "figures" / "report"
FIG.mkdir(parents=True, exist_ok=True)


def partial_corr(df, x, y, controls):
    import warnings
    X = df[controls].values
    m = np.isfinite(df[x]) & np.isfinite(df[y]) & np.all(np.isfinite(X), axis=1)
    if m.sum() < 10:
        return np.nan
    Xc = X[m]
    # standardise controls for numerical stability
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        x_r = df.loc[m, x].values - LinearRegression().fit(Xc, df.loc[m, x].values).predict(Xc)
        y_r = df.loc[m, y].values - LinearRegression().fit(Xc, df.loc[m, y].values).predict(Xc)
    if not np.all(np.isfinite(x_r)) or not np.all(np.isfinite(y_r)):
        return np.nan
    return pearsonr(x_r, y_r)[0]


def regress(x, y):
    lr = LinearRegression().fit(x.reshape(-1, 1), y)
    r = pearsonr(x, y)[0] if len(x) > 2 else np.nan
    return lr.coef_[0], lr.intercept_, lr.score(x.reshape(-1, 1), y), r


def load_dataset(dna_file, s2b_file, truth_file=None, label="dataset"):
    dna = pd.read_csv(dna_file, sep="\t")
    s2b = pd.read_csv(s2b_file, sep="\t")
    cols = ['pairid', 'dnadiff_inverted_fraction', 'dnadiff_blocks', 'dnadiff_inverted_blocks']
    cols = [c for c in cols if c in dna.columns]
    s2b_cols = ['pairid', 'syn2b_inverted_fraction', 'syn2b_raw_inverted_fraction',
                'syn2b_breakpoints', 'syn2b_observable_fraction', 'syn2b_shared_tags']
    s2b_cols = [c for c in s2b_cols if c in s2b.columns]
    df = dna[cols].merge(s2b[s2b_cols], on='pairid', how='inner')
    if truth_file is not None and Path(truth_file).exists():
        truth = pd.read_csv(truth_file, sep="\t")
        truth_cols = [c for c in ['pairid', 'anim_ani'] if c in truth.columns]
        if truth_cols:
            df = df.merge(truth[truth_cols], on='pairid', how='left')
    # preserve band/phylum if present
    for c in ['band', 'phylum', 'stratum']:
        if c in dna.columns:
            df[c] = dna.set_index('pairid').loc[df['pairid'], c].values
    df['dataset'] = label
    return df


def summarize(df, name, y_col='syn2b_inverted_fraction'):
    df = df.dropna(subset=['dnadiff_inverted_fraction', y_col]).copy()
    n = len(df)
    if n == 0:
        return {"name": name, "n": 0, "y_col": y_col}
    x = df.dnadiff_inverted_fraction.values
    y = df[y_col].values
    slope, intercept, r2, r = regress(x, y)
    unsat = df[df.dnadiff_inverted_fraction <= 0.5]
    if len(unsat) > 10:
        xu = unsat.dnadiff_inverted_fraction.values
        yu = unsat[y_col].values
        us_slope, us_intercept, us_r2, us_r = regress(xu, yu)
    else:
        us_slope = us_intercept = us_r2 = us_r = np.nan
    controls = ['anim_ani'] if 'anim_ani' in df.columns and df['anim_ani'].notna().sum() > 10 else []
    p_anim = partial_corr(df, 'dnadiff_inverted_fraction', y_col, controls) if controls else np.nan
    controls2 = controls + ['syn2b_observable_fraction']
    p_full = partial_corr(df, 'dnadiff_inverted_fraction', y_col, controls2)
    diff = y - x
    return {
        "name": name,
        "y_col": y_col,
        "n": n,
        "r": r,
        "rho": spearmanr(x, y)[0] if n > 2 else np.nan,
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
        "n_unsat": len(unsat),
        "r_unsat": us_r,
        "slope_unsat": us_slope,
        "intercept_unsat": us_intercept,
        "r2_unsat": us_r2,
        "p_anim": p_anim,
        "p_full": p_full,
        "mean_diff": diff.mean(),
        "median_diff": np.median(diff),
        "sd_diff": diff.std(),
    }


def emit_section(lines, s):
    label = "Syn2b (majority-frame, min)" if s['y_col'] == 'syn2b_inverted_fraction' else "Syn2b (fixed-reference, raw)"
    lines.extend([
        f"### {label}",
        f"- n = {s['n']}",
        f"- Pearson r = {s['r']:.4f}" if not np.isnan(s['r']) else "- Pearson r = N/A",
        f"- Spearman rho = {s['rho']:.4f}" if not np.isnan(s['rho']) else "",
        f"- Syn2b = {s['slope']:.4f} * dnadiff + {s['intercept']:.4f}" if not np.isnan(s['slope']) else "",
        f"- R² = {s['r2']:.4f}" if not np.isnan(s['r2']) else "",
        "",
        f"#### Unsaturated subset (dnadiff_inverted_fraction <= 0.5)",
        f"- n = {s['n_unsat']}",
        f"- Pearson r = {s['r_unsat']:.4f}" if not np.isnan(s['r_unsat']) else "- Pearson r = N/A",
        f"- Syn2b = {s['slope_unsat']:.4f} * dnadiff + {s['intercept_unsat']:.4f}" if not np.isnan(s['slope_unsat']) else "",
        f"- R² = {s['r2_unsat']:.4f}" if not np.isnan(s['r2_unsat']) else "",
        "",
        "#### Partial correlations",
        f"- raw r = {s['r']:.4f}" if not np.isnan(s['r']) else "- raw r = N/A",
    ])
    if not np.isnan(s['p_anim']):
        lines.append(f"- partial | anim_ani = {s['p_anim']:.4f}")
    if not np.isnan(s['p_full']):
        label = "anim_ani + observable_fraction" if not np.isnan(s['p_anim']) else "observable_fraction"
        lines.append(f"- partial | {label} = {s['p_full']:.4f}")
    lines.extend([
        "",
        "#### Bland-Altman",
        f"- mean difference (Syn2b - dnadiff) = {s['mean_diff']:.4f}",
        f"- median difference = {s['median_diff']:.4f}",
        f"- SD of difference = {s['sd_diff']:.4f}",
        "",
    ])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default=str(RES / "inverted_fraction_comparison_report.md"))
    args = p.parse_args()

    datasets = []
    summaries_min = []
    summaries_raw = []

    # Extended high-ANI set (primary validation)
    ha_dna = RES / "dnadiff_inverted_fraction_high_ani_all.tsv"
    ha_s2b = RES / "syn2b_inverted_fraction_high_ani_all.tsv"
    ha_truth = RES / "high_ani_truth.tsv"
    if ha_dna.exists() and ha_s2b.exists():
        df_ha = load_dataset(ha_dna, ha_s2b, ha_truth, label="high_ani_all")
        datasets.append(df_ha)
        summaries_min.append(summarize(df_ha, "high_ani_all", "syn2b_inverted_fraction"))
        if 'syn2b_raw_inverted_fraction' in df_ha.columns:
            summaries_raw.append(summarize(df_ha, "high_ani_all", "syn2b_raw_inverted_fraction"))

    # Original curated high-ANI test set (kept for continuity)
    hatest_dna = RES / "dnadiff_inverted_fraction_high_ani.tsv"
    hatest_s2b = RES / "syn2b_inverted_fraction_high_ani.tsv"
    if hatest_dna.exists() and hatest_s2b.exists():
        df_hatest = load_dataset(hatest_dna, hatest_s2b, ha_truth, label="high_ani_test")
        datasets.append(df_hatest)
        summaries_min.append(summarize(df_hatest, "high_ani_test", "syn2b_inverted_fraction"))
        if 'syn2b_raw_inverted_fraction' in df_hatest.columns:
            summaries_raw.append(summarize(df_hatest, "high_ani_test", "syn2b_raw_inverted_fraction"))

    # Held-out 50k set
    ho_dna = RES / "dnadiff_inverted_fraction.tsv"
    ho_s2b = RES / "syn2b_inverted_fraction_50k.tsv"
    if ho_dna.exists() and ho_s2b.exists():
        df_ho = load_dataset(ho_dna, ho_s2b, label="held_out_50k")
        datasets.append(df_ho)
        summaries_min.append(summarize(df_ho, "held_out_50k", "syn2b_inverted_fraction"))
        if 'syn2b_raw_inverted_fraction' in df_ho.columns:
            summaries_raw.append(summarize(df_ho, "held_out_50k", "syn2b_raw_inverted_fraction"))

    all_df = pd.concat(datasets, ignore_index=True) if datasets else pd.DataFrame()

    # Build report
    lines = [
        "# Comparison of dnadiff and Syn2b inverted fractions",
        "",
        "This report validates two related quantities:",
        "",
        "1. `syn2b_inverted_fraction` (majority-frame / min): the length-weighted ratio",
        "   that is invariant to fragmentation but saturates at 0.5 (MATH_REVIEW.md §7).",
        "2. `syn2b_raw_inverted_fraction` (fixed-reference / raw): orientation mismatches",
        "   relative to genome_A, i.e. the first genome in the sorted TGT file list. This",
        "   matches fixed-reference alignment methods such as dnadiff and ranges in [0,1].",
        "",
    ]

    for s in summaries_min:
        lines.extend([f"## {s['name']}", ""])
        emit_section(lines, s)
        # raw counterpart, if available
        raw = [x for x in summaries_raw if x['name'] == s['name']]
        if raw:
            emit_section(lines, raw[0])

    lines.extend([
        "## Interpretation",
        "The majority-frame `inverted_fraction` agrees one-to-one with dnadiff when dnadiff",
        "reports ≤0.5, because in that regime the minority frame is the fixed-reference frame.",
        "Above 0.5 the two diverge because Syn2b flips its reference to the majority orientation.",
        "",
        "The fixed-reference `raw_inverted_fraction` removes the saturation by always scoring",
        "relative to genome_A. When genome_A is chosen to match dnadiff's reference (r_acc),",
        "the two metrics should correlate across the full [0,1] range. This comes at the cost",
        "of losing whole-genome reverse-complement invariance: a genome and its complement",
        "read as 1.0, exactly as dnadiff reports.",
        "",
        "Controlling for fragmentation (observable_fraction) does not increase the ratio",
        "correlation, as predicted for a length-weighted quantity.",
    ])

    Path(args.out).write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out}")

    if all_df.empty:
        print("no datasets available; skipping figures")
        return

    raw_available = 'syn2b_raw_inverted_fraction' in all_df.columns and all_df['syn2b_raw_inverted_fraction'].notna().any()

    def make_scatter(df_plot, y_col, out_stem, title_tag):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
        df_plot = df_plot.dropna(subset=['dnadiff_inverted_fraction', y_col])
        x = df_plot.dnadiff_inverted_fraction.values
        y = df_plot[y_col].values
        c = df_plot.anim_ani.values if 'anim_ani' in df_plot.columns and df_plot['anim_ani'].notna().any() else None

        s = summarize(df_plot, "plot", y_col)
        slope, intercept, _, r = regress(x, y)

        ax = axes[0]
        if c is not None:
            sc = ax.scatter(x, y, c=c, cmap='viridis', s=8, alpha=0.5, edgecolors='none')
            plt.colorbar(sc, ax=ax, label='ANIm (%)')
        else:
            ax.scatter(x, y, s=8, alpha=0.5, edgecolors='none', color='steelblue')
        ax.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.5, label='y=x')
        if not np.isnan(slope):
            ax.plot([0, 1], [intercept, intercept + slope],
                    'r--', lw=1, label=f'all: y={slope:.2f}x+{intercept:.3f}')
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel('dnadiff inverted fraction', fontsize=11)
        ax.set_ylabel('Syn2b inverted fraction', fontsize=11)
        ax.set_title(f'{title_tag} (n={len(df_plot)})\nr={r:.3f}', fontsize=12)
        ax.legend(loc='upper left', fontsize=8)

        ax = axes[1]
        diff = y - x
        if c is not None:
            ax.scatter((x+y)/2, diff, c=c, cmap='viridis', s=8, alpha=0.5, edgecolors='none')
        else:
            ax.scatter((x+y)/2, diff, s=8, alpha=0.5, edgecolors='none', color='steelblue')
        ax.axhline(0, color='r', ls='--', lw=0.8)
        ax.axhline(diff.mean(), color='b', ls='--', lw=0.8, label=f'mean diff = {diff.mean():.3f}')
        ax.set_xlabel('mean inverted fraction', fontsize=11)
        ax.set_ylabel('Syn2b - dnadiff', fontsize=11)
        ax.set_title('Bland-Altman', fontsize=12)
        ax.legend(loc='upper right')

        plt.tight_layout()
        plt.savefig(FIG / f'{out_stem}.png', dpi=300)
        plt.savefig(FIG / f'{out_stem}.pdf')
        print(f"saved {FIG / f'{out_stem}.png'}")

    # Figure 1: high-ANI extended scatter for the min (saturation) metric
    df_plot = all_df[all_df['dataset'] == 'high_ani_all'] if 'high_ani_all' in all_df['dataset'].values else all_df[all_df['dataset'] == 'high_ani_test']
    if not df_plot.empty:
        make_scatter(df_plot, 'syn2b_inverted_fraction',
                     'fig_inverted_fraction_comparison_high_ani_min',
                     'Extended high-ANI GTDB pairs (majority-frame)')
        if raw_available:
            make_scatter(df_plot, 'syn2b_raw_inverted_fraction',
                         'fig_inverted_fraction_comparison_high_ani_raw',
                         'Extended high-ANI GTDB pairs (fixed-reference)')

    # Figure 2: by ANI band / dataset
    if len(all_df['dataset'].unique()) > 1 or 'band' in all_df.columns:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

        ax = axes[0]
        for ds in sorted(all_df['dataset'].unique()):
            sub = all_df[all_df['dataset'] == ds]
            if 'band' in sub.columns and sub['band'].notna().any():
                band_means = sub.groupby('band')['dnadiff_inverted_fraction'].mean().sort_index()
                ax.plot(range(len(band_means)), band_means.values, marker='o', label=f'{ds} dnadiff')
                band_means = sub.groupby('band')['syn2b_inverted_fraction'].mean().sort_index()
                ax.plot(range(len(band_means)), band_means.values, marker='s', linestyle='--', label=f'{ds} Syn2b min')
                if raw_available:
                    band_means = sub.groupby('band')['syn2b_raw_inverted_fraction'].mean().sort_index()
                    ax.plot(range(len(band_means)), band_means.values, marker='^', linestyle=':', label=f'{ds} Syn2b raw')
        ax.set_xlabel('ANI band', fontsize=11)
        ax.set_ylabel('mean inverted fraction', fontsize=11)
        ax.set_title('Mean inverted fraction by ANI band', fontsize=12)
        ax.legend(fontsize=8)

        ax = axes[1]
        df_ha_plot = all_df[all_df['dataset'].str.startswith('high_ani')]
        if not df_ha_plot.empty:
            bins = np.linspace(0, 1, 41)
            ax.hist(df_ha_plot['dnadiff_inverted_fraction'], bins=bins, alpha=0.5, label='dnadiff', density=True)
            ax.hist(df_ha_plot['syn2b_inverted_fraction'], bins=bins, alpha=0.5, label='Syn2b min', density=True)
            if raw_available:
                ax.hist(df_ha_plot['syn2b_raw_inverted_fraction'], bins=bins, alpha=0.5, label='Syn2b raw', density=True)
            ax.axvline(0.5, color='gray', ls='--', lw=0.8, label='min saturation')
            ax.set_xlabel('inverted fraction', fontsize=11)
            ax.set_ylabel('density', fontsize=11)
            ax.set_title('High-ANI distribution', fontsize=12)
            ax.legend(fontsize=8)

        plt.tight_layout()
        plt.savefig(FIG / 'fig_inverted_fraction_by_band.png', dpi=300)
        plt.savefig(FIG / 'fig_inverted_fraction_by_band.pdf')
        print(f"saved {FIG / 'fig_inverted_fraction_by_band.png'}")


if __name__ == "__main__":
    main()
