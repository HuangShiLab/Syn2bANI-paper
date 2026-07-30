#!/usr/bin/env python3
"""Merge Syn2bANI + skani outputs and plot ANI vs synteny_score per species."""
import argparse, sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def norm_pair(a, b):
    return tuple(sorted([a, b]))


def load_syn2bani(path):
    df = pd.read_csv(path, sep='\t')
    # Use calibrated ANI when available, otherwise rate-heterogeneous ani
    ani_col = 'ani_cal' if 'ani_cal' in df.columns else 'ani'
    df = df.rename(columns={ani_col: 'ani_syn2bani'})
    # Syn2bANI output uses the first contig ID; after our normalisation this equals the isolate
    df['pair'] = [norm_pair(q, r) for q, r in zip(df['query'], df['reference'])]
    df = df[['pair', 'ani_syn2bani', 'synteny_score', 'breakpoint_count']].copy()
    # Drop self-comparisons
    df = df[df['pair'].apply(lambda x: x[0] != x[1])]
    return df


def load_skani(path):
    df = pd.read_csv(path, sep='\t')
    # skani column names vary slightly; normalise
    colmap = {}
    for c in df.columns:
        low = c.lower().replace(' ', '_')
        if 'query' in low:
            colmap[c] = 'query_path'
        elif 'reference' in low and 'ani' not in low:
            colmap[c] = 'ref_path'
        elif 'ani' in low:
            colmap[c] = 'ani_skani'
    df = df.rename(columns=colmap)
    df['query'] = df['query_path'].apply(lambda x: Path(x).stem)
    df['reference'] = df['ref_path'].apply(lambda x: Path(x).stem)
    df['pair'] = [norm_pair(q, r) for q, r in zip(df['query'], df['reference'])]
    df = df[['pair', 'ani_skani']].copy()
    df = df[df['pair'].apply(lambda x: x[0] != x[1])]
    return df


def annotate_species(row, species_order):
    for sp in species_order:
        if row['pair'][0].startswith(sp.split('_')[0][:2].upper() + '_') or \
           row['pair'][1].startswith(sp.split('_')[0][:2].upper() + '_'):
            return sp
    return 'unknown'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--syn2bani-dir', required=True)
    parser.add_argument('--skani-dir', required=True)
    parser.add_argument('--metadata-dir', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    species = [
        'Neisseria_gonorrhoeae',
        'Escherichia_coli_hypermutator',
        'Helicobacter_pylori',
        'Streptomyces_rimosus',
    ]

    all_frames = []
    for sp in species:
        syn_path = Path(args.syn2bani_dir) / f'syn2bani_{sp}.tsv'
        ska_path = Path(args.skani_dir) / f'skani_{sp}.tsv'
        if not syn_path.exists():
            print(f'WARN: missing {syn_path}', file=sys.stderr)
            continue
        syn = load_syn2bani(syn_path)
        if ska_path.exists():
            ska = load_skani(ska_path)
            merged = pd.merge(syn, ska, on='pair', how='outer')
        else:
            merged = syn.copy()
            merged['ani_skani'] = float('nan')
        merged['species'] = sp

        # Add host/participant info for H. pylori
        if sp == 'Helicobacter_pylori':
            meta = pd.read_csv(Path(args.metadata_dir) / f'samples_{sp}.tsv', sep='\t')
            iso_to_host = dict(zip(meta['isolate'], meta['host'].astype(str)))
            merged['host'] = merged['pair'].apply(lambda p: iso_to_host.get(p[0]) or iso_to_host.get(p[1], 'unknown'))

        all_frames.append(merged)

    df = pd.concat(all_frames, ignore_index=True)
    # Unpack pair into columns for the table
    df['query'] = df['pair'].apply(lambda x: x[0])
    df['reference'] = df['pair'].apply(lambda x: x[1])
    df.drop(columns=['pair']).to_csv(outdir / 'merged_ani_synteny.tsv', sep='\t', index=False)

    # Scatter: ANI (Syn2bANI calibrated) vs synteny_score, faceted by species
    sns.set_theme(style='whitegrid')
    g = sns.FacetGrid(df, col='species', col_wrap=2, sharex=False, sharey=False,
                      height=4, aspect=1.1)
    g.map_dataframe(sns.scatterplot, x='ani_syn2bani', y='synteny_score',
                    alpha=0.6, edgecolor=None, s=40)
    g.set_axis_labels('Syn2bANI ANI (%)', 'Syn2bANI synteny score')
    g.set_titles(col_template='{col_name}')
    for ax in g.axes.flat:
        ax.axhline(0.955, color='red', ls='--', lw=0.8, label='SynTracker cutoff')
    g.add_legend()
    plt.tight_layout()
    out_png = outdir / 'ani_vs_synteny_syntracker_species.png'
    plt.savefig(out_png, dpi=300)
    print(f'Saved {out_png}')

    # Spearman correlation per species
    corr = df.groupby('species').apply(
        lambda g: pd.Series({
            'n_pairs': len(g),
            'rho_ani_synteny': g['ani_syn2bani'].corr(g['synteny_score'], method='spearman'),
            'rho_skani_synteny': g['ani_skani'].corr(g['synteny_score'], method='spearman') if g['ani_skani'].notna().sum() > 2 else float('nan'),
            'mean_synteny': g['synteny_score'].mean(),
            'std_synteny': g['synteny_score'].std(),
            'mean_ani': g['ani_syn2bani'].mean(),
        })
    ).reset_index()
    corr.to_csv(outdir / 'correlation_summary.tsv', sep='\t', index=False)
    print(corr.to_string(index=False))


if __name__ == '__main__':
    main()
