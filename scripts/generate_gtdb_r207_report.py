#!/usr/bin/env python3
"""Generate GTDB-R207 summary statistics report."""
import argparse
from pathlib import Path

import pandas as pd
import numpy as np


def parse_taxonomy(tax_str: str) -> dict:
    """Parse GTDB taxonomy string into dict."""
    parts = [p.strip() for p in tax_str.split(';')]
    keys = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    return {k: parts[i] if i < len(parts) else 'unknown'
            for i, k in enumerate(keys)}


def load_metadata(meta_file: Path) -> pd.DataFrame:
    """Load GTDB metadata TSV."""
    cols = [
        'accession', 'checkm_completeness', 'checkm_contamination',
        'genome_size', 'gc_percentage', 'gtdb_taxonomy',
        'gtdb_representative', 'mimag_high_quality',
        'mimag_medium_quality', 'mimag_low_quality',
        'contig_count', 'n50_contigs'
    ]
    df = pd.read_csv(meta_file, sep='\t', usecols=cols, low_memory=False)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bac-metadata', required=True)
    parser.add_argument('--ar-metadata', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    bac = load_metadata(Path(args.bac_metadata))
    ar = load_metadata(Path(args.ar_metadata))
    df = pd.concat([bac, ar], ignore_index=True)

    # Clean accessions (strip GB_/RS_ prefixes if present)
    df['accession_clean'] = df['accession'].str.replace(r'^(GB_|RS_)', '', regex=True)

    # Parse taxonomy
    tax_df = df['gtdb_taxonomy'].apply(parse_taxonomy).apply(pd.Series)
    df = pd.concat([df, tax_df], axis=1)

    # Quality category
    df['quality_category'] = 'unknown'
    df.loc[df['mimag_high_quality'] == 't', 'quality_category'] = 'high'
    df.loc[df['mimag_medium_quality'] == 't', 'quality_category'] = 'medium'
    df.loc[df['mimag_low_quality'] == 't', 'quality_category'] = 'low'

    # Representative flag
    df['is_representative'] = df['gtdb_representative'] == 't'

    lines = []
    lines.append('# GTDB-R207 Summary Statistics Report')
    lines.append('')
    lines.append(f'Generated from: `{args.bac_metadata}` and `{args.ar_metadata}`')
    lines.append('')

    # Overall counts
    lines.append('## Genome Counts')
    lines.append('')
    lines.append(f'- **Total genomes**: {len(df):,}')
    lines.append(f'- **Bacteria**: {len(bac):,}')
    lines.append(f'- **Archaea**: {len(ar):,}')
    lines.append(f'- **GTDB representatives**: {df["is_representative"].sum():,}')
    lines.append('')

    # Domain distribution
    lines.append('## Domain Distribution')
    lines.append('')
    domain_counts = df['domain'].value_counts()
    for domain, count in domain_counts.items():
        lines.append(f'- {domain}: {count:,} ({count/len(df)*100:.2f}%)')
    lines.append('')

    # Phylum distribution (top 20)
    lines.append('## Top 20 Phyla')
    lines.append('')
    phylum_counts = df['phylum'].value_counts().head(20)
    for phylum, count in phylum_counts.items():
        lines.append(f'- {phylum}: {count:,} ({count/len(df)*100:.2f}%)')
    lines.append('')

    # Quality stats
    lines.append('## Quality Metrics')
    lines.append('')
    lines.append('### CheckM Completeness')
    desc = df['checkm_completeness'].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    lines.append(f'- count: {int(desc["count"]):,}')
    for stat in ['mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max']:
        lines.append(f'- {stat}: {desc[stat]:.2f}%')
    lines.append('')

    lines.append('### CheckM Contamination')
    desc = df['checkm_contamination'].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    lines.append(f'- count: {int(desc["count"]):,}')
    for stat in ['mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max']:
        lines.append(f'- {stat}: {desc[stat]:.2f}%')
    lines.append('')

    lines.append('### MIMAG Quality Category')
    qc = df['quality_category'].value_counts()
    for cat, count in qc.items():
        lines.append(f'- {cat}: {count:,} ({count/len(df)*100:.2f}%)')
    lines.append('')

    # Genome size
    lines.append('## Genome Size (Mb)')
    lines.append('')
    df['genome_size_mb'] = df['genome_size'] / 1e6
    desc = df['genome_size_mb'].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    lines.append(f'- count: {int(desc["count"]):,}')
    for stat in ['mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max']:
        lines.append(f'- {stat}: {desc[stat]:.3f} Mb')
    lines.append('')

    # GC content
    lines.append('## GC Content')
    lines.append('')
    desc = df['gc_percentage'].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    lines.append(f'- count: {int(desc["count"]):,}')
    for stat in ['mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max']:
        lines.append(f'- {stat}: {desc[stat]:.2f}%')
    lines.append('')

    # Contig count
    lines.append('## Contig Count')
    lines.append('')
    df['contig_count'] = pd.to_numeric(df['contig_count'], errors='coerce')
    desc = df['contig_count'].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    lines.append(f'- count: {int(desc["count"]):,}')
    for stat in ['mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max']:
        lines.append(f'- {stat}: {desc[stat]:.1f}')
    lines.append('')

    # N50
    lines.append('## Contig N50 (kb)')
    lines.append('')
    df['n50_contigs'] = pd.to_numeric(df['n50_contigs'], errors='coerce')
    df['n50_contigs_kb'] = df['n50_contigs'] / 1000
    desc = df['n50_contigs_kb'].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    lines.append(f'- count: {int(desc["count"]):,}')
    for stat in ['mean', 'std', 'min', '5%', '25%', '50%', '75%', '95%', 'max']:
        lines.append(f'- {stat}: {desc[stat]:.2f} kb')
    lines.append('')

    # Representatives by domain
    lines.append('## Representatives by Domain')
    lines.append('')
    rep_by_domain = df[df['is_representative']]['domain'].value_counts()
    for domain, count in rep_by_domain.items():
        lines.append(f'- {domain}: {count:,}')
    lines.append('')

    # Representatives by quality
    lines.append('## Representatives by Quality Category')
    lines.append('')
    rep_qc = df[df['is_representative']]['quality_category'].value_counts()
    for cat, count in rep_qc.items():
        lines.append(f'- {cat}: {count:,}')
    lines.append('')

    # Top species
    lines.append('## Top 20 Species (by genome count)')
    lines.append('')
    sp_counts = df['species'].value_counts().head(20)
    for sp, count in sp_counts.items():
        lines.append(f'- {sp}: {count:,}')
    lines.append('')

    output_path = Path(args.output)
    output_path.write_text('\n'.join(lines))
    print(f'Report saved: {output_path}')


if __name__ == '__main__':
    main()
