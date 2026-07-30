#!/usr/bin/env python3
"""Parse SynTracker Nature Biotech 2024 supplementary tables into clean metadata."""
import sys, argparse, re
from pathlib import Path

# Ensure local openpyxl is available
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / '.pydeps'))
import pandas as pd


def clean_df(df, expected_header=None):
    """Drop leading blank rows/cols and make first non-empty row the header.

    If expected_header is a set of column names, locate the first row that
    contains them and use it as the header.
    """
    # Drop rows that are completely blank
    df = df.dropna(how='all').reset_index(drop=True)
    # Drop columns that are completely blank
    df = df.dropna(axis=1, how='all').reset_index(drop=True)
    if df.empty:
        return df

    if expected_header:
        header_row = None
        for i, row in df.iterrows():
            cells = {str(c).strip() for c in row if pd.notna(c)}
            if expected_header <= cells:
                header_row = i
                break
        if header_row is None:
            header_row = 0
        df.columns = df.iloc[header_row]
        df = df.iloc[header_row + 1:].reset_index(drop=True)
    else:
        header_row = 0
        df.columns = df.iloc[header_row]
        df = df.iloc[1:].reset_index(drop=True)

    # Strip whitespace from column names
    df.columns = [str(c).strip() if c is not None else '' for c in df.columns]
    return df


def parse_single_column_sra(xl, sheet_name, species, prefix):
    """Parse sheets that are a single column with a title row followed by SRA runs."""
    df = xl.parse(sheet_name, header=None)
    df = df.dropna(how='all').reset_index(drop=True)
    # Drop leading title row(s) until we hit an SRA-style accession
    rows = []
    for _, row in df.iterrows():
        val = str(row.iloc[0]).strip() if pd.notna(row.iloc[0]) else ''
        if re.match(r'^[SED]RR\d+$', val):
            rows.append(val)
    df = pd.DataFrame({'sra_run': rows})
    df['species'] = species
    df['isolate'] = [f"{prefix}_{str(i+1).zfill(2)}" for i in range(len(df))]
    return df[['species', 'isolate', 'sra_run']]


def parse_table_s2(xl):
    return parse_single_column_sra(xl, 'Table_s2', 'Neisseria_gonorrhoeae', 'NG')


def parse_table_s3(xl):
    df = clean_df(xl.parse('table_s3'), expected_header={'Isolate', 'Mouse', 'BioSample', 'SRA', 'day'})
    # Expected columns: Isolate, Mouse, BioSample, SRA, day
    df = df.rename(columns={
        'Isolate': 'isolate',
        'Mouse': 'mouse',
        'BioSample': 'biosample',
        'SRA': 'sra_run',
        'day': 'day',
    })
    for col in ['isolate', 'mouse', 'biosample', 'sra_run', 'day']:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    df['species'] = 'Escherichia_coli_hypermutator'
    return df[['species', 'isolate', 'mouse', 'day', 'biosample', 'sra_run']]


def parse_table_s4(xl):
    df = clean_df(xl.parse('table_s4'), expected_header={'Isolate', 'SRA', 'host2'})
    # Expected columns: Isolate, SRA, host2
    df = df.rename(columns={
        'Isolate': 'isolate',
        'SRA': 'sra_run',
        'host2': 'host',
    })
    for col in ['isolate', 'sra_run', 'host']:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    df['species'] = 'Helicobacter_pylori'
    return df[['species', 'isolate', 'host', 'sra_run']]


def parse_table_s5(xl):
    return parse_single_column_sra(xl, 'table_s5', 'Streptomyces_rimosus', 'SR')


def references():
    return pd.DataFrame([
        {'species': 'Neisseria_gonorrhoeae', 'ref_accession': 'GCF_900087635.2'},
        {'species': 'Escherichia_coli_hypermutator', 'ref_accession': 'NC_000913.3'},
        {'species': 'Helicobacter_pylori', 'ref_accession': 'CP032479.1'},
        {'species': 'Streptomyces_rimosus', 'ref_accession': 'GCF_000331185.2'},
    ])


def make_pairs(samples, group_col=None):
    """Make within-species pair list. If group_col is given, only pairs within the same group."""
    pairs = []
    if group_col and group_col in samples.columns:
        for _, g in samples.groupby(group_col):
            ids = g['isolate'].tolist()
            for i in range(len(ids)):
                for j in range(i+1, len(ids)):
                    pairs.append((ids[i], ids[j]))
    else:
        ids = samples['isolate'].tolist()
        for i in range(len(ids)):
            for j in range(i+1, len(ids)):
                pairs.append((ids[i], ids[j]))
    return pd.DataFrame(pairs, columns=['query', 'reference'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--xlsx', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    xl = pd.ExcelFile(args.xlsx)

    ng = parse_table_s2(xl)
    ec = parse_table_s3(xl)
    hp = parse_table_s4(xl)
    sr = parse_table_s5(xl)

    ng.to_csv(outdir / 'samples_Neisseria_gonorrhoeae.tsv', sep='\t', index=False)
    ec.to_csv(outdir / 'samples_Escherichia_coli_hypermutator.tsv', sep='\t', index=False)
    hp.to_csv(outdir / 'samples_Helicobacter_pylori.tsv', sep='\t', index=False)
    sr.to_csv(outdir / 'samples_Streptomyces_rimosus.tsv', sep='\t', index=False)
    references().to_csv(outdir / 'references.tsv', sep='\t', index=False)

    make_pairs(ng).to_csv(outdir / 'pairs_Neisseria_gonorrhoeae.tsv', sep='\t', index=False)
    make_pairs(ec, group_col='mouse').to_csv(outdir / 'pairs_Escherichia_coli_hypermutator.tsv', sep='\t', index=False)
    make_pairs(hp, group_col='host').to_csv(outdir / 'pairs_Helicobacter_pylori.tsv', sep='\t', index=False)
    make_pairs(sr).to_csv(outdir / 'pairs_Streptomyces_rimosus.tsv', sep='\t', index=False)

    print(f"Wrote metadata to {outdir}")
    for sp, df in [('Neisseria_gonorrhoeae', ng), ('Escherichia_coli_hypermutator', ec),
                   ('Helicobacter_pylori', hp), ('Streptomyces_rimosus', sr)]:
        n_pairs = len(pd.read_csv(outdir / f'pairs_{sp}.tsv', sep='\t'))
        print(f"  {sp}: {len(df)} isolates -> {n_pairs} pairs")


if __name__ == '__main__':
    main()
