#!/usr/bin/env python3
"""
sample_gtdb_r207_pairs_v2.py
Stratified sampling of genome pairs from the downloaded GTDB-R207 collection.

Produces a TSV with columns:
  query, reference, label, q_domain, r_domain, q_phylum, r_phylum

Labels:
  high   = same species (ANI typically >95%)
  mid_high = same genus, different species
  mid    = same phylum, different genus
  low    = different phylum

Usage:
  python3 scripts/sample_gtdb_r207_pairs_v2.py \
    --genomes ~/data/gtdb-r207/genomes_all \
    --bac-metadata ~/data/gtdb-r207/bac120_metadata_r207.tsv \
    --ar-metadata ~/data/gtdb-r207/ar53_metadata_r207.tsv \
    --output results/pairs_gtdb_r207.tsv \
    --n-per-label 250 \
    --seed 42
"""

import argparse
import csv
import random
from collections import defaultdict
from pathlib import Path
from typing import Optional


def parse_taxonomy(tax_str: str) -> dict:
    """Parse GTDB taxonomy string into dict."""
    parts = [p.strip() for p in tax_str.split(';')]
    keys = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    return {k: parts[i] if i < len(parts) else 'unknown'
            for i, k in enumerate(keys)}


def load_records(genomes_dir: Path, bac_meta: Path, ar_meta: Path) -> list[dict]:
    """Load metadata or taxonomy files and filter to genomes that exist on disk."""
    records = []
    for meta_file in [bac_meta, ar_meta]:
        with open(meta_file) as f:
            first_line = f.readline()
            f.seek(0)
            # GTDB taxonomy files have no header: accession\ttaxonomy
            is_taxonomy_file = not first_line.startswith('accession')
            if is_taxonomy_file:
                reader = csv.DictReader(f, delimiter='\t',
                                        fieldnames=['accession', 'gtdb_taxonomy'])
            else:
                reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                acc = row.get('accession', '')
                if not acc:
                    continue
                acc = acc.replace('GB_', '').replace('RS_', '')
                fna_path = genomes_dir / f'{acc}.fna'
                if not fna_path.exists():
                    continue
                tax = parse_taxonomy(row.get('gtdb_taxonomy', ''))
                records.append({
                    'accession': acc,
                    'file': str(fna_path),
                    **tax,
                })
    return records


def sample_pairs(records: list[dict], n_per_label: int, seed: int) -> list[dict]:
    """Sample stratified pairs."""
    rng = random.Random(seed)

    by_species = defaultdict(list)
    by_genus = defaultdict(list)
    by_phylum = defaultdict(list)
    for r in records:
        by_species[r['species']].append(r)
        by_genus[r['genus']].append(r)
        by_phylum[r['phylum']].append(r)

    pairs = []
    pair_keys = set()

    def add_pair(q: dict, r: dict, label: str) -> bool:
        key = tuple(sorted([q['accession'], r['accession']]))
        if key not in pair_keys and q['accession'] != r['accession']:
            pair_keys.add(key)
            pairs.append({
                'query': q['accession'],
                'reference': r['accession'],
                'label': label,
                'q_domain': q['domain'],
                'r_domain': r['domain'],
                'q_phylum': q['phylum'],
                'r_phylum': r['phylum'],
                'q_genus': q['genus'],
                'r_genus': r['genus'],
                'q_species': q['species'],
                'r_species': r['species'],
            })
            return True
        return False

    # high: same species
    multi_species = {sp: g for sp, g in by_species.items() if len(g) >= 2}
    count = 0
    items = list(multi_species.items())
    rng.shuffle(items)
    for sp, genomes in items:
        rng.shuffle(genomes)
        for i in range(min(len(genomes), 4)):
            for j in range(i + 1, min(len(genomes), 5)):
                if add_pair(genomes[i], genomes[j], 'high'):
                    count += 1
                    if count >= n_per_label:
                        break
            if count >= n_per_label:
                break
        if count >= n_per_label:
            break
    print(f'high pairs: {count}')

    # mid_high: same genus, different species
    multi_genus = {g: list({r['species'] for r in genomes})
                   for g, genomes in by_genus.items()
                   if len({r['species'] for r in genomes}) >= 2}
    count = 0
    items = list(multi_genus.items())
    rng.shuffle(items)
    for g, species_list in items:
        rng.shuffle(species_list)
        for _ in range(min(5, len(species_list) * (len(species_list) - 1) // 2)):
            s1, s2 = rng.sample(species_list, 2)
            g1 = rng.choice([r for r in by_genus[g] if r['species'] == s1])
            g2 = rng.choice([r for r in by_genus[g] if r['species'] == s2])
            if add_pair(g1, g2, 'mid_high'):
                count += 1
                if count >= n_per_label:
                    break
        if count >= n_per_label:
            break
    print(f'mid_high pairs: {count}')

    # mid: same phylum, different genus
    count = 0
    items = list(by_phylum.items())
    rng.shuffle(items)
    for phylum, genomes in items:
        genera = list({r['genus'] for r in genomes})
        if len(genera) >= 2:
            rng.shuffle(genera)
            for _ in range(min(5, len(genera) * (len(genera) - 1) // 2)):
                g1, g2 = rng.sample(genera, 2)
                r1 = rng.choice([r for r in genomes if r['genus'] == g1])
                r2 = rng.choice([r for r in genomes if r['genus'] == g2])
                if add_pair(r1, r2, 'mid'):
                    count += 1
                    if count >= n_per_label:
                        break
        if count >= n_per_label:
            break
    print(f'mid pairs: {count}')

    # low: different phylum
    phyla = list(by_phylum.keys())
    count = 0
    while count < n_per_label:
        p1, p2 = rng.sample(phyla, 2)
        r1 = rng.choice(by_phylum[p1])
        r2 = rng.choice(by_phylum[p2])
        if add_pair(r1, r2, 'low'):
            count += 1
    print(f'low pairs: {count}')

    return pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genomes', required=True)
    parser.add_argument('--bac-metadata', required=True)
    parser.add_argument('--ar-metadata', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--n-per-label', type=int, default=250)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    genomes_dir = Path(args.genomes)
    records = load_records(genomes_dir, Path(args.bac_metadata), Path(args.ar_metadata))
    print(f'Total genomes available: {len(records)}')

    pairs = sample_pairs(records, args.n_per_label, args.seed)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'query', 'reference', 'label',
            'q_domain', 'r_domain', 'q_phylum', 'r_phylum',
            'q_genus', 'r_genus', 'q_species', 'r_species'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(pairs)

    print(f'\nSaved {len(pairs)} pairs to {out_path}')
    print('Label distribution:')
    from collections import Counter
    for label, count in Counter(p['label'] for p in pairs).items():
        print(f'  {label}: {count}')


if __name__ == '__main__':
    main()
