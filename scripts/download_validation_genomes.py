#!/usr/bin/env python3
"""Download oral/gut validation genomes from NCBI RefSeq using datasets CLI."""
import argparse
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple


# species_name: (category, max_to_download)
SPECIES = [
    # Oral
    ("Streptococcus mutans", "oral", 5),
    ("Porphyromonas gingivalis", "oral", 5),
    ("Fusobacterium nucleatum", "oral", 5),
    ("Aggregatibacter actinomycetemcomitans", "oral", 5),
    ("Streptococcus sanguinis", "oral", 5),
    # Gut
    ("Bacteroides fragilis", "gut", 5),
    ("Faecalibacterium prausnitzii", "gut", 5),
    ("Akkermansia muciniphila", "gut", 5),
    ("Bifidobacterium longum", "gut", 5),
    ("Roseburia intestinalis", "gut", 5),
]


def run_datasets_summary(name: str, datasets_path: str) -> List[dict]:
    cmd = [
        datasets_path, "summary", "genome", "taxon", name,
        "--tax-exact-match",
        "--assembly-source", "refseq",
        "--mag", "exclude",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error summarizing taxid {taxid}: {result.stderr}", file=sys.stderr)
        return []
    data = json.loads(result.stdout)
    return data.get("reports", [])


def pick_accessions(reports: List[dict], n: int, seed: int) -> List[str]:
    """Pick up to n accessions, preferring complete/chromosome assemblies."""
    # Sort by assembly level quality: complete > chromosome > scaffold > contig
    level_order = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    entries = []
    for r in reports:
        acc = r.get("accession")
        if not acc:
            continue
        level = r.get("assembly_info", {}).get("assembly_level", "Contig")
        entries.append((level_order.get(level, 9), acc))
    entries.sort()
    # Use deterministic sampling within each quality tier
    rng = random.Random(seed)
    selected = []
    for level, group in itertools.groupby(entries, key=lambda x: x[0]):
        group = list(group)
        rng.shuffle(group)
        selected.extend([acc for _, acc in group])
        if len(selected) >= n:
            break
    return selected[:n]


def download_genome(accession: str, out_dir: Path, datasets_path: str, retries: int = 3) -> bool:
    dst = out_dir / f"{accession}.fna"
    if dst.exists():
        return True
    tmp_zip = out_dir / f"{accession}.zip"
    for attempt in range(retries):
        cmd = [
            datasets_path, "download", "genome", "accession", accession,
            "--include", "genome",
            "--filename", str(tmp_zip),
            "--no-progressbar",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error downloading {accession} attempt {attempt+1}: {result.stderr}", file=sys.stderr)
            continue
        # Unzip
        unzip_dir = out_dir / f"{accession}_tmp"
        subprocess.run(["unzip", "-o", str(tmp_zip), "-d", str(unzip_dir)], capture_output=True)
        # Find .fna file
        fna_files = list(unzip_dir.rglob("*_genomic.fna"))
        if not fna_files:
            print(f"No .fna found for {accession}", file=sys.stderr)
            shutil.rmtree(unzip_dir, ignore_errors=True)
            tmp_zip.unlink(missing_ok=True)
            continue
        src = fna_files[0]
        shutil.move(str(src), str(dst))
        # Cleanup
        shutil.rmtree(unzip_dir, ignore_errors=True)
        tmp_zip.unlink(missing_ok=True)
        return True
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--datasets', default='/group/aos_shihuang/conda/bin/datasets')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genomes_dir = out_dir / 'genomes'
    genomes_dir.mkdir(exist_ok=True)

    manifest = []
    for name, category, n in SPECIES:
        print(f"\n=== {name} ===")
        reports = run_datasets_summary(name, args.datasets)
        print(f"Found {len(reports)} RefSeq genomes")
        if not reports:
            continue
        accessions = pick_accessions(reports, n, args.seed + hash(name) % 10000)
        print(f"Selected {len(accessions)}: {accessions}")
        for acc in accessions:
            ok = download_genome(acc, genomes_dir, args.datasets)
            if ok:
                manifest.append({
                    'accession': acc,
                    'species': name,
                    'category': category,
                    'file': str(genomes_dir / f"{acc}.fna"),
                })
                print(f"Downloaded {acc}")
            else:
                print(f"Failed {acc}")

    manifest_df = pd.DataFrame(manifest)
    manifest_path = out_dir / 'manifest.tsv'
    manifest_df.to_csv(manifest_path, sep='\t', index=False)
    print(f"\nManifest saved: {manifest_path}")
    print(f"Total downloaded: {len(manifest_df)}")


if __name__ == '__main__':
    import itertools
    import pandas as pd
    main()
