#!/usr/bin/env python3
"""Download genomic FASTA files for a list of assembly accessions."""
import csv
import json
import subprocess
import sys
import time
import urllib.request
import urllib.error
import zipfile
import shutil
from pathlib import Path

WORKDIR = Path(__file__).resolve().parent
GENOMES_DIR = WORKDIR / "genomes"
FAIL_LOG = WORKDIR / "download_failures.log"


def accession_to_ftp_path(acc):
    """GCF_000008865.2 -> /genomes/all/GCF/000/008/865/GCF_000008865.2/"""
    prefix = acc[:3]  # GCF or GCA
    digits = acc.split("_")[1].split(".")[0]
    parts = [digits[i:i+3] for i in range(0, len(digits), 3)]
    return f"/genomes/all/{prefix}/{'/'.join(parts)}/{acc}/"


def ftp_fasta_url(acc):
    return f"https://ftp.ncbi.nlm.nih.gov{accession_to_ftp_path(acc)}{acc}_genomic.fna.gz"


def run_datasets(acc):
    zip_path = GENOMES_DIR / f"{acc}.zip"
    cmd = [
        "datasets", "download", "genome", "accession", acc,
        "--filename", str(zip_path), "--no-progressbar"
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    extract_dir = GENOMES_DIR / f"{acc}_unzip"
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)
    # find genomic fasta
    fna_files = list(extract_dir.rglob("*_genomic.fna"))
    if not fna_files:
        raise FileNotFoundError(f"No *_genomic.fna found in {extract_dir}")
    shutil.move(str(fna_files[0]), str(GENOMES_DIR / f"{acc}.fna"))
    shutil.rmtree(extract_dir)
    zip_path.unlink()


def ftp_fallback(acc):
    url = ftp_fasta_url(acc)
    out = GENOMES_DIR / f"{acc}.fna.gz"
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            with out.open("wb") as fh:
                fh.write(r.read())
        subprocess.run(["gunzip", "-f", str(out)], check=True)
    except Exception:
        if out.exists():
            out.unlink()
        raise


def main():
    GENOMES_DIR.mkdir(exist_ok=True)
    with (WORKDIR / "accession_map.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        accessions = [row["assembly_acc"] for row in reader if row.get("assembly_acc")]

    failures = []
    for i, acc in enumerate(accessions, 1):
        fna = GENOMES_DIR / f"{acc}.fna"
        if fna.exists() and fna.stat().st_size > 0:
            print(f"[{i}/{len(accessions)}] {acc}: already present", file=sys.stderr)
            continue
        print(f"[{i}/{len(accessions)}] {acc}: downloading...", file=sys.stderr)
        try:
            run_datasets(acc)
            print(f"[{i}/{len(accessions)}] {acc}: datasets OK", file=sys.stderr)
        except Exception as e:
            print(f"[{i}/{len(accessions)}] {acc}: datasets failed ({e}), trying FTP...", file=sys.stderr)
            try:
                ftp_fallback(acc)
                print(f"[{i}/{len(accessions)}] {acc}: FTP OK", file=sys.stderr)
            except Exception as e2:
                failures.append((acc, f"datasets: {e}; ftp: {e2}"))
                print(f"[{i}/{len(accessions)}] {acc}: FAILED", file=sys.stderr)
        time.sleep(0.2)

    with FAIL_LOG.open("w") as fh:
        for acc, err in failures:
            fh.write(f"{acc}\t{err}\n")
    print(f"Downloaded {len(accessions)-len(failures)}/{len(accessions)}. Failures: {len(failures)}", file=sys.stderr)


if __name__ == "__main__":
    main()
