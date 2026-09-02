#!/usr/bin/env python3
"""Download FDA-ARGOS S. aureus genome assemblies using NCBI datasets."""

import csv
import shutil
import subprocess
import sys
from pathlib import Path

WORK_DIR = Path("/Volumes/MoneyCat/Data/fda_argos_staphylococcus_aureus")
GENOMES_DIR = WORK_DIR / "genomes"
FAIL_LOG = WORK_DIR / "download_failures.log"


def download_one(accession):
    zip_path = WORK_DIR / f"{accession}.zip"
    out_fna = GENOMES_DIR / f"{accession}.fna"
    if out_fna.exists() and out_fna.stat().st_size > 0:
        return True, "already exists"

    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        accession,
        "--filename",
        str(zip_path),
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
    except subprocess.CalledProcessError as e:
        return False, f"datasets download failed: {e.stderr or e.stdout}"
    except subprocess.TimeoutExpired:
        return False, "datasets download timed out"

    extract_dir = WORK_DIR / f"extract_{accession}"
    try:
        shutil.unpack_archive(zip_path, extract_dir)
    except Exception as e:
        return False, f"unzip failed: {e}"

    fna_files = list(extract_dir.rglob("*_genomic.fna"))
    if not fna_files:
        return False, "no genomic FASTA found in archive"
    if len(fna_files) > 1:
        return False, f"multiple genomic FASTAs found: {fna_files}"

    shutil.move(str(fna_files[0]), str(out_fna))
    zip_path.unlink(missing_ok=True)
    shutil.rmtree(extract_dir, ignore_errors=True)
    return True, "ok"


def main():
    GENOMES_DIR.mkdir(parents=True, exist_ok=True)
    metadata = WORK_DIR / "assembly_metadata.tsv"
    accessions = []
    with open(metadata) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            accessions.append(row["assembly_acc"])

    failures = []
    success = 0
    for i, acc in enumerate(accessions, 1):
        print(f"[{i}/{len(accessions)}] Downloading {acc} ...", flush=True)
        ok, msg = download_one(acc)
        if ok:
            success += 1
            print(f"  -> {msg}")
        else:
            failures.append((acc, msg))
            print(f"  -> FAILED: {msg}", file=sys.stderr)

    with open(FAIL_LOG, "w") as fh:
        for acc, msg in failures:
            fh.write(f"{acc}\t{msg}\n")

    print(f"\nDownloaded {success}/{len(accessions)} genomes")
    print(f"Failures logged to {FAIL_LOG}")


if __name__ == "__main__":
    main()
