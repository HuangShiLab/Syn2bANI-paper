#!/usr/bin/env python3
"""Stream representative OAPGC FASTAs out of gzipped Zenodo tarballs.

The public archive contains all 149,921 HQMAG FASTAs, while the census uses
only the 99,215 non-null strain representatives (minus singleton SGBs).
This script makes one streaming pass over each tar.gz, extracting only
matching members, so it avoids expanding the redundant MAG set.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import tarfile
from pathlib import Path


def genome_stem(member_name: str) -> str | None:
    name = Path(member_name).name
    for suffix in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return None


def md5(path: Path, expected: str | None = None) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    digest = h.hexdigest()
    if expected and digest != expected:
        raise ValueError(f"MD5 mismatch for {path}: got {digest}, expected {expected}")
    return digest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archives", nargs="+", required=True,
                    help="OAPGC_genomes_part*.tar.gz, in any order")
    ap.add_argument("--expected-md5", nargs="*", default=[],
                    help="expected MD5 for each archive, same order")
    ap.add_argument("--representatives", required=True,
                    help="two-column genome_id<TAB>cluster list")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if len(args.expected_md5) not in (0, len(args.archives)):
        raise SystemExit("expected-md5 must be empty or match archives count")

    wanted = {}
    with open(args.representatives) as fh:
        for line in fh:
            genome, cluster = line.rstrip("\n").split("\t")[:2]
            wanted[genome] = cluster
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    done = out / ".extraction_done.tsv"
    if done.exists():
        extracted = {line.split("\t")[0] for line in done.read_text().splitlines()}
    else:
        extracted = set()

    todo = set(wanted) - extracted
    print(f"representatives={len(wanted)} already_extracted={len(extracted)} "
          f"todo={len(todo)}", flush=True)
    expected = dict(zip(args.archives, args.expected_md5))
    with done.open("a") as donefh:
        for archive in args.archives:
            archive_path = Path(archive)
            if not archive_path.is_file():
                raise FileNotFoundError(archive_path)
            if archive in expected:
                print(f"MD5 {archive}", flush=True)
                md5(archive_path, expected[archive])
            found = []
            with tarfile.open(archive_path, "r:gz") as tf:
                for member in tf:
                    if not member.isfile():
                        continue
                    genome = genome_stem(member.name)
                    if genome is None:
                        continue
                    if genome in todo:
                        source = tf.extractfile(member)
                        if source is None:
                            raise RuntimeError(f"cannot read tar member {member.name}")
                        tmp = out / f".{genome}.fna.tmp"
                        try:
                            with gzip.open(source, "rb") as fin, tmp.open("wb") as fout:
                                shutil.copyfileobj(fin, fout, length=8 * 1024 * 1024)
                            tmp.replace(out / f"{genome}.fna")
                        finally:
                            if tmp.exists():
                                tmp.unlink()
                        found.append(genome)
                        if len(found) % 1000 == 0:
                            print(f"{archive_path.name}: {len(found):,} files", flush=True)
            for genome in found:
                donefh.write(f"{genome}\t{wanted[genome]}\n")
            donefh.flush()
            todo.difference_update(found)
            print(f"DONE {archive_path.name}: extracted={len(found):,} remaining={len(todo):,}",
                  flush=True)
    if todo:
        missing = sorted(todo)[:20]
        raise SystemExit(f"missing {len(todo)} representatives; examples: {missing}")
    print(f"COMPLETE: {len(wanted):,} representative FASTAs in {out}", flush=True)


if __name__ == "__main__":
    main()
