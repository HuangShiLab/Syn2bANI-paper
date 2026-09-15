#!/usr/bin/env python3
"""Fetch genomic FASTA files for NCBI assembly accessions without the
`datasets` CLI: resolve the assembly name through the Datasets v2 REST API,
then pull `<acc>_<name>_genomic.fna.gz` from the NCBI FTP mirror.

Usage:
    python3 fetch_ncbi_fna.py --accessions list.txt --outdir genomes/

`list.txt` has one accession per line (GCF_/GCA_); a TSV with an
`assembly_acc` column is also accepted.
"""
import argparse
import csv
import gzip
import json
import os
import shutil
import sys
import time
import urllib.request

API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/dataset_report"
API_ALL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/dataset_report"
FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/{pre}/{a}/{b}/{c}/{acc}_{name}/{acc}_{name}_genomic.fna.gz"


def read_accessions(path):
    with open(path) as fh:
        first = fh.readline()
        fh.seek(0)
        if "\t" in first:
            return [r["assembly_acc"] for r in csv.DictReader(fh, delimiter="\t") if r.get("assembly_acc")]
        return [l.strip() for l in fh if l.strip()]


def resolve(acc):
    """Return (accession to download, assembly name).

    The per-accession endpoint returns an empty report for a superseded
    assembly version (56 of the 122 FDA-ARGOS *S. aureus* accessions are
    superseded). The batch endpoint with `all_assemblies` still reports those
    and names their replacement in `current_accession`, which is what the FTP
    mirror carries.
    """
    with urllib.request.urlopen(API.format(acc=acc), timeout=120) as r:
        reports = (json.load(r).get("reports") or [])
    if not reports:
        body = json.dumps({"accessions": [acc],
                           "filters": {"assembly_version": "all_assemblies"}}).encode()
        req = urllib.request.Request(API_ALL, data=body,
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=120) as r:
            reports = (json.load(r).get("reports") or [])
        if not reports:
            raise RuntimeError("no dataset report")
        cur = reports[0].get("current_accession") or acc
        if cur != acc:
            with urllib.request.urlopen(API.format(acc=cur), timeout=120) as r:
                cur_reports = (json.load(r).get("reports") or [])
            if cur_reports:
                return cur, cur_reports[0]["assembly_info"]["assembly_name"].replace(" ", "_")
    return acc, reports[0]["assembly_info"]["assembly_name"].replace(" ", "_")


def fetch(acc, outdir):
    out = os.path.join(outdir, f"{acc}.fna")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return "present"
    acc_dl, name = resolve(acc)
    digits = acc_dl.split("_")[1].split(".")[0]
    url = FTP.format(pre=acc_dl[:3], a=digits[0:3], b=digits[3:6], c=digits[6:9],
                     acc=acc_dl, name=name)
    tmp = out + ".gz"
    with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh)
    with gzip.open(tmp, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(tmp)
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    accs = read_accessions(args.accessions)
    failed = []
    for i, acc in enumerate(accs, 1):
        try:
            status = fetch(acc, args.outdir)
            print(f"[{i}/{len(accs)}] {acc}: {status}", file=sys.stderr)
        except Exception as e:
            failed.append(acc)
            print(f"[{i}/{len(accs)}] {acc}: FAILED ({e})", file=sys.stderr)
        time.sleep(0.34)  # NCBI API rate limit without a key: 3 requests/s
    print(f"{len(accs) - len(failed)}/{len(accs)} fetched; failed: {failed}", file=sys.stderr)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
