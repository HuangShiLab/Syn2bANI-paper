#!/usr/bin/env python3
"""Build a merged metadata table from accession_map.tsv and NCBI Datasets summaries."""
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

WORKDIR = Path(__file__).resolve().parent
OUT = WORKDIR / "metadata.tsv"
FAIL_LOG = WORKDIR / "metadata_failures.log"


def biosample_attrs(attrs):
    d = {}
    for attr in attrs:
        d[attr.get("name", "")] = attr.get("value", "")
    return d


def fetch_summary(acc):
    cmd = ["datasets", "summary", "genome", "accession", acc]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "datasets failed")
    data = json.loads(proc.stdout)
    reports = data.get("reports", [])
    if not reports:
        raise ValueError("empty reports")
    r = reports[0]
    assembly_info = r.get("assembly_info", {})
    biosample = assembly_info.get("biosample", {})
    attrs = biosample_attrs(biosample.get("attributes", []))
    return {
        "assembly_acc": acc,
        "assembly_name": assembly_info.get("assembly_name", ""),
        "organism": r.get("organism", {}).get("organism_name", ""),
        "status": assembly_info.get("assembly_status", ""),
        "strain": attrs.get("strain", biosample.get("strain", "")),
        "country": attrs.get("geo_loc_name", ""),
        "isolation_source": attrs.get("isolation_source", ""),
        "collection_date": attrs.get("collection_date", ""),
        "serotype": attrs.get("serotype", ""),
        "biosample": biosample.get("accession", ""),
    }


def main():
    with (WORKDIR / "accession_map.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    summaries = []
    failures = []
    for i, row in enumerate(rows, 1):
        acc = row["assembly_acc"]
        print(f"[{i}/{len(rows)}] metadata for {acc}", file=sys.stderr)
        try:
            summary = fetch_summary(acc)
        except Exception as e:
            failures.append((acc, str(e)))
            summary = {
                "assembly_acc": acc, "assembly_name": "", "organism": "",
                "status": "", "strain": "", "country": "", "isolation_source": "",
                "collection_date": "", "serotype": "", "biosample": ""
            }
        summary["nucleotide_acc"] = row["nucleotide_acc"]
        summaries.append(summary)
        time.sleep(0.2)

    fieldnames = [
        "nucleotide_acc", "assembly_acc", "assembly_name", "organism", "status",
        "strain", "country", "isolation_source", "collection_date", "serotype", "biosample"
    ]
    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    if failures:
        with FAIL_LOG.open("w") as fh:
            for acc, err in failures:
                fh.write(f"{acc}\t{err}\n")
    print(f"Metadata written to {OUT}. Failures: {len(failures)}", file=sys.stderr)


if __name__ == "__main__":
    main()
