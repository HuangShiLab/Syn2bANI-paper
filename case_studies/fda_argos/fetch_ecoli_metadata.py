#!/usr/bin/env python3
"""Fetch BioSample attributes for the selected FDA-ARGOS E. coli accessions."""
import csv
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

OUTDIR = Path(__file__).resolve().parent
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def http_get(url, params, retries=4):
    qs = urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(f"{url}?{qs}", timeout=120) as r:
                return r.read()
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def main():
    accs = [l.strip() for l in open(OUTDIR / "acc_ecoli.txt") if l.strip()]
    meta = {}
    for r in csv.DictReader(open(OUTDIR / "all_argos_assemblies.tsv"), delimiter="\t"):
        meta[r["assembly_acc"]] = r

    biosamples = [meta[a]["biosample"] for a in accs if meta.get(a, {}).get("biosample")]
    print(f"{len(accs)} accessions, {len(biosamples)} biosamples")

    records = {}
    for i in range(0, len(biosamples), 50):
        chunk = biosamples[i:i + 50]
        data = http_get(f"{BASE}/efetch.fcgi", {"db": "biosample", "id": ",".join(chunk)})
        root = ET.fromstring(data)
        for bs in root.findall("BioSample"):
            acc = bs.get("accession", "")
            attrs = {}
            for attr in bs.findall(".//Attribute"):
                name = attr.get("harmonized_name") or attr.get("attribute_name", "")
                attrs[name] = (attr.text or "").strip()
            records[acc] = attrs
        print(f"  fetched {len(records)}/{len(biosamples)}")
        time.sleep(0.4)

    # attribute coverage
    cov = Counter()
    for attrs in records.values():
        for k in attrs:
            cov[k] += 1
    print("\nAttribute coverage (top 40):")
    for k, c in cov.most_common(40):
        print(f"  {c:4d}  {k}")

    keep_cols = ["strain", "serotype", "serovar", "isolation_source", "host",
                 "geo_loc_name", "collection_date", "host_disease",
                 "antibiogram", "antimicrobial_resistance", "resistance_phenotype",
                 "pathotype", "isolate", "collected_by", "sample_type"]
    rows = []
    for a in accs:
        m = meta.get(a, {})
        attrs = records.get(m.get("biosample", ""), {})
        row = {
            "assembly_acc": a,
            "organism": m.get("organism", "").split(" (")[0],
            "strain": m.get("strain", "") or attrs.get("strain", ""),
            "biosample": m.get("biosample", ""),
            "assembly_level": m.get("assembly_level", ""),
        }
        for c in keep_cols:
            row[c] = attrs.get(c, "")
        rows.append(row)

    fieldnames = ["assembly_acc", "organism", "strain", "biosample",
                  "assembly_level"] + keep_cols
    with open(OUTDIR / "accessions_used.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {OUTDIR / 'accessions_used.tsv'} ({len(rows)} rows)")

    # dump full attribute set for reference
    all_cols = sorted(cov)
    with open(OUTDIR / "biosample_attributes_full.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["biosample"] + all_cols)
        for bs, attrs in records.items():
            w.writerow([bs] + [attrs.get(c, "") for c in all_cols])
    print(f"Wrote {OUTDIR / 'biosample_attributes_full.tsv'}")


if __name__ == "__main__":
    main()
