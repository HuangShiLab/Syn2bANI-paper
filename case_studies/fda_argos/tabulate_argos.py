#!/usr/bin/env python3
"""Tabulate all FDA-ARGOS (PRJNA231221) assemblies per organism via NCBI eutils."""
import json
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
            wait = 2 ** attempt
            print(f"  retry {attempt+1} after error: {e} (sleep {wait}s)")
            time.sleep(wait)


def get_uids():
    data = http_get(f"{BASE}/esearch.fcgi", {
        "db": "assembly", "term": "PRJNA231221[bioproject]", "retmax": 5000})
    root = ET.fromstring(data)
    ids = [e.text for e in root.findall(".//Id")]
    print(f"esearch returned {len(ids)} UIDs")
    return ids


def get_summaries(uids, batch=200):
    docs = []
    for i in range(0, len(uids), batch):
        chunk = uids[i:i + batch]
        data = http_get(f"{BASE}/esummary.fcgi", {
            "db": "assembly", "id": ",".join(chunk), "retmode": "json"})
        js = json.loads(data)
        for uid in js["result"]["uids"]:
            docs.append(js["result"][uid])
        print(f"  summarized {len(docs)}/{len(uids)}")
        time.sleep(0.4)
    return docs


def main():
    uids = get_uids()
    docs = get_summaries(uids)

    rows = []
    for d in docs:
        biosource = d.get("biosource", {}) or {}
        strain = ""
        for infra in biosource.get("infraspecieslist", []) or []:
            if infra.get("sub_type") == "strain":
                strain = infra.get("sub_value", "")
        rows.append({
            "assembly_acc": d.get("assemblyaccession", ""),
            "assembly_name": d.get("assemblyname", ""),
            "organism": d.get("organism", ""),
            "strain": strain,
            "biosample": d.get("biosampleaccn", ""),
            "seq_rel_date": d.get("seqreleasedate", ""),
            "submission": d.get("submissiondate", ""),
            "assembly_level": d.get("assemblystatus", ""),
        })

    counts = Counter(r["organism"] for r in rows)
    print("\nTop organisms:")
    for org, n in counts.most_common(40):
        print(f"{n:5d}  {org}")

    with open(OUTDIR / "all_argos_assemblies.tsv", "w") as fh:
        fh.write("\t".join(["assembly_acc", "assembly_name", "organism", "strain",
                            "biosample", "seq_rel_date", "submission",
                            "assembly_level"]) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[k]).replace("\t", " ") for k in
                               ["assembly_acc", "assembly_name", "organism",
                                "strain", "biosample", "seq_rel_date",
                                "submission", "assembly_level"]) + "\n")
    print(f"\nWrote {OUTDIR / 'all_argos_assemblies.tsv'} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
