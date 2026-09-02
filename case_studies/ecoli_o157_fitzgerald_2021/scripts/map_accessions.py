#!/usr/bin/env python3
"""Map nucleotide accessions (chromosome/plasmid) to NCBI assembly accessions."""
import time
import json
import csv
import sys
from pathlib import Path
import urllib.request
import urllib.error

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ACCESSIONS = [
    "CP018237", "CP015832", "CP018239", "CP018241", "CP018243", "CP018245",
    "CP015831", "CP018247", "CP018252", "CP018250", "NC_011353", "CP008957",
    "NC_002695", "CP010304", "NC_013008", "CP038414", "CP038402", "CP038372",
    "CP038366", "CP038360", "CP038357", "CP038353", "CP038349", "CP038344",
    "CP038416", "CP038342", "CP038309", "CP039834", "CP038333", "CP038292",
    "CP038290", "CP038302", "CP038300", "CP038346", "CP038355", "CP038351",
    "CP038282", "CP038339", "CP038319", "CP038284", "CP062749", "CP062746",
    "CP062744", "CP062742", "CP062736", "CP062733", "CP062731", "CP062780",
    "CP062729", "CP062727", "CP062725", "CP062723", "CP062721", "CP062719",
    "CP062717", "CP062715", "CP062713", "CP062711", "CP062708", "CP062705",
    "CP062702", "CP062700", "CP062778", "CP062774", "CP062771", "CP062769",
    "CP062766", "CP062761", "CP062758", "CP062755", "CP062752", "CP062763",
    "CP062782", "CP062739"
]


def get_json(url, retries=5):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (429, 502, 503, 504):
                time.sleep(2 ** attempt)
                continue
            raise
        except Exception:
            time.sleep(2 ** attempt)
            continue
    raise RuntimeError(f"Failed to fetch {url}")


def nucleotide_to_assembly_uid(nuc_acc):
    # Search the nucleotide database first (handles aliases)
    search_url = f"{BASE}/esearch.fcgi?db=nuccore&term={nuc_acc}&retmode=json"
    data = get_json(search_url)
    ids = data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None, f"no nuccore id for {nuc_acc}"
    uid = ids[0]
    # Link to assembly
    link_url = f"{BASE}/elink.fcgi?dbfrom=nuccore&db=assembly&id={uid}&retmode=json"
    link = get_json(link_url)
    uids = []
    for item in link.get("linksets", []):
        for linkset in item.get("linksetdbs", []):
            if linkset.get("dbto") == "assembly":
                uids.extend(linkset.get("links", []))
    if not uids:
        return None, f"no assembly link for {nuc_acc} (nuccore {uid})"
    return uids[0], None


def assembly_summary(uid):
    url = f"{BASE}/esummary.fcgi?db=assembly&id={uid}&retmode=json"
    data = get_json(url)
    result = data.get("result", {})
    docsum = result.get(str(uid), {})
    return {
        "assembly_acc": docsum.get("assemblyaccession"),
        "assembly_name": docsum.get("assemblyname"),
        "strain": docsum.get("organism", "").split("str.")[-1].strip() if "str." in docsum.get("organism", "") else "",
        "organism": docsum.get("organism"),
        "status": docsum.get("assemblystatus"),
    }


def main():
    out = Path("accession_map.tsv")
    rows = []
    failures = []
    for i, nuc_acc in enumerate(ACCESSIONS, 1):
        print(f"[{i}/{len(ACCESSIONS)}] {nuc_acc}", file=sys.stderr)
        try:
            uid, err = nucleotide_to_assembly_uid(nuc_acc)
            if err:
                failures.append((nuc_acc, err))
                continue
            summary = assembly_summary(uid)
            rows.append({
                "nucleotide_acc": nuc_acc,
                **summary,
            })
        except Exception as e:
            failures.append((nuc_acc, str(e)))
        time.sleep(0.35)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["nucleotide_acc", "assembly_acc", "assembly_name", "strain", "organism", "status"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if failures:
        with Path("map_failures.log").open("w") as fh:
            for acc, err in failures:
                fh.write(f"{acc}\t{err}\n")
        print(f"Mapping complete. Success={len(rows)}, Failures={len(failures)}", file=sys.stderr)
    else:
        print(f"Mapping complete. Success={len(rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()
