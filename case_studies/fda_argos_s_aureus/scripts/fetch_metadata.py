#!/usr/bin/env python3
"""Fetch FDA-ARGOS Staphylococcus aureus assembly metadata from NCBI."""

import csv
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

OUTDIR = Path("/Volumes/MoneyCat/Data/fda_argos_staphylococcus_aureus")
OUTDIR.mkdir(parents=True, exist_ok=True)

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

TERM = "Staphylococcus aureus[Organism] AND PRJNA231221[bioproject]"


def get_uids():
    params = {"db": "assembly", "term": TERM, "retmax": 200}
    r = requests.get(ESEARCH_URL, params=params, timeout=120)
    r.raise_for_status()
    root = ET.fromstring(r.text)
    ids = [elem.text for elem in root.findall(".//Id")]
    print(f"esearch returned {len(ids)} UIDs")
    return ids


def get_assembly_summaries(uids, batch=50):
    summaries = []
    for i in range(0, len(uids), batch):
        batch_ids = uids[i : i + batch]
        params = {"db": "assembly", "id": ",".join(batch_ids)}
        r = requests.get(ESUMMARY_URL, params=params, timeout=120)
        r.raise_for_status()
        root = ET.fromstring(r.text)
        for ds in root.findall(".//DocumentSummary"):
            summaries.append(ds)
        time.sleep(0.4)
    return summaries


def extract_strain(biosource):
    for infraspecie in biosource.findall(".//Infraspecie"):
        sub_type = infraspecie.find("Sub_type")
        if sub_type is not None and sub_type.text == "strain":
            sub_value = infraspecie.find("Sub_value")
            if sub_value is not None:
                return sub_value.text
    return ""


def fetch_biosamples(biosample_ids, batch=50):
    records = {}
    for i in range(0, len(biosample_ids), batch):
        batch_ids = biosample_ids[i : i + batch]
        params = {"db": "biosample", "id": ",".join(batch_ids)}
        r = requests.get(EFETCH_URL, params=params, timeout=120)
        r.raise_for_status()
        root = ET.fromstring(r.text)
        for bs in root.findall("BioSample"):
            acc = bs.get("accession", "")
            attrs = {}
            for attr in bs.findall(".//Attribute"):
                name = attr.get("harmonized_name") or attr.get("attribute_name", "")
                attrs[name] = attr.text or ""
            records[acc] = attrs
        time.sleep(0.4)
    return records


def main():
    uids = get_uids()
    summaries = get_assembly_summaries(uids)

    rows = []
    biosample_ids = []
    for ds in summaries:
        acc = ds.findtext("AssemblyAccession", "")
        name = ds.findtext("AssemblyName", "")
        organism = ds.findtext("Organism", "")
        biosample = ds.findtext("BioSampleAccn", "")
        biosource = ds.find("Biosource")
        strain = extract_strain(biosource) if biosource is not None else ""
        rows.append(
            {
                "assembly_acc": acc,
                "assembly_name": name,
                "organism": organism,
                "strain": strain,
                "biosample": biosample,
                "country": "",
                "isolation_source": "",
                "collection_date": "",
            }
        )
        if biosample:
            biosample_ids.append(biosample)

    print(f"Collected {len(rows)} assemblies; fetching {len(biosample_ids)} BioSamples")
    biosample_attrs = fetch_biosamples(biosample_ids)

    for row in rows:
        attrs = biosample_attrs.get(row["biosample"], {})
        row["country"] = attrs.get("geo_loc_name", "")
        row["isolation_source"] = attrs.get("isolation_source", "")
        row["collection_date"] = attrs.get("collection_date", "")

    out_tsv = OUTDIR / "assembly_metadata.tsv"
    fieldnames = [
        "assembly_acc",
        "assembly_name",
        "organism",
        "strain",
        "biosample",
        "country",
        "isolation_source",
        "collection_date",
    ]
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {out_tsv}")


if __name__ == "__main__":
    main()
