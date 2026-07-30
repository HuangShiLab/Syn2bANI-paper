#!/usr/bin/env python3
"""Build a download manifest for SynTracker isolate FASTQs from ENA."""
import argparse, csv, sys, time, urllib.request, urllib.error
from pathlib import Path

ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = "run_accession,fastq_ftp,fastq_bytes,submitted_ftp,library_layout"


def fetch_manifest(accessions):
    """Query ENA filereport for a batch of run accessions."""
    # ENA expects repeated 'accession=' parameters, not a comma-separated list.
    acc_params = '&'.join(f"accession={a}" for a in accessions)
    url = f"{ENA_URL}?{acc_params}&result=read_run&fields={FIELDS}"
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:
            text = resp.read().decode('utf-8')
    except urllib.error.URLError as e:
        raise RuntimeError(f"ENA query failed: {e}")
    rows = []
    for line in text.strip().splitlines()[1:]:
        cols = line.split('\t')
        if len(cols) < 4:
            continue
        run = cols[0]
        fastq_ftp = cols[1] if len(cols) > 1 else ''
        fastq_bytes = cols[2] if len(cols) > 2 else ''
        submitted_ftp = cols[3] if len(cols) > 3 else ''
        layout = cols[4] if len(cols) > 4 else ''
        rows.append({
            'run_accession': run,
            'fastq_ftp': fastq_ftp,
            'fastq_bytes': fastq_bytes,
            'submitted_ftp': submitted_ftp,
            'library_layout': layout,
        })
    return rows


def build_manifest(sample_tsv, batch_size=200):
    """Read sample TSVs and merge with ENA metadata."""
    samples = []
    for p in sample_tsv:
        with open(p, newline='') as fh:
            rdr = csv.DictReader(fh, delimiter='\t')
            for row in rdr:
                samples.append(row)

    run_to_sample = {}
    for s in samples:
        run = s['sra_run'].strip()
        run_to_sample[run] = s

    accessions = list(run_to_sample.keys())
    ena_rows = []
    # ENA filereport appears to honour only the first accession when multiple
    # are supplied, so query one run at a time to stay robust.
    for i, acc in enumerate(accessions):
        if i % 10 == 0:
            print(f"Querying ENA: {i}/{len(accessions)}", file=sys.stderr)
        try:
            ena_rows.extend(fetch_manifest([acc]))
        except RuntimeError as e:
            print(f"WARN: {e}", file=sys.stderr)
        time.sleep(0.2)

    out = []
    for e in ena_rows:
        run = e['run_accession']
        s = run_to_sample.get(run, {})
        urls = [u.strip() for u in e['fastq_ftp'].split(';') if u.strip()]
        if not urls and e['submitted_ftp']:
            urls = [u.strip() for u in e['submitted_ftp'].split(';') if u.strip()]
        out.append({
            'isolate': s.get('isolate', run),
            'species': s.get('species', ''),
            'sra_run': run,
            'library_layout': e['library_layout'],
            'url1': urls[0] if len(urls) > 0 else '',
            'url2': urls[1] if len(urls) > 1 else '',
            'bytes1': e['fastq_bytes'].split(';')[0] if e['fastq_bytes'] else '',
            'bytes2': e['fastq_bytes'].split(';')[1] if e['fastq_bytes'] and ';' in e['fastq_bytes'] else '',
        })
    return out


def main():
    parser = argparse.ArgumentParser(description="Build ENA FASTQ manifest")
    parser.add_argument('--samples', nargs='+', required=True, help='sample TSVs')
    parser.add_argument('--out', required=True, help='output manifest TSV')
    args = parser.parse_args()

    manifest = build_manifest(args.samples)
    fieldnames = ['isolate', 'species', 'sra_run', 'library_layout', 'url1', 'url2', 'bytes1', 'bytes2']
    with open(args.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        w.writerows(manifest)
    print(f"Wrote {len(manifest)} rows to {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
