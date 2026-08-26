#!/bin/bash
# s1_prepare_genomes.sh — download H. pylori assemblies from an accession list.
# Usage: s1_prepare_genomes.sh <accessions.txt> <outdir>
# accessions.txt: one GCA_/GCF_ accession per line (version optional).
set -euo pipefail

ACC=$1
OUT=$2
mkdir -p "$OUT"

while read -r acc; do
    [ -z "$acc" ] && continue
    fasta="$OUT/${acc}.fna"
    if [ -s "$fasta" ]; then
        continue
    fi
    # Resolve the FTP path via NCBI datasets v2 API (no API key needed at low rate).
    url=$(curl -s "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${acc}/download_summary" \
        | python3 -c "import json,sys; d=json.load(sys.stdin); print(d['reports'][0]['ftp_path'])" 2>/dev/null || true)
    if [ -z "${url:-}" ] || [ "$url" = "None" ]; then
        echo "WARN: cannot resolve $acc" >&2
        continue
    fi
    base=$(basename "$url")
    curl -s "${url}/${base}_genomic.fna.gz" -o "${fasta}.gz" \
        && gunzip -f "${fasta}.gz" \
        || { echo "WARN: download failed $acc" >&2; rm -f "${fasta}.gz"; }
done < "$ACC"

echo "downloaded: $(ls "$OUT"/*.fna 2>/dev/null | wc -l) / $(grep -c . "$ACC")"
