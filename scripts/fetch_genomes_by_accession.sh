#!/bin/bash
# fetch_genomes_by_accession.sh <accessions.txt> <outdir>
#
# Resolve GCF_/GCA_ accessions to NCBI FTP paths via assembly_summary.txt and
# download <acc>.fna into outdir. Idempotent: existing non-empty files are kept,
# so it can be re-run to fill gaps after a partial fetch.
#
# assembly_summary.txt is the authoritative accession -> path map, which is why
# this resolves rather than constructing URLs: the assembly *name* is part of the
# path and cannot be derived from the accession.
#
# If NCBI is unreachable from the compute nodes (it has been blocked here
# before), fetch the two summary files from a machine that can reach it, drop
# them in $SUMMARY_DIR, and re-run -- the download loop below only needs the
# ftp_path column. The established fallback for the genomes themselves is ENA,
# converting GCF_ to GCA_.
set -uo pipefail

ACCS=${1:?usage: fetch_genomes_by_accession.sh <accessions.txt> <outdir>}
OUT=${2:?usage: fetch_genomes_by_accession.sh <accessions.txt> <outdir>}
SUMMARY_DIR=${SUMMARY_DIR:-$OUT/.summary}
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes

mkdir -p "$OUT" "$SUMMARY_DIR"

for db in refseq genbank; do
    f=$SUMMARY_DIR/assembly_summary_$db.txt
    if [ ! -s "$f" ]; then
        echo "[fetch] downloading assembly_summary_$db.txt"
        wget -q -O "$f" "$NCBI/$db/bacteria/assembly_summary.txt" || {
            echo "[fetch] FAILED to get assembly_summary_$db.txt -- see header note"
            rm -f "$f"
        }
    fi
done

MAP=$SUMMARY_DIR/acc2path.tsv
if [ ! -s "$MAP" ]; then
    # column 1 = assembly accession, column 20 = ftp_path
    awk -F'\t' '!/^#/ && $20 != "na" {print $1"\t"$20}' \
        "$SUMMARY_DIR"/assembly_summary_*.txt > "$MAP"
fi
echo "[fetch] $(wc -l < "$MAP") accessions in the path map"

ok=0; skip=0; fail=0
while read -r acc; do
    [ -z "$acc" ] && continue
    dest=$OUT/${acc}.fna
    if [ -s "$dest" ]; then skip=$((skip+1)); continue; fi
    path=$(awk -F'\t' -v a="$acc" '$1==a {print $2; exit}' "$MAP")
    if [ -z "$path" ]; then
        echo "[fetch] no path for $acc"; fail=$((fail+1)); continue
    fi
    url="$path/$(basename "$path")_genomic.fna.gz"
    if wget -q -O "$dest.gz" "$url" && gunzip -f "$dest.gz"; then
        ok=$((ok+1))
    else
        echo "[fetch] failed $acc ($url)"; rm -f "$dest.gz" "$dest"; fail=$((fail+1))
    fi
done < "$ACCS"

echo "[fetch] downloaded $ok, already present $skip, failed $fail -> $OUT"
