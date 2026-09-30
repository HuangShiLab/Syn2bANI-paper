#!/usr/bin/env bash
# Download-waiting, extraction, and local OAPGC structural census driver.
# Safe to rerun: extraction and census tasks are checkpointed.
set -euo pipefail

RAW="${RAW:-/Users/macstudio/Downloads/OAPGC/raw}"
PLAN="${PLAN:-/Users/macstudio/Downloads/OAPGC/census_plan}"
GENOMES="${GENOMES:-/Users/macstudio/Downloads/OAPGC/genomes}"
WORK="${WORK:-/Users/macstudio/Downloads/OAPGC/census}"
RESULTS="${RESULTS:-/Users/macstudio/Downloads/OAPGC/census_results}"
REPO="${REPO:-/Users/macstudio/Downloads/Syn2bANI-paper}"
SYN2B="${SYN2B:-/Users/macstudio/Downloads/Syn2b/target/release/syn2b}"
CORES="${CORES:-12}"
DOWNLOAD_PID="${DOWNLOAD_PID:-$RAW/aria2.pid}"
ANNOTATION="${ANNOTATION:-$RAW/cluster.annotation.xlsx}"

cd "$REPO"
mkdir -p "$RAW" "$RESULTS"

if [[ -f "$DOWNLOAD_PID" ]]; then
    pid="$(cat "$DOWNLOAD_PID")"
    echo "waiting for downloader pid=$pid"
    while kill -0 "$pid" 2>/dev/null; do sleep 300; done
fi

python3 - <<'PY'
import json
from pathlib import Path
raw = Path('/Users/macstudio/Downloads/OAPGC/raw')
wanted = [('19212794', 'OAPGC_genomes_part1.tar.gz'),
          ('19212794', 'OAPGC_genomes_part2.tar.gz'),
          ('19220405', 'OAPGC_genomes_part3.tar.gz'),
          ('19220405', 'OAPGC_genomes_part4.tar.gz'),
          ('19220405', 'OAPGC_genomes_part5.tar.gz')]
with (raw / 'expected_md5.tsv').open('w') as out:
    for rid, key in wanted:
        rec = json.load(open(raw / f'zenodo_{rid}.json'))
        f = next(x for x in rec['files'] if x['key'] == key)
        out.write(f"{raw / key}\t{f['checksum'].split(':', 1)[1]}\n")
PY

archives=()
md5s=()
while IFS="$(printf '\t')" read -r archive checksum; do
    archives+=("$archive")
    md5s+=("$checksum")
done < "$RAW/expected_md5.tsv"

echo "extracting representative FASTAs"
python3 scripts/oapgc_census/extract_oapgc_representatives.py \
  --archives "${archives[@]}" \
  --representatives "$PLAN/representative_genomes.txt" \
  --out "$GENOMES"

echo "preparing census manifest and tasks"
python3 scripts/oapgc_census/prepare_oapgc_census.py \
  --annotation "$ANNOTATION" \
  --genome-root "$GENOMES" \
  --workdir "$WORK" \
  --batch-size 250 --node-cores "$CORES" --wall-hours 72

echo "running structural census"
python3 scripts/gtdb_census/census_worker.py \
  --workdir "$WORK" \
  --tasks "$WORK/tasks.jsonl" \
  --syn2b "$SYN2B" \
  --ensure-filename-genome-id \
  --tmpdir /tmp \
  --workers "$CORES" \
  --max-mega 4 \
  --soft-gb 250 --hard-gb 400 --min-free-gb 50 \
  --drain

echo "merging and QC"
python3 scripts/gtdb_census/merge_census.py \
  --workdir "$WORK" --out "$RESULTS" \
  --output-name oapgc_within_species_sv.tsv.gz

echo "OAPGC_LOCAL_CENSUS_DONE"
