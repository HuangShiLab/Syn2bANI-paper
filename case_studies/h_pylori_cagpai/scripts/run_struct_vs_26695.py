#!/usr/bin/env python3
"""Run syn2bani struct --bed for all cohort genomes vs hp26695 locally in parallel."""
import os
import subprocess
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count

S2B = Path('/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani')
GENOMES = Path('/Volumes/MoneyCat/Data/song_2026_hpylori/genomes')
REF = Path('/Users/macstudio/Downloads/Syn2bANI-paper/data/cagpai_pilot/hp26695.fna')
OUT = Path('/Volumes/MoneyCat/Data/song_2026_hpylori/struct_vs_26695_filtered')
OUT.mkdir(parents=True, exist_ok=True)

def run_one(fna):
    gid = fna.stem
    bed = OUT / f"{gid}.vs_hp26695.bed"
    if bed.exists() and bed.stat().st_size > 0:
        return gid, 'skipped'
    cmd = [str(S2B), 'struct', '--bed', '--circular', 'NC_000915.1',
           str(fna), str(REF), '-o', str(bed)]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=300)
        return gid, 'ok'
    except subprocess.TimeoutExpired:
        return gid, 'timeout'
    except subprocess.CalledProcessError as e:
        return gid, f'error:{e.stderr.decode(errors="ignore")[:80]}'

if __name__ == '__main__':
    fnas = sorted(GENOMES.glob('*.fna'))
    print(f"Found {len(fnas)} genomes; using {min(cpu_count(), 8)} workers", file=sys.stderr)
    with Pool(min(cpu_count(), 8)) as pool:
        results = pool.map(run_one, fnas)
    ok = sum(1 for _, s in results if s in ('ok', 'skipped'))
    print(f"Done: {ok}/{len(fnas)} ok or skipped")
    for gid, status in results:
        if status not in ('ok', 'skipped'):
            print(gid, status)
