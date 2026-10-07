#!/usr/bin/env python3
"""Chunked GTDB ANI completion using skani search."""
import argparse
import csv
import gzip
import json
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


def run(cmd, stdout=None, stderr=None):
    with open(stdout, "w") if stdout else open(os.devnull, "w") as fo, \
         open(stderr, "w") if stderr else open(os.devnull, "w") as fe:
        r = subprocess.run(cmd, stdout=fo, stderr=fe, text=True)
    if r.returncode != 0:
        tail = ""
        if stderr and Path(stderr).exists():
            tail = "\n".join(Path(stderr).read_text(errors="replace").splitlines()[-30:])
        raise RuntimeError(f"command failed ({r.returncode}): {' '.join(map(str, cmd))}\n{tail}")


def parse_and_write(rawfile, acc_for_stem, normfh):
    n = 0
    with open(rawfile) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ra = Path(row["Ref_file"]).stem
            qb = Path(row["Query_file"]).stem
            aa = acc_for_stem.get(ra, ra)
            bb = acc_for_stem.get(qb, qb)
            if aa == bb:
                continue
            if bb < aa:
                aa, bb = bb, aa
            normfh.write(f"{aa}\t{bb}\t{row['ANI']}\t{row['Align_fraction_ref']}\t{row['Align_fraction_query']}\n")
            n += 1
    return n


def run_cluster(task):
    skani, cl, accs, manifest, tmpdir, normdir, a = task
    safe = cl.replace("/", "_")
    cdir = tmpdir / safe
    cdir.mkdir(parents=True, exist_ok=True)
    paths = [manifest[x] for x in accs if x in manifest and os.path.isfile(manifest[x])]
    if len(paths) < 2:
        shutil.rmtree(cdir, ignore_errors=True)
        return cl, 0

    sketchdb = cdir / "sketchdb"
    allpaths = cdir / "all_paths.txt"
    allpaths.write_text("\n".join(paths) + "\n")
    run([skani, "sketch", "-t", str(max(1, a.search_threads)), "-l", str(allpaths), "-o", str(sketchdb)])

    chunks = []
    for i in range(0, len(accs), a.chunk_size):
        qacc = accs[i:i + a.chunk_size]
        qpath = cdir / f"query_{i:08d}.txt"
        qpath.write_text("\n".join(manifest[x] for x in qacc) + "\n")
        raw = cdir / f"search_{i:08d}.raw"
        chunks.append((i, qpath, raw, qacc))

    normpath = normdir / f"{safe}.norm"
    acc_for_stem = {Path(p).stem: Path(p).stem for p in paths}
    total = 0
    with open(normpath, "w") as normfh:
        def launch(chunk):
            i, qpath, raw, qacc = chunk
            cmd = [skani, "search", "-t", str(max(1, a.search_threads)),
                   "-d", str(sketchdb), "--ql", str(qpath),
                   "--min-af", str(a.min_af), "-s", str(a.screen),
                   "-o", str(raw)]
            return subprocess.Popen(cmd, stdout=subprocess.DEVNULL,
                                    stderr=open(cdir / f"search_{i:08d}.err", "w"),
                                    text=True)

        running = []
        for chunk in chunks:
            running.append((chunk, launch(chunk)))
            while len(running) >= max(1, a.search_threads // 4):
                time.sleep(15)
                failed = False
                still = []
                for x, p in running:
                    code = p.poll()
                    if code is None:
                        still.append((x, p))
                    elif code != 0:
                        tail = Path(cdir / f"search_{x[0]:08d}.err").read_text(errors="replace")
                        raise RuntimeError(f"skani search failed for {cl} chunk {x[0]}:\n{tail[-4000:]}")
                running = still
                if not running:
                    break
        for x, p in running:
            code = p.wait()
            if code != 0:
                tail = Path(cdir / f"search_{x[0]:08d}.err").read_text(errors="replace")
                raise RuntimeError(f"skani search failed for {cl} chunk {x[0]}:\n{tail[-4000:]}")

        for i, qpath, raw, qacc in chunks:
            total += parse_and_write(raw, acc_for_stem, normfh)

    shutil.rmtree(cdir, ignore_errors=True)
    return cl, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--skani", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--norm-dir", required=True)
    ap.add_argument("--tmpdir", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--search-threads", type=int, default=4)
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--min-af", type=float, default=0.0)
    ap.add_argument("--screen", type=float, default=0.0)
    a = ap.parse_args()

    root = Path(a.workdir)
    out = Path(a.output)
    normdir = Path(a.norm_dir)
    tmpdir = Path(a.tmpdir)
    normdir.mkdir(parents=True, exist_ok=True)
    tmpdir.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)

    manifest = json.load(open(root / "manifest.json"))
    done = {p.name[:-11] for p in (root / "ani").glob("*.ani.tsv.gz")}
    clusters = []
    for p in sorted((root / "cluster_accessions").glob("*.txt")):
        cl = p.name[:-4]
        if cl not in done:
            clusters.append((cl, p.read_text().split()))
    print(f"missing clusters={len(clusters)}", flush=True)

    tasks = [(a.skani, cl, accs, manifest, tmpdir, normdir, a)
             for cl, accs in clusters]
    with mp.Pool(max(1, a.workers)) as pool:
        for cl, n in pool.imap_unordered(run_cluster, tasks):
            print(f"done {cl}: {n} ANI rows", flush=True)

    final_tmp = str(out) + ".unsorted"
    sort_tmp = tmpdir / "sort"
    sort_tmp.mkdir(parents=True, exist_ok=True)
    normfiles = " ".join(f'"{p}"' for p in sorted(normdir.glob("*.norm")))
    cmd = (
        f"{{ printf 'genome_A\\tgenome_B\\tskani_ani\\n'; "
        f"cat {normfiles}; }} | "
        f"LC_ALL=C sort -t $'\\t' -k1,1 -k2,2 -u -S 32G "
        f"--parallel={max(1, min(16, a.search_threads))} -T {sort_tmp} | "
        f"pigz -p {max(1, min(16, a.search_threads))} > {out}"
    )
    subprocess.run(["bash", "-lc", cmd], check=True)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
