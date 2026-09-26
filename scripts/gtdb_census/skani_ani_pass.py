#!/usr/bin/env python3
"""ANI pass for the census: skani triangle within every species cluster.

Complements census_worker.py (which computes the structural columns) with the
ANI column. Idempotent: clusters with an existing output are skipped, so the
pass can be resubmitted freely and can run concurrently with the census array.

Output: <workdir>/ani/<cluster>.ani.tsv.gz  with rows  a  b  ani
Run AFTER prepare_workdir.py (needs cluster_accessions/ and manifest.json).
"""
import argparse
import gzip
import json
import os
import subprocess
from pathlib import Path


def run_cluster(skani, cl_key, accs, manifest, outdir):
    paths = []
    for acc in accs:
        p = manifest.get(acc)
        if p:
            paths.append(p)
    if len(paths) < 2:
        return cl_key, 0
    out = outdir / f"{cl_key}.ani.tsv.gz"
    if out.exists():
        return cl_key, -1
    res = subprocess.run([skani, "triangle"] + paths, check=True,
                         capture_output=True, text=True)
    order, ani = [], {}
    for i, line in enumerate(res.stdout.splitlines()):
        if i == 0 or not line.strip():
            continue
        f = line.split("\t")
        stem = Path(f[0]).stem
        order.append(stem)
        for j, v in enumerate(f[1:]):
            ani[(order[j], stem)] = v
    stem2acc = {Path(manifest[a]).stem: a for a in accs if a in manifest}
    n = 0
    tmp = out.with_name(out.name + ".tmp")
    with gzip.open(tmp, "wt") as fh:
        fh.write("genome_A\tgenome_B\tskani_ani\n")
        for (sa, sb), v in ani.items():
            fh.write(f"{stem2acc.get(sa, sa)}\t{stem2acc.get(sb, sb)}\t{v}\n")
            n += 1
    # A killed worker can leave a partial gzip; publication must be atomic so
    # an idempotent rerun does not mistake it for a completed cluster.
    os.replace(tmp, out)
    return cl_key, n


def main():
    import multiprocessing as mp
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--skani", required=True)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    root = Path(args.workdir)
    outdir = root / "ani"
    outdir.mkdir(exist_ok=True)
    manifest = json.loads((root / "manifest.json").read_text())
    jobs = []
    for f in sorted((root / "cluster_accessions").glob("*.txt")):
        cl = f.name[:-4]
        accs = f.read_text().split()
        if len(accs) >= 2:
            jobs.append((cl, accs))
    print(f"{len(jobs)} clusters")
    with mp.Pool(args.workers) as pool:
        for cl, n in pool.starmap(run_cluster,
                                  [(args.skani, cl, accs, manifest, outdir)
                                   for cl, accs in jobs]):
            if n >= 0:
                print(f"{cl}: {n} ANI rows", flush=True)
            if n == 0:
                print(f"WARN {cl}: fewer than 2 resolvable genomes", flush=True)


if __name__ == "__main__":
    main()
