#!/usr/bin/env python3
"""ANI pass for the census: skani triangle within every species cluster.

Output: <workdir>/ani/<cluster>.ani.tsv.gz with genome_A, genome_B, skani_ani.
The parser is streaming so clusters with thousands of genomes do not require a
dense Python matrix in memory.
"""
import argparse
import gzip
import json
import os
import subprocess
from pathlib import Path


def _starmap_runner(args):
    return run_cluster(*args)


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

    list_path = outdir / f".{cl_key}.list"
    err_path = outdir / f".{cl_key}.stderr"
    tmp = out.with_name(out.name + ".tmp")
    tmp.unlink(missing_ok=True)
    err_path.unlink(missing_ok=True)
    list_path.write_text("\n".join(paths) + "\n")
    cmd = [skani, "triangle", "-t", "1", "-l", str(list_path)]
    written = 0
    proc = None
    try:
        with err_path.open("w") as errfh, gzip.open(tmp, "wt") as fh:
            fh.write("genome_A\tgenome_B\tskani_ani\n")
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=errfh, text=True, bufsize=1)
            header = proc.stdout.readline()
            if not header:
                raise RuntimeError("skani produced no triangle header")
            n = int(header.strip())

            # skani emits: N; path0; path1 value0; path2 value0 value1; ...
            line = proc.stdout.readline()
            if not line:
                raise RuntimeError("missing genome path 0")
            order = [Path(line.split("\t")[0]).stem]
            acc_for_stem = {Path(manifest[a]).stem: a
                            for a in accs if a in manifest}

            row = 1
            for line in proc.stdout:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                row_stem = Path(fields[0]).stem
                values = fields[1:]
                if len(values) != row:
                    raise RuntimeError(
                        f"expected {row} values on triangle row {row}, got {len(values)}")
                for j, value in enumerate(values):
                    a = acc_for_stem.get(order[j], order[j])
                    b = acc_for_stem.get(row_stem, row_stem)
                    if a == b:
                        continue
                    fh.write(f"{a}\t{b}\t{value}\n")
                    written += 1
                order.append(row_stem)
                row += 1
            if row != n:
                raise RuntimeError(f"expected {n} triangle rows, got {row}")
            ret = proc.wait()
            if ret != 0:
                raise subprocess.CalledProcessError(ret, cmd)
        os.replace(tmp, out)
        return cl_key, written
    except Exception as exc:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass
        if isinstance(exc, subprocess.CalledProcessError):
            tail = err_path.read_text().splitlines()[-20:]
            raise RuntimeError("; ".join(tail)) from exc
        raise
    finally:
        if proc is not None and proc.poll() is None:
            proc.kill()
            proc.wait()
        try:
            list_path.unlink()
        except FileNotFoundError:
            pass
        if out.exists():
            err_path.unlink(missing_ok=True)


def main():
    import multiprocessing as mp
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--skani", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--large-threshold", type=int, default=1000)
    ap.add_argument("--large-workers", type=int, default=2)
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
    jobs.sort(key=lambda x: -len(x[1]))
    large = [(cl, accs) for cl, accs in jobs
             if len(accs) >= args.large_threshold]
    small = [(cl, accs) for cl, accs in jobs
             if len(accs) < args.large_threshold]
    print(f"{len(jobs)} clusters; {len(large)} large; {len(small)} small",
          flush=True)

    large_workers = max(1, min(args.large_workers, args.workers))
    with mp.Pool(large_workers) as pool:
        for cl, n in pool.imap_unordered(_starmap_runner,
                                         task_args(large), chunksize=1):
            if n >= 0:
                print(f"LARGE {cl}: {n} ANI rows", flush=True)
            if n == 0:
                print(f"WARN {cl}: fewer than 2 resolvable genomes", flush=True)

    with mp.Pool(args.workers) as pool:
        for cl, n in pool.imap_unordered(_starmap_runner,
                                         task_args(small), chunksize=4):
            if n >= 0:
                print(f"{cl}: {n} ANI rows", flush=True)
            if n == 0:
                print(f"WARN {cl}: fewer than 2 resolvable genomes", flush=True)


if __name__ == "__main__":
    main()
