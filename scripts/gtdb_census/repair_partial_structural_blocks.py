#!/usr/bin/env python3
"""Plan or archive-and-retry partial GTDB structural block-pair tasks.

The input is a block-pair coverage table produced by
``diagnose_cluster_blockpairs.py``. In plan mode this writes the original
census-task JSON schema for every partial block pair. In execute mode it also
moves the existing partial output and its ``done`` checkpoint into a
quarantine directory and removes any stale lock, so the standard census worker
can safely rerun exactly those tasks.
"""

import argparse
import json
import shutil
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--cluster", required=True)
    ap.add_argument("--workdir", default="")
    ap.add_argument("--output", required=True)
    ap.add_argument("--archive-dir", default="")
    ap.add_argument("--execute", action="store_true")
    args = ap.parse_args()

    tasks = []
    with open(args.coverage, newline="") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        required = {"block_i", "block_j", "expected_pairs", "observed_pairs", "status"}
        if not required.issubset(ix):
            raise SystemExit(f"missing coverage columns: {sorted(required - set(ix))}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header) or f[ix["status"]] != "partial":
                continue
            expected = int(f[ix["expected_pairs"]])
            observed = int(f[ix["observed_pairs"]])
            if observed >= expected:
                continue
            tasks.append({
                "cluster": args.cluster,
                "block_i": int(f[ix["block_i"]]),
                "block_j": int(f[ix["block_j"]]),
                "block_size": int(f.get(ix["block_size"], 500)) if "block_size" in ix else 500,
                "est_core_h": 0.0,
                "repair_missing_pairs": expected - observed,
                "repair_expected_pairs": expected,
                "repair_observed_pairs": observed,
            })

    tmp = Path(str(args.output) + ".tmp")
    with tmp.open("w") as out:
        for task in tasks:
            out.write(json.dumps(task, sort_keys=True) + "\n")
    tmp.replace(args.output)

    archived = []
    if args.execute:
        if not args.workdir or not args.archive_dir:
            raise SystemExit("--execute requires --workdir and --archive-dir")
        root = Path(args.workdir)
        archive = Path(args.archive_dir)
        archive.mkdir(parents=True)
        for task in tasks:
            task_id = f"{task['cluster'].replace('/', '_')}|{task['block_i']}x{task['block_j']}"
            output = root / "outputs" / f"{task_id}.tsv.gz"
            done = root / "checkpoints" / f"{task_id}.done"
            lock = root / "locks" / task_id
            for src in (output, done):
                if src.exists():
                    dst = archive / src.name
                    if dst.exists():
                        raise SystemExit(f"refusing to overwrite archived file: {dst}")
                    shutil.move(str(src), str(dst))
                    archived.append(str(dst))
            try:
                lock.rmdir()
            except FileNotFoundError:
                pass
            except OSError as exc:
                raise SystemExit(f"could not remove stale lock {lock}: {exc}") from exc

    print(json.dumps({
        "cluster": args.cluster,
        "partial_tasks": len(tasks),
        "missing_pairs": sum(x["repair_missing_pairs"] for x in tasks),
        "output": str(args.output),
        "executed": args.execute,
        "archived_files": len(archived),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
