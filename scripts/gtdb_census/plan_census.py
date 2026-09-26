#!/usr/bin/env python3
"""Build the GTDB within-species census task plan.

Reads the GTDB R207 taxonomy table, enumerates species clusters, splits each
cluster into query batches, and emits one JSON task per (cluster, query batch)
plus a summary sized for a FEW LONG SLURM JOBS (not many short ones).

A task's worker builds a temp TGT directory containing the query batch plus the
whole cluster and keeps only rows whose genome_A is in the query batch, so
tasks partition the output rows exactly and any task can be re-run safely.

Outputs (in --outdir):
  tasks.jsonl        one JSON task per line, cost-descending
  job_plan_summary.md  core-hour budget split across --jobs long jobs
"""
import argparse
import gzip
import json
from pathlib import Path

MS_PER_ORDERED_PAIR = 9.0     # measured: 4.36 s / 484 ordered pairs (22-genome bench)
MS_PER_DIGEST = 42.0          # measured: E. coli K-12, 4-enzyme panel
CORE_H = 3600.0 * 1000.0      # ms per core-hour


def load_clusters(taxonomy_tsv):
    clusters = {}
    op = gzip.open if taxonomy_tsv.endswith(".gz") else open
    with op(taxonomy_tsv, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            acc, tax = parts[0], parts[1]
            if acc.startswith(("GB_", "RS_")):
                acc = acc[3:]
            sp = next((p[3:] for p in tax.split(";") if p.startswith("s__")), "")
            if sp:
                clusters.setdefault(sp, []).append(acc)
    for accs in clusters.values():
        accs.sort()
    return clusters


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--taxonomy", required=True, help="accession_taxonomy_r207.tsv.gz")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--batch-size", type=int, default=500,
                    help="query genomes per task (block-pair tasks: dir = block_i+block_j)")
    ap.add_argument("--jobs", type=int, default=6,
                    help="number of long jobs the budget is split across in the summary")
    ap.add_argument("--node-cores", type=int, default=16)
    ap.add_argument("--wall-hours", type=float, default=96.0)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    clusters = load_clusters(args.taxonomy)

    tasks = []
    digest_h = 0.0
    n_genomes_multi = 0
    for sp, accs in clusters.items():
        n = len(accs)
        if n < 2:
            continue
        n_genomes_multi += n
        digest_h += n * MS_PER_DIGEST / CORE_H
        # block-pair tasks: dir = block_i + block_j, so the tool computes
        # exactly (i+j) size-ordered pairs -- the whole-cluster directory of
        # the earlier query-batch design made every task pay n**2, inflating
        # mega clusters by n/batch (E. coli ~134x).
        nb = (n + args.batch_size - 1) // args.batch_size
        for bi in range(nb):
            for bj in range(bi, nb):
                ni = min(args.batch_size, n - bi * args.batch_size)
                nj = min(args.batch_size, n - bj * args.batch_size)
                size = ni + nj
                cost_h = size * size * MS_PER_ORDERED_PAIR / CORE_H
                tasks.append({
                    "cluster": sp, "n_genomes": n,
                    "block_i": bi, "block_j": bj, "block_size": args.batch_size,
                    "est_core_h": round(cost_h, 4),
                })

    tasks.sort(key=lambda t: -t["est_core_h"])
    with open(out / "tasks.jsonl", "w") as fh:
        for t in tasks:
            fh.write(json.dumps(t) + "\n")

    total_h = sum(t["est_core_h"] for t in tasks) + digest_h
    job_cap = args.node_cores * args.wall_hours
    n_jobs_needed = max(1, int(-(-total_h // job_cap)))

    big = tasks[:8]
    lines = [
        "# GTDB within-species census — task plan",
        "",
        f"- clusters with >=2 genomes: {sum(1 for a in clusters.values() if len(a) >= 2):,}",
        f"- genomes in those clusters: {n_genomes_multi:,}",
        f"- tasks (query batches of {args.batch_size}): {len(tasks):,}",
        f"- digest budget: {digest_h:,.0f} core-hours",
        f"- synteny budget: {sum(t['est_core_h'] for t in tasks):,.0f} core-hours",
        f"- total: {total_h:,.0f} core-hours  (at {MS_PER_ORDERED_PAIR:g} ms/ordered pair)",
        f"- capacity per long job: {args.node_cores} cores x {args.wall_hours:g} h = {job_cap:,.0f} core-hours",
        f"- long jobs needed (including ~20% I/O headroom): {max(args.jobs, int(-(-total_h * 1.2 // job_cap)))}",
        "",
        "## Largest tasks (schedule first; workers claim cost-descending anyway)",
        "",
        "| cluster | genomes | task batch | est core-h |",
        "|---|---:|---|---:|",
    ]
    for t in big:
        lines.append(f"| {t['cluster']} | {t['n_genomes']:,} | "
                     f"blocks {t['block_i']}x{t['block_j']} | {t['est_core_h']:.1f} |")
    lines += [
        "",
        "Execution policy: submit a FEW LONG array jobs (e.g. `--array=0-%d%%%d`)," % (max(args.jobs, n_jobs_needed) - 1, max(args.jobs, n_jobs_needed)),
        "each running one `census_worker.py` with `--workers = node cores`. Tasks are",
        "claimed atomically via lock directories on the shared filesystem, so jobs can",
        "start at different times, die and be resubmitted without double work. A space",
        "auditor inside each worker watches quota during the run (see README).",
        "",
    ]
    (out / "job_plan_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
