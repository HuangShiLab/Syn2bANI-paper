#!/usr/bin/env python3
"""Prepare a Syn2b within-species census for the HROM genome catalog.

HROM already groups conspecific genomes by representative species
(`HROM_Genome_XXXX_N.fna` belongs to cluster `HROM_Genome_XXXX`), so this
script does not infer taxonomy.  It writes the inputs consumed by
``scripts/gtdb_census/census_worker.py`` and ``skani_ani_pass.py``:

- taxonomy.tsv
- cluster_accessions/<cluster>.txt
- manifest.json
- checkpoints/task_counts.json
- tasks.jsonl
- job_plan_summary.md
- resource_estimate.tsv
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

MS_PER_ORDERED_PAIR = 9.0
MS_PER_DIGEST = 42.0


def make_tasks(cluster_accessions, batch_size):
    tasks = []
    unique_pairs = 0
    for cl, accs in sorted(cluster_accessions.items()):
        n = len(accs)
        unique_pairs += n * (n - 1) // 2
        nb = (n + batch_size - 1) // batch_size
        for bi in range(nb):
            for bj in range(bi, nb):
                ni = min(batch_size, n - bi * batch_size)
                nj = min(batch_size, n - bj * batch_size)
                # The tool cost scales with the number of genomes in the two
                # blocks; diagonal blocks deliberately include self checks and
                # are filtered when the matrix is written.
                est_core_h = (ni + nj) ** 2 * MS_PER_ORDERED_PAIR / 3_600_000
                tasks.append({
                    "cluster": cl,
                    "n_genomes": n,
                    "block_i": bi,
                    "block_j": bj,
                    "block_size": batch_size,
                    "est_core_h": round(est_core_h, 6),
                })
    tasks.sort(key=lambda t: (-t["est_core_h"], t["cluster"],
                              t["block_i"], t["block_j"]))
    return tasks, unique_pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", required=True,
                    help="HROM_Conspecific-genomes-metadata.tsv")
    ap.add_argument("--genome-root", required=True,
                    help="HROM_nonredundant_genomes directory")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--batch-size", type=int, default=250)
    ap.add_argument("--node-cores", type=int, default=16)
    ap.add_argument("--wall-hours", type=float, default=24.0)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--verify-paths", action="store_true",
                    help="stat every expected FASTA (slow on Lustre)")
    ap.add_argument("--limit-clusters", type=int, default=0,
                    help="smoke test: retain the N largest clusters only")
    args = ap.parse_args()

    genome_root = Path(args.genome_root)
    root = Path(args.workdir)
    (root / "cluster_accessions").mkdir(parents=True, exist_ok=True)
    (root / "checkpoints").mkdir(parents=True, exist_ok=True)

    # cluster -> accession -> (bp, quality)
    clusters = {}
    missing = []
    rows = 0
    with open(args.metadata, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # header has one trailing blank name
        for row in reader:
            if not row:
                continue
            genome, cluster = row[0], row[-1]
            if not cluster:
                raise ValueError(f"missing HROM species cluster for {genome}")
            path = genome_root / cluster / f"{genome}.fna"
            if args.verify_paths and not path.is_file():
                missing.append(genome)
                continue
            clusters.setdefault(cluster, {})[genome] = (int(row[4]), row[8])
            rows += 1

    if args.limit_clusters > 0:
        keep = set(sorted(clusters,
                          key=lambda c: (-len(clusters[c]), c))[:args.limit_clusters])
        clusters = {c: clusters[c] for c in keep}
    multi = {c: dict(a) for c, a in clusters.items() if len(a) >= 2}
    accession_lists = {c: sorted(a) for c, a in multi.items()}
    retained_rows = sum(len(a) for a in multi.values())

    manifest = {
        acc: str(genome_root / cl / f"{acc}.fna")
        for cl, accs in accession_lists.items() for acc in accs
    }
    for cl, accs in accession_lists.items():
        (root / "cluster_accessions" / f"{cl}.txt").write_text(
            "\n".join(accs) + "\n")

    tasks, unique_pairs = make_tasks(accession_lists, args.batch_size)
    with (root / "tasks.jsonl").open("w") as fh:
        for task in tasks:
            fh.write(json.dumps(task) + "\n")
    counts = Counter(task["cluster"] for task in tasks)
    with (root / "checkpoints" / "task_counts.json").open("w") as fh:
        for cl in sorted(counts):
            fh.write(f"{cl}\t{counts[cl]}\n")

    with (root / "taxonomy.tsv").open("w") as tax, \
         (root / "hrom_genome_metadata.tsv").open("w") as meta:
        tax.write("accession\ttaxonomy\n")
        meta.write("genome\tcluster\tlength_bp\tquality\n")
        for cl, accs in multi.items():
            for acc, (bp, quality) in sorted(accs.items()):
                tax.write(f"{acc}\ts__{cl}\n")
                meta.write(f"{acc}\t{cl}\t{bp}\t{quality}\n")

    (root / "manifest.json").write_text(json.dumps(manifest, indent=2))
    (root / "missing_genomes.txt").write_text("\n".join(missing) + "\n")

    selected_bp = sum(bp for accs in multi.values() for bp, _ in accs.values())
    structure_h = sum(task["est_core_h"] for task in tasks)
    digest_h = len(manifest) * MS_PER_DIGEST / 3_600_000
    capacity_h = args.jobs * args.node_cores * args.wall_hours

    summary = [
        "# HROM within-species census plan",
        "",
        f"- metadata rows retained: {retained_rows:,}",
        f"- multi-genome clusters: {len(multi):,}",
        f"- genomes in census: {len(manifest):,}",
        f"- unique unordered pairs: {unique_pairs:,}",
        f"- selected sequence bp: {selected_bp:,} ({selected_bp / 1e9:.1f} GB)",
        f"- tasks at block size {args.batch_size}: {len(tasks):,}",
        f"- structural estimate: {structure_h:,.0f} core-h "
        f"({MS_PER_ORDERED_PAIR:g} ms per ordered-pair equivalent)",
        f"- digest estimate: {digest_h:,.0f} core-h",
        f"- requested capacity: {capacity_h:,.0f} core-h "
        f"({args.jobs} jobs x {args.node_cores} cores x {args.wall_hours:g} h)",
        f"- missing FASTAs (--verify-paths only): {len(missing):,}",
        "",
        "## Largest clusters",
        "",
        "| cluster | genomes | unordered pairs |",
        "|---|---:|---:|",
    ]
    for cl in sorted(accession_lists,
                     key=lambda c: (-len(accession_lists[c]), c))[:20]:
        n = len(accession_lists[cl])
        summary.append(f"| {cl} | {n:,} | {n * (n - 1) // 2:,} |")
    (root / "job_plan_summary.md").write_text("\n".join(summary) + "\n")

    with (root / "resource_estimate.tsv").open("w") as fh:
        fh.write("metric\tvalue\n")
        for key, value in [
            ("metadata_rows_retained", retained_rows),
            ("multi_genome_clusters", len(multi)),
            ("census_genomes", len(manifest)),
            ("unique_pairs", unique_pairs),
            ("selected_bp", selected_bp),
            ("tasks", len(tasks)),
            ("structure_core_h_est", round(structure_h, 2)),
            ("digest_core_h_est", round(digest_h, 2)),
            ("requested_capacity_core_h", capacity_h),
        ]:
            fh.write(f"{key}\t{value}\n")

    print(f"clusters={len(multi):,} genomes={len(manifest):,} "
          f"pairs={unique_pairs:,} tasks={len(tasks):,} "
          f"structure_est={structure_h:,.0f} core-h -> {root}")


if __name__ == "__main__":
    main()
