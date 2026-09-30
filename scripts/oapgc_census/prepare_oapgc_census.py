#!/usr/bin/env python3
"""Prepare an OAPGC within-species structural census.

The public cluster annotation has 149,921 HQMAG rows.  A non-null
``Represent(STRAINs)`` value marks the 99,215 strain-level non-redundant
genomes used by OAPGC.  This script keeps those representative rows and groups
them by SGB (species-level cluster).  It emits the same inputs consumed by
``scripts/gtdb_census/census_worker.py`` and ``skani_ani_pass.py``.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import time
from collections import Counter
from pathlib import Path

MS_PER_ORDERED_PAIR = 9.0       # measured HROM/GTDB planning constant
MS_PER_DIGEST = 0.042           # 42 ms per genome
GZIPPED_BYTES_PER_PAIR = 30.0   # empirical HROM structural output


def make_tasks(cluster_accessions: dict[str, list[str]], batch_size: int):
    tasks = []
    unique_pairs = 0
    for cluster, accessions in sorted(cluster_accessions.items()):
        n = len(accessions)
        unique_pairs += n * (n - 1) // 2
        nblocks = (n + batch_size - 1) // batch_size
        for block_i in range(nblocks):
            for block_j in range(block_i, nblocks):
                ni = min(batch_size, n - block_i * batch_size)
                nj = min(batch_size, n - block_j * batch_size)
                # Diagonal blocks include ordered self-checks, hence (ni+nj)^2
                # is a conservative but empirical planning cost.
                est_core_h = (ni + nj) ** 2 * MS_PER_ORDERED_PAIR / 3_600_000
                tasks.append({
                    "cluster": cluster,
                    "n_genomes": n,
                    "block_i": block_i,
                    "block_j": block_j,
                    "block_size": batch_size,
                    "est_core_h": round(est_core_h, 6),
                })
    tasks.sort(key=lambda t: (-t["est_core_h"], t["cluster"],
                              t["block_i"], t["block_j"]))
    return tasks, unique_pairs


def human_gb(nbytes: float) -> str:
    return f"{nbytes / 1e9:.2f} GB"


def index_fasta_files(genome_root: Path):
    by_stem = {}
    n_files = 0
    for dirpath, dirnames, filenames in os.walk(genome_root):
        dirnames.sort()
        for filename in sorted(filenames):
            if Path(filename).suffix.lower() not in {".fna", ".fa", ".fasta"}:
                continue
            path = Path(dirpath) / filename
            stem = path.stem
            n_files += 1
            if stem in by_stem:
                raise ValueError(f"duplicate FASTA stem {stem}: "
                                 f"{by_stem[stem]} and {path}")
            by_stem[stem] = str(path)
    return by_stem, n_files


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotation", required=True,
                    help="OAPGC cluster.annotation.xlsx")
    ap.add_argument("--genome-root",
                    help="directory containing extracted OAPGC FASTAs; if "
                         "omitted, write plans but not manifest.json")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--batch-size", type=int, default=250)
    ap.add_argument("--node-cores", type=int, default=16)
    ap.add_argument("--wall-hours", type=float, default=48.0)
    ap.add_argument("--limit-clusters", type=int, default=0,
                    help="smoke test: retain N largest multi-genome clusters")
    ap.add_argument("--max-cluster-size", type=int, default=0,
                    help="optional sensitivity cap; 0 keeps all genomes")
    args = ap.parse_args()

    from openpyxl import load_workbook

    annotation = Path(args.annotation)
    if not annotation.is_file():
        raise FileNotFoundError(annotation)
    root = Path(args.workdir)
    (root / "cluster_accessions").mkdir(parents=True, exist_ok=True)
    (root / "checkpoints").mkdir(parents=True, exist_ok=True)
    (root / "manifests").mkdir(parents=True, exist_ok=True)

    wb = load_workbook(annotation, read_only=True, data_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = ws.iter_rows(values_only=True)
    header = [str(x).strip() if x is not None else "" for x in next(rows)]
    try:
        i_genome = header.index("HQMAGs")
        i_sgb = header.index("SGBs")
        i_strain_rep = next(i for i, x in enumerate(header)
                            if x.replace(" ", "").lower() == "represent(strains)")
    except (ValueError, StopIteration) as e:
        raise ValueError(f"unexpected annotation columns: {header}") from e

    selected = {}
    all_rows = 0
    for values in rows:
        if not values or i_genome >= len(values):
            continue
        genome = values[i_genome]
        sgb = values[i_sgb]
        strain_rep = values[i_strain_rep] if i_strain_rep < len(values) else None
        if genome is None or sgb is None or strain_rep is None:
            continue
        genome, sgb, strain_rep = str(genome).strip(), str(sgb).strip(), str(strain_rep).strip()
        if not genome or not sgb or genome != strain_rep:
            continue
        selected.setdefault(sgb, []).append(genome)
        all_rows += 1

    multi = {cluster: sorted(set(genomes))
             for cluster, genomes in selected.items()
             if len(set(genomes)) >= 2}
    if args.max_cluster_size and args.max_cluster_size > 0:
        multi = {cluster: accessions[:args.max_cluster_size]
                 for cluster, accessions in multi.items()}
    if args.limit_clusters and args.limit_clusters > 0:
        keep = set(sorted(multi, key=lambda c: (-len(multi[c]), c))
                   [:args.limit_clusters])
        multi = {cluster: multi[cluster] for cluster in sorted(keep)}

    fasta_index = None
    missing = []
    if args.genome_root:
        fasta_index, n_fasta = index_fasta_files(Path(args.genome_root))
        for accessions in multi.values():
            missing.extend(acc for acc in accessions if acc not in fasta_index)
        if missing:
            example = "\n".join(missing[:20])
            raise FileNotFoundError(
                f"{len(missing)} representative FASTAs absent; examples:\n{example}")
        manifest = {acc: fasta_index[acc]
                    for accessions in multi.values() for acc in accessions}
        with (root / "manifest.json").open("w") as fh:
            json.dump(manifest, fh, indent=2, sort_keys=True)
        with (root / "manifests" / "selected_manifest.tsv").open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["genome_id", "cluster", "fasta_path"])
            for cluster, accessions in sorted(multi.items()):
                for acc in accessions:
                    w.writerow([acc, cluster, fasta_index[acc]])

    accession_lists = {cluster: sorted(accessions)
                       for cluster, accessions in multi.items()}
    for cluster, accessions in accession_lists.items():
        safe = cluster.replace("/", "_")
        (root / "cluster_accessions" / f"{safe}.txt").write_text(
            "\n".join(accessions) + "\n")

    tasks, unique_pairs = make_tasks(accession_lists, args.batch_size)
    with (root / "tasks.jsonl").open("w") as fh:
        for task in tasks:
            fh.write(json.dumps(task, sort_keys=True) + "\n")

    counts = Counter(task["cluster"] for task in tasks)
    with (root / "checkpoints" / "task_counts.json").open("w") as fh:
        for cluster in sorted(counts):
            fh.write(f"{cluster}\t{counts[cluster]}\n")

    retained = sum(len(x) for x in accession_lists.values())
    with (root / "genome_metadata.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome_id", "cluster", "in_census", "fasta_path"])
        fasta_lookup = fasta_index or {}
        for cluster, accessions in sorted(selected.items()):
            census_set = set(accession_lists.get(cluster, []))
            for acc in sorted(set(accessions)):
                w.writerow([acc, cluster, int(acc in census_set),
                            fasta_lookup.get(acc, "")])

    # Output the flat list consumed by extract_oapgc_representatives.py.
    with (root / "representative_genomes.txt").open("w") as fh:
        for cluster, accessions in sorted(accession_lists.items()):
            for acc in accessions:
                fh.write(f"{acc}\t{cluster}\n")

    core_h = sum(t["est_core_h"] for t in tasks)
    digest_h = retained * MS_PER_DIGEST / 3_600_000
    node_h = (core_h + digest_h) / max(1, args.node_cores)
    output_gb = unique_pairs * GZIPPED_BYTES_PER_PAIR / 1e9
    cluster_sizes = sorted((len(x) for x in accession_lists.values()), reverse=True)

    with (root / "resource_estimate.tsv").open("w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["annotation_hqmag_rows", 149921])
        w.writerow(["strain_representative_rows", all_rows])
        w.writerow(["multi_genome_clusters", len(accession_lists)])
        w.writerow(["genomes_in_census", retained])
        w.writerow(["unique_within_species_pairs", unique_pairs])
        w.writerow(["structural_tasks", len(tasks)])
        w.writerow(["largest_cluster", cluster_sizes[0] if cluster_sizes else 0])
        w.writerow(["estimated_syn2b_core_hours", round(core_h, 2)])
        w.writerow(["estimated_digest_core_hours", round(digest_h, 4)])
        w.writerow(["node_cores", args.node_cores])
        w.writerow(["estimated_node_hours", round(node_h, 2)])
        w.writerow(["estimated_output_gzip_gb", round(output_gb, 2)])
        w.writerow(["recommended_wall_hours", args.wall_hours])

    top = "\n".join(f"| {cluster} | {len(accession_lists[cluster])} |"
                    for cluster in sorted(accession_lists,
                                          key=lambda c: (-len(accession_lists[c]), c))[:20])
    with (root / "job_plan_summary.md").open("w") as fh:
        fh.write(f"""# OAPGC within-species census plan

Generated: `{time.strftime('%F %T')}`

## Selection

- annotation rows: **{149921:,}** HQMAGs
- non-redundant strain representatives: **{all_rows:,}**
- multi-genome SGB clusters: **{len(accession_lists):,}**
- genomes in structural census: **{retained:,}**
- unique within-species pairs: **{unique_pairs:,}**

## Cost

- batch size: **{args.batch_size}**
- structural tasks: **{len(tasks):,}**
- Syn2b estimate: **{core_h:,.1f}** core-h
- digestion estimate: **{digest_h:,.3f}** core-h
- at {args.node_cores} local cores: **{node_h:,.1f}** node-h
- gzipped structural output estimate: **{output_gb:,.2f} GB**

This estimate assumes OAPGC MAG lengths and assembly fragmentation are broadly
similar to HROM. Run the first 20 largest clusters as a smoke test before the
full run; revise constants if measured core-h/pair differs.

## Largest clusters

| SGB | genomes |
|---|---:|
{top}
""")

    print(json.dumps({
        "workdir": str(root),
        "annotation_rows": 149921,
        "strain_representatives": all_rows,
        "clusters": len(accession_lists),
        "genomes": retained,
        "pairs": unique_pairs,
        "tasks": len(tasks),
        "estimated_core_h": round(core_h + digest_h, 2),
        "estimated_node_h": round(node_h, 2),
        "estimated_output_gzip_gb": round(output_gb, 2),
        "missing_fastas": len(missing),
    }, indent=2))


if __name__ == "__main__":
    main()
