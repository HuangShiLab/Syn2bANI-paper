#!/usr/bin/env python3
"""Run syn2bani ani --strata-out on a specific pair list, grouping by reference.

syn2bani ani --ql/--rl computes the Cartesian product, so this script groups
pairs by reference and runs one ani invocation per reference with exactly the
queries that should be compared to it.

Usage as SLURM array (example):
    # 1. Split pair list into chunks
    python3 run_strata_by_reference.py split \
        results/matrix_gtdb_r207_100k_v8_final.tsv \
        --genome-dir /lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all \
        --outdir results/strata_ref_chunks \
        --n-chunks 25

    # 2. Submit array (1-25)
    sbatch scripts/strata_by_reference_array.slurm

    # 3. Merge
    python3 run_strata_by_reference.py merge \
        results/strata_ref_chunks \
        --strata-out results/gtdb_r207_100k_strata.tsv \
        --ani-out results/gtdb_r207_100k_11enzyme.tsv
"""
import argparse
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd


def split_pairs(pairs_tsv, outdir, genome_dir, n_chunks, suffix=".fna"):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(pairs_tsv, sep="\t")
    n = len(m)
    chunk_size = (n + n_chunks - 1) // n_chunks
    for i in range(n_chunks):
        s = i * chunk_size
        e = min(s + chunk_size, n)
        chunk = m.iloc[s:e].copy()
        chunk["query_path"] = chunk["query"].apply(lambda x: os.path.join(genome_dir, x + suffix))
        chunk["ref_path"] = chunk["reference"].apply(lambda x: os.path.join(genome_dir, x + suffix))
        chunk.to_csv(outdir / f"pairs_{i+1:03d}.tsv", sep="\t", index=False)
        print(f"chunk {i+1}: {len(chunk)} pairs")


def run_chunk(chunk_tsv, syn2bani, panel, threads_per_run, max_workers, out_strata, out_ani, calibrate=False):
    df = pd.read_csv(chunk_tsv, sep="\t")
    by_ref = defaultdict(list)
    for _, row in df.iterrows():
        by_ref[row["ref_path"]].append((row["query_path"], row["query"]))

    chunk_tag = Path(chunk_tsv).stem
    tmp_dir = Path(out_strata).parent / ".tmp" / chunk_tag
    tmp_dir.mkdir(parents=True, exist_ok=True)
    do_calibrate = calibrate

    def run_one(ref_path, items):
        # Write temporary query list and reference list for this reference
        q_paths = [x[0] for x in items]
        q_names = [x[1] for x in items]
        tag = os.path.basename(ref_path).replace(".fna", "")
        ql_file = tmp_dir / f"ql_{tag}.txt"
        rl_file = tmp_dir / f"rl_{tag}.txt"
        ref_ani = tmp_dir / f"ani_{tag}.tsv"
        ref_strata = tmp_dir / f"strata_{tag}.tsv"
        with open(ql_file, "w") as fh:
            fh.write("\n".join(q_paths) + "\n")
        with open(rl_file, "w") as fh:
            fh.write(ref_path + "\n")
        cmd = [
            syn2bani, "ani",
            "--ql", str(ql_file),
            "--rl", str(rl_file),
            "-e", panel,
            "-o", str(ref_ani),
            "--strata-out", str(ref_strata),
            "-p", "-t", str(threads_per_run),
        ]
        if do_calibrate:
            cmd.append("--calibrate")
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            return (ref_path, len(items), None, ref_ani, ref_strata)
        except subprocess.CalledProcessError as e:
            return (ref_path, len(items), e.stderr[:500], None, None)

    # Run ani invocations in parallel within this task
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as exe:
        futures = {exe.submit(run_one, ref, items): ref for ref, items in by_ref.items()}
        for fut in as_completed(futures):
            results.append(fut.result())

    # Merge per-reference outputs into chunk output files
    with open(out_strata, "w") as sfh, open(out_ani, "w") as afh:
        strata_header = "query\treference\tenzyme\ttag_len\tbody_len\tn_miss\thist\n"
        ani_header = "query\treference\tani\tani_uniform\taf_query\taf_reference\tstd_err\n"
        sfh.write(strata_header)
        afh.write(ani_header)
        for ref_path, nq, err, ref_ani, ref_strata in results:
            if err:
                print(f"ERROR {ref_path}: {err}", file=sys.stderr)
                continue
            print(f"OK {ref_path}: {nq} queries")
            if ref_strata and ref_strata.exists():
                with open(ref_strata) as fh:
                    _ = fh.readline()  # skip header
                    for line in fh:
                        sfh.write(line)
            if ref_ani and ref_ani.exists():
                with open(ref_ani) as fh:
                    _ = fh.readline()  # skip header
                    for line in fh:
                        afh.write(line)

    success = sum(1 for _, _, e, _, _ in results if e is None)
    print(f"chunk {chunk_tsv}: {success}/{len(results)} references OK")


def merge_chunks(chunk_dir, strata_out, ani_out, prefix=""):
    chunk_dir = Path(chunk_dir)
    strata_files = sorted(chunk_dir.glob(f"{prefix}strata_*.tsv"))
    ani_files = sorted(chunk_dir.glob(f"{prefix}ani_*.tsv"))

    if not strata_files or not ani_files:
        print("No chunk outputs found", file=sys.stderr)
        return 1

    # Merge strata
    with open(strata_out, "w") as out:
        header_written = False
        for f in strata_files:
            with open(f) as fh:
                header = fh.readline()
                if not header_written:
                    out.write(header)
                    header_written = True
                for line in fh:
                    out.write(line)
    print(f"merged strata -> {strata_out}")

    # Merge ani
    with open(ani_out, "w") as out:
        header_written = False
        for f in ani_files:
            with open(f) as fh:
                header = fh.readline()
                if not header_written:
                    out.write(header)
                    header_written = True
                for line in fh:
                    out.write(line)
    print(f"merged ani -> {ani_out}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("split", help="split pair TSV into chunks grouped by reference")
    sp.add_argument("pairs_tsv")
    sp.add_argument("--outdir", default="results/strata_ref_chunks")
    sp.add_argument("--genome-dir", required=True)
    sp.add_argument("--n-chunks", type=int, default=25)
    sp.add_argument("--suffix", default=".fna")

    rp = sub.add_parser("run", help="run one chunk")
    rp.add_argument("chunk_tsv")
    rp.add_argument("--syn2bani", default="/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani")
    rp.add_argument("--panel", default="BcgI,AlfI,AloI,FalI,BplI,Bsp24I,PpiI,PsrI,BsaXI,CjeI,CjePI")
    rp.add_argument("--threads-per-run", type=int, default=4)
    rp.add_argument("--max-workers", type=int, default=8)
    rp.add_argument("--strata-out", required=True)
    rp.add_argument("--ani-out", required=True)
    rp.add_argument("--calibrate", action="store_true", help="pass --calibrate to syn2bani ani")

    mp = sub.add_parser("merge", help="merge chunk outputs")
    mp.add_argument("chunk_dir")
    mp.add_argument("--strata-out", default="results/gtdb_r207_100k_strata.tsv")
    mp.add_argument("--ani-out", default="results/gtdb_r207_100k_11enzyme.tsv")
    mp.add_argument("--prefix", default="",
                    help="filename prefix to select chunks, e.g. 'cal_'")

    args = ap.parse_args()

    if args.cmd == "split":
        split_pairs(args.pairs_tsv, args.outdir, args.genome_dir, args.n_chunks, args.suffix)
    elif args.cmd == "run":
        run_chunk(args.chunk_tsv, args.syn2bani, args.panel, args.threads_per_run,
                  args.max_workers, args.strata_out, args.ani_out, args.calibrate)
    elif args.cmd == "merge":
        return merge_chunks(args.chunk_dir, args.strata_out, args.ani_out, args.prefix)
    return 0


if __name__ == "__main__":
    sys.exit(main())
