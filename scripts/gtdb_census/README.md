# GTDB within-species structural census

All within-species pairs of GTDB R207 (711,020,841 pairs across 317,542
genomes; see `results/census/job_plan_summary.md`), computed with the Rust
Syn2b binary. Design constraints from the lab HPC policy:

- **A few long jobs, not many short ones.** One small SLURM array
  (`census_long_job.slurm`, 6 x 16 cores x 96 h) is the entire submission.
  Work is handed out by atomic task claims on the shared filesystem, so jobs
  can start late, die and be resubmitted without lost or duplicated work.
- **Space is audited during the run.** A background auditor in every job
  measures the workdir (`du`) and filesystem free space every 5 minutes.
  Above the soft limit it reclaims TGTs of finished clusters (re-digest costs
  42 ms/genome); at the hard limit — or below the free-space floor — workers
  stop claiming tasks until space recovers. All thresholds are CLI flags.

## Pipeline

```
plan_census.py      taxonomy -> tasks.jsonl (query batches, cost-descending)
                    + job_plan_summary.md
prepare_workdir.py  on the HPC: cluster accession lists, manifest.json
                    (accession -> FASTA), task counts, missing-genome report
census_worker.py    the long job: N processes claim tasks, digest missing
                    TGTs, run `syn2b synteny` on batch directories, keep
                    query-batch rows, gzip, checkpoint; space auditor runs
                    throughout; mega tasks memory-gated (--max-mega)
merge_census.py     outputs/*.tsv.gz -> one table + species_summary.tsv
                    + census_qc.md (expected vs actual pairs, duplicates)
```

Task shape: a task is (species cluster, query batch of ~100 genomes). Its
worker builds a temp directory with the batch's TGTs plus the whole cluster's,
runs `syn2b synteny` once, and keeps only rows whose `genome_A` is in the
batch (intra-batch rows keep lexicographic order). This makes tasks exactly
partition the output rows, idempotent, and safely re-runnable.

## Commands

```bash
# 1. plan (local, any machine)
python3 scripts/gtdb_census/plan_census.py \
    --taxonomy data/gtdb_metadata/accession_taxonomy_r207.tsv.gz \
    --outdir results/census --batch-size 100 --jobs 6

# 2. prepare (HPC, where the genomes live)
python3 scripts/gtdb_census/prepare_workdir.py \
    --taxonomy ... --genome-dir /lustre/.../gtdb_r207 \
    --tasks results/census/tasks.jsonl --workdir $WORKDIR

# 3. run (a few long jobs; resubmit the same array to resume)
sbatch scripts/gtdb_census/census_long_job.slurm

# 4. merge + QC
python3 scripts/gtdb_census/merge_census.py \
    --workdir $WORKDIR --out results/census
```

## Space budget (measured)

| item | size |
|---|---:|
| TGT store, 317,542 genomes (~65 B/tag, ~350 KB/genome) | ~110 GB |
| task outputs, 711M rows x ~70 B, gzipped | ~12-15 GB |
| uncompressed task outputs before merge | ~50 GB |
| total workdir ceiling (hard limit) | 320 GB default |

Default thresholds: soft 250 GB, hard 320 GB, min free 50 GB — all flags
(`--soft-gb/--hard-gb/--min-free-gb`) to be tuned to the actual quota before
submission.

## QC gates

- self-comparison rows never emitted (`genome_A != genome_B` filter)
- duplicate pairs suppressed and counted in `census_qc.md`
- per-cluster expected C(n,2) vs actual rows (missing genomes surfaced)
- controls to run once on the HPC before the array: a genome vs a renamed copy
  of itself = 0 junctions; EDL933 vs Sakai = 2 (recorded in
  `results/census/controls.md`)
