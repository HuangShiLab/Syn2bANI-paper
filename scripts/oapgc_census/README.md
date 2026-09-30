# OAPGC within-species census

OAPGC publishes 149,921 HQMAGs, but its annotation marks **99,215
strain-level non-redundant representatives** with a non-null
`Represent(STRAINs)` field. Grouping those representatives by `SGBs` gives
1,488 multi-genome clusters and 30,560,247 unique within-species pairs.

## Local feasibility

The measured/planning estimate on 16 Apple Silicon cores is:

- structural tasks: 2,120 at batch size 250;
- Syn2b cost: about 406 core-h;
- elapsed estimate: about 25 h (excluding ANI and I/O overhead);
- gzipped structural output: about 0.9 GB;
- extracted representative FASTAs: likely 60–100 GB;
- local free space at planning time: about 1.2 TB.

A local run is therefore feasible. The main uncertainty is whether MAG assembly
fragmentation changes the 9 ms/ordered-pair planning constant; a largest-20
cluster smoke test should be run first.

A partial-tar smoke test using 100 representative FASTAs successfully generated
all expected 58 pairs and 67 checkpoints on the local Syn2b binary. This
confirms gzip extraction, genome-ID handling, and the census worker workflow;
it is not a throughput benchmark.

## Workflow

```bash
# 1. Plan from annotation (already run)
python3 scripts/oapgc_census/prepare_oapgc_census.py \
  --annotation /Users/macstudio/Downloads/OAPGC/raw/cluster.annotation.xlsx \
  --workdir /Users/macstudio/Downloads/OAPGC/census_plan \
  --batch-size 250 --node-cores 16 --wall-hours 48

# 2. Extract only representative FASTAs from the five genome tarballs
python3 scripts/oapgc_census/extract_oapgc_representatives.py \
  --archives /Users/macstudio/Downloads/OAPGC/raw/OAPGC_genomes_part{1..5}.tar.gz \
  --representatives /Users/macstudio/Downloads/OAPGC/census_plan/representative_genomes.txt \
  --out /Users/macstudio/Downloads/OAPGC/genomes

# 3. Re-run preparation with genome-root to create manifest.json
python3 scripts/oapgc_census/prepare_oapgc_census.py \
  --annotation /Users/macstudio/Downloads/OAPGC/raw/cluster.annotation.xlsx \
  --genome-root /Users/macstudio/Downloads/OAPGC/genomes \
  --workdir /Users/macstudio/Downloads/OAPGC/census \
  --batch-size 250 --node-cores 16 --wall-hours 48

# 4. Smoke test: first 20 largest clusters
# Use the largest-20 output/checkpoint directories, then merge and audit before full run.

# 5. Run all tasks locally
python3 scripts/gtdb_census/census_worker.py \
  --workdir /Users/macstudio/Downloads/OAPGC/census \
  --tasks /Users/macstudio/Downloads/OAPGC/census/tasks.jsonl \
  --syn2b /Users/macstudio/Downloads/Syn2b/target/release/syn2b \
  --ensure-filename-genome-id --workers 12 --max-mega 4 \
  --soft-gb 250 --hard-gb 350 --min-free-gb 50 --drain

# 6. Merge structural output and QC
python3 scripts/gtdb_census/merge_census.py \
  --workdir /Users/macstudio/Downloads/OAPGC/census \
  --out results/census/oapgc
```

Keep raw FASTAs and tarballs outside Git.
