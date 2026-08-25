# Syn2bANI HPC Scaling Plan

> Goal: Benchmark Syn2bANI against skani/FastANI on datasets larger than the GTDB-R207 1k-pair sample (e.g., 10k–100k pairs, 1,000+ genomes).

---

## 1. Target Datasets

| Dataset | Size | Pairs | Purpose |
|---------|------|-------|---------|
| GTDB-R207 full representatives | 65k genomes, ~200 GB | All-vs-all ≈ 2.1B pairs (infeasible) | Use sketch-based pre-filter or sample |
| GTDB-R207 10k stratified pairs | Subset of 65k | 10,000 | Direct accuracy comparison |
| GTDB-R207 100k pairs | Subset of 65k | 100,000 | Throughput / scaling test |
| RefSeq complete bacteria | ~30k genomes, ~50 GB | Sampled 10k–100k pairs | Cross-database validation |
| EBI MGnify/UHGG MAGs | 4k–100k MAGs | Sampled pairs | Fragmented-genome stress test |

---

## 2. Recommended HPC Resources

### For 10k–100k pairwise comparisons

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU cores | 32 | 128–256 |
| Memory | 64 GB | 256 GB |
| Storage (scratch) | 500 GB | 2 TB |
| Wall time (10k pairs) | 4 h | 1 h |
| Wall time (100k pairs) | 40 h | 8 h |

### Bottlenecks

- **FastANI (pyfastani)**: single-pair indexing is CPU-bound; parallelize by pair.
- **skani**: very fast per pair, but 100k pairs still need many cores.
- **Syn2bANI**: digestion + matching per pair; most time spent in digestion. The memchr optimization is critical here.
- **I/O**: reading 200 GB of FASTA files repeatedly is expensive. Strategy:
  - Convert all genomes to `.s2ba` sketches once (Syn2bANI `sketch` / `db build`).
  - Run `syn2bANI db search` or in-memory pairwise comparison from sketches.

---

## 3. Job-Array Workflow (SLURM)

The workflow is split into four dependent steps. FastANI is only run on a stratified subset; Syn2bANI + skani are run on all pairs.

### Phase 1: Pre-compute sketches (one-time, optional)

`scripts/slurm_sketch_array.sh` pre-computes Syn2bANI sketches for every genome. This is optional but speeds up repeated pairwise runs.

```bash
#!/bin/bash -l
#SBATCH --job-name=s2b_sketch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-1000%50
```

### Phase 2: Pairwise Syn2bANI + skani on all pairs

`scripts/slurm_matrix_array.sh` runs Syn2bANI and skani on pair chunks.

```bash
#!/bin/bash -l
#SBATCH --job-name=s2b_matrix
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --array=1-100%10
```

### Phase 3: Merge + sample FastANI subset

`scripts/slurm_merge_phase1.sh` merges the chunk TSVs, then samples a stratified subset for FastANI reference.

### Phase 4: FastANI on the stratified subset

`scripts/slurm_fastani_array.sh` runs FastANI only on the sampled pairs.

```bash
#!/bin/bash -l
#SBATCH --job-name=s2b_fastani
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=1-100%20
```

### Phase 5: Final merge, train, and evaluate

`scripts/slurm_train_final.sh` merges the FastANI results into the full matrix, trains the GBRT v4 model, and evaluates.

### Dependency chain

```bash
MATRIX=$(sbatch scripts/slurm_matrix_array.sh | awk '{print $4}')
MERGE=$(sbatch --dependency=afterok:$MATRIX scripts/slurm_merge_phase1.sh | awk '{print $4}')
FASTANI=$(sbatch --dependency=afterok:$MERGE scripts/slurm_fastani_array.sh | awk '{print $4}')
sbatch --dependency=afterok:$FASTANI scripts/slurm_train_final.sh
```

---

## 4. Checkpoint & Resume Strategy

The `run_benchmark_matrix_v2.py` script already supports resume: if `--output` exists, completed `(query, reference)` pairs are skipped. This allows:

- Restarting after wall-time limits.
- Running partial chunks on preemptible nodes.
- Merging chunks into a final matrix.

For HPC, split the pair list into chunks of 1,000–5,000 pairs and run each as a separate job-array task.

---

## 5. Performance Monitoring

For each job, capture:

```bash
/usr/bin/time -v python3 run_benchmark_matrix_v2.py ... > job.log 2>&1
```

Key metrics to record:
- Elapsed wall time
- Maximum resident set size (peak RAM)
- CPU utilization
- Pairs completed per hour
- Failures / timeouts per tool

---

## 6. Comparison Design

| Comparison | Command | Notes |
|------------|---------|-------|
| Syn2bANI raw | `syn2bani dist q.fna r.fna -e BcgI --raw-features --min-af 0.0` | Baseline without GBRT |
| Syn2bANI GBRT v4 | `syn2bani dist q.fna r.fna -e BcgI --raw-features --min-af 0.0` | With embedded v4 clean model |
| skani | `skani dist q.fna r.fna` | Fast ANI approximation |
| FastANI | `pyfastani` or `fastANI` binary | Reference truth |

For 100k+ pairs, run Syn2bANI + skani on every pair and compute FastANI only on a stratified subset (e.g., 10–20k pairs) as reference, because FastANI is the slowest step. Stratification uses `scripts/sample_fastani_subset.py` on the merged skani matrix.

### Model training on HPC data

After merging the FastANI reference into the full matrix, train the production GBRT v4 with only inference-time features:

```bash
python3 scripts/train_gbrt_v3.py \
  --matrix results/matrix_gtdb_r207_100k.tsv \
  --output results/gbrt_model_v4_100k.json \
  --report results/gbrt_v4_100k_report.txt \
  --reference fastani_ani
```

Copy the resulting JSON to `~/Downloads/Syn2bANI/gbrt_model_v4.json` and run `cargo build --release` to embed it.

---

## 7. Expected Outputs

| File | Description |
|------|-------------|
| File | Description |
|------|-------------|
| `results/matrix_gtdb_r207_100k_skani.tsv` | Unified Syn2bANI + skani matrix |
| `results/pairs_fastani_subset_100k.tsv` | Stratified FastANI subset pairs |
| `results/matrix_gtdb_r207_100k.tsv` | Final matrix with FastANI reference merged |
| `results/evaluation_gtdb_r207_100k.json` | MAE/RMSE/Pearson by label and phylum |
| `results/gbrt_model_v4_100k.json` | Production GBRT v4 model |
| `figures/gtdb_r207_scatter.png` | Scatter plots |
| `figures/gtdb_r207_error_by_label.png` | Error by pair type |
| `figures/gtdb_r207_phylum_error.png` | Per-phylum MAE |
| `results/HPC_RUNTIME_REPORT.md` | Wall time, memory, throughput |

---

## 8. Local Mac Studio vs HPC

| Task | Mac Studio | HPC |
|------|-----------|-----|
| Algorithm dev & unit tests | ✅ Ideal | Possible but less convenient |
| 1k–5k pair benchmarks | ✅ Fast | Overhead not worth it |
| 10k–100k pair benchmarks | ⚠️ Feasible but slow | ✅ Recommended |
| Full GTDB all-vs-all | ❌ Infeasible | ✅ Required |

---

## 9. HKU HPC2021 Deployment Notes

Actual deployment on `hpc2021.hku.hk` uses the `amd` partition and `/lustre1/g/aos_shihuang/` storage.

### Setup on HPC2021

1. **Login node setup** (compute nodes have no internet):
   ```bash
   ssh shihuang@hpc2021.hku.hk
   cd /lustre1/g/aos_shihuang/Syn2bANI-paper
   # Install conda env + skani + fastANI
   bash scripts/hpc_setup_env.sh
   # Build Syn2bANI
   bash scripts/hpc_build.sh
   # Download GTDB-R207 genomes (login node only)
   bash scripts/hpc_download_gtdb_r207.sh
   ```

2. **Extract and sample** (can be submitted as SLURM once the tar is on disk):
   ```bash
   sbatch scripts/hpc_extract_and_sample.sh
   ```

3. **Submit the full workflow**:
   ```bash
   bash scripts/submit_hpc_workflow.sh
   ```

### Important HPC2021 constraints

- Compute nodes **do not have internet access**. All downloads and `cargo build` must run on the login node.
- The conda env is created at `/home/shihuang/.conda/envs/syn2bani`; SLURM scripts activate it explicitly.
- Wall time for `amd` partition jobs is set conservatively (4–8 h for array tasks, 4 h for merge/train).

## 10. Next Steps for HPC Execution

1. Confirm access to a cluster with SLURM.
2. Transfer the optimized `syn2bani` binary, `skani`, and `fastANI` binaries to `/scratch`.
3. Stage GTDB-R207 genomes on scratch storage.
4. Create pair chunks with `scripts/split_pair_chunks.py`.
5. Submit the dependency chain:
   ```bash
   MATRIX=$(sbatch scripts/slurm_matrix_array.sh | awk '{print $4}')
   MERGE=$(sbatch --dependency=afterok:$MATRIX scripts/slurm_merge_phase1.sh | awk '{print $4}')
   FASTANI=$(sbatch --dependency=afterok:$MERGE scripts/slurm_fastani_array.sh | awk '{print $4}')
   sbatch --dependency=afterok:$FASTANI scripts/slurm_train_final.sh
   ```
6. Compare throughput (pairs/hour) and memory vs skani/FastANI.
7. Copy the best `gbrt_model_v4_*.json` to `~/Downloads/Syn2bANI/gbrt_model_v4.json` and rebuild.

---

*Plan version: 2026-07-20*
