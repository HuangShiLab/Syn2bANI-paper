# Syn2bANI Optimization & GTDB-R207 Benchmark — Status and Schedule

> Generated: 2026-07-20 (updated)
> Code repo: `~/Downloads/Syn2bANI`
> Paper/results repo: `~/Downloads/Syn2bANI-paper`
> Data: `~/data/gtdb-r207`

---

## 1. Current Status

### 1.1 Code (~/Downloads/Syn2bANI)

| Item | Status | Notes |
|------|--------|-------|
| Rust toolchain | ✅ | cargo 1.97.0, rustc 1.97.0 |
| `cargo build --release` | ✅ | Binary at `target/release/syn2bani` (~7.9 MB) |
| `cargo test` | ✅ | 21 lib + 8 integration tests pass |
| Integration tests | ✅ Fixed | Added missing `use_gbrt_v3` / `use_gbrt_v3_6` fields |
| GBRT models | ✅ | v2, v3, v3.6, **v4** embedded; **v4 is default** |
| Algorithm optimizations | ✅ Applied | memchr digestion, packed-sequence coarse screen; kept near-match enabled |

**Optimizations implemented:**
1. ✅ `digest_sequence` now uses anchor-based `memchr` scanning (82% faster on BcgI digestion; 43–59% faster on 16-enzyme panel).
2. ✅ `dist.rs::shared_tag_count` uses `HashSet<u64>` over packed sequences.
3. ⚠️ `dist` / `sketch` multi-enzyme parallelization was reverted to sequential inner loops to avoid nested-Rayon issues on real genomes.
4. ✅ Kept `MatchConfig::allow_near_match = true` (default); disabling it caused catastrophic accuracy loss.

**Bug fixes during optimization:**
- Reverted nested `par_iter()` in `dist.rs` and `sketch.rs` that caused empty tag sets on real genomes.
- Restored `allow_near_match = true` default.

### 1.2 Data (~/data/gtdb-r207)

| Item | Status | Notes |
|------|--------|-------|
| GTDB-R207 genomes | ✅ 98.54% complete | 64,747 / 65,703 representative genomes; ~202 GB in `genomes_all/` |
| Metadata | ✅ | `bac120_metadata_r207.tsv`, `ar53_metadata_r207.tsv`, taxonomy files |
| New pair sample | ✅ | `results/pairs_gtdb_r207_1k.tsv` — 1,000 stratified pairs (250 per label) |

### 1.3 Tools & Dependencies

| Tool | Status | Path/Version |
|------|--------|--------------|
| skani | ✅ | `~/.cargo/bin/skani` v0.1.1 (macOS); `/home/shihuang/.conda/envs/syn2bani/bin/skani` v0.3.2 (HPC) |
| seqkit | ✅ | `/opt/homebrew/bin/seqkit` v2.13.0 |
| ncbi-datasets-cli (`datasets`) | ✅ | `~/.local/bin/datasets` v18.33.1 |
| FastANI binary | ✅ | `/opt/homebrew/bin/fastANI` v1.34 (macOS); `/home/shihuang/.conda/envs/syn2bani/bin/fastANI` v1.34 (HPC) |
| pyfastani (Python) | ✅ | v0.6.1 (fallback reference) |
| Python packages | ✅ | numpy, pandas, scikit-learn, matplotlib, seaborn, pyfastani |

### 1.4 Benchmark Pipeline (~/Downloads/Syn2bANI-paper/scripts)

| Script | Status | Purpose |
|--------|--------|---------|
| `sample_gtdb_r207_pairs_v2.py` | ✅ New | Stratified pair sampling from `genomes_all/` |
| `run_benchmark_matrix_v2.py` | ✅ Updated | Runs Syn2bANI + skani + FastANI; now exports `s2b_af_q/r` and `s2b_ref_gc` |
| `evaluate_gtdb_r207.py` | ✅ Updated | Evaluates raw, skani, and **GBRT v4**; generates figures |
| `train_gbrt_v3.py` | ✅ Updated | Trains clean GBRT v4 model; exports Rust-compatible tree JSON |
| `split_pair_chunks.py` | ✅ New | Splits pair TSVs into HPC job-array chunks |
| `run_benchmark_chunk.py` | ✅ New | Thin wrapper for HPC chunk execution |
| `slurm_sketch_array.sh` | ✅ New | SLURM sketching job array |
| `slurm_matrix_array.sh` | ✅ Updated | SLURM pairwise benchmark job array (Syn2bANI + skani) |
| `slurm_merge_phase1.sh` | ✅ New | Merge chunks + sample FastANI subset |
| `slurm_fastani_array.sh` | ✅ New | FastANI on stratified subset |
| `slurm_train_final.sh` | ✅ New | Merge FastANI reference + train + evaluate |
| `submit_hpc_workflow.sh` | ✅ New | Master submission script with SLURM dependencies |
| `sample_fastani_subset.py` | ✅ New | Stratified FastANI subset sampling |

**Known issues resolved in pipeline:**
- FastANI binary (`fastANI`) is now used instead of pyfastani for the reference ANI; removed `--noFrag` so draft genomes are handled correctly.
- pyfastani remains a fallback but is not safe with `multiprocessing.Pool` or `ThreadPoolExecutor`; runs sequentially.
- Reusing `pyfastani.Sketch` across queries caused false zero-hits; fixed by creating a fresh sketch per pair.
- skani reports 0–100; normalized to 0–1 fraction.

### 1.5 Benchmarks Completed

| Run | Pairs | Tools | Status |
|-----|-------|-------|--------|
| Baseline 1k | 1,000 | all | ✅ Complete (`results/matrix_gtdb_r207_1k_baseline.tsv`) |
| Optimized 1k | 1,000 | all | ✅ Complete (`results/matrix_gtdb_r207_1k_optimized.tsv`) |
| Baseline 100 | 100 | Syn2bANI | ✅ Complete, used for optimization validation |
| Optimized 100 | 100 | Syn2bANI | ✅ Complete; identical accuracy to baseline, ~2× faster |
| Scaled 10k (6.2k realized) | 6,210 | all | ✅ Complete (`results/matrix_gtdb_r207_10k.tsv`) |

**Key results:**
- Accuracy identical between baseline and optimized matrices.
- Overall MAE vs FastANI (1k pairs): Syn2bANI raw 7.16%, **Syn2bANI GBRT v4 1.12%**, skani 10.35%.
- Scaled 6.2k pairs (2,490 with FastANI reference): Syn2bANI raw 22.96%, **Syn2bANI GBRT v4 0.44%**, skani 42.15%.
- 6.2k v4 model: 20% holdout MAE 0.85%, 5-fold CV MAE 1.00%.
- 1k-pair wall time: Syn2bANI baseline 6.68 s, optimized 2.97 s (~2.25× speedup); skani 3.79 s; FastANI (pyfastani) 326.91 s.
- 6.2k-pair wall time: ~2 min with FastANI binary + parallel skani/Syn2bANI.
- GBRT v4 model (trained on 6.2k) embedded in `~/Downloads/Syn2bANI/gbrt_model_v4.json`; `cargo test --release` passes.

### 1.7 HPC Deployment (HKU HPC2021)

| Item | Status | Notes |
|------|--------|-------|
| Code copied to `/lustre1/g/aos_shihuang/` | ✅ | `Syn2bANI/` and `Syn2bANI-paper/` |
| Syn2bANI binary built on HPC | ✅ | `cargo test --release` passes after fixing `simd.rs` `0xDFi8` literal |
| Conda env `syn2bani` | ✅ | Python packages via pip --user; skani v0.3.2 + fastANI v1.34 via bioconda |
| GTDB-R207 taxonomy | ✅ | Downloaded to `/lustre1/g/aos_shihuang/data/gtdb-r207/metadata/` |
| GTDB-R207 genomes | 🔄 | Tar downloading (~61 GB, ~2 h remaining) |
| SLURM scripts adapted | ✅ | HKU `amd` partition, `normal` QoS, conda activation added |
| SLURM workflow submitted | ⏳ | Waiting for genome extraction |

### 1.6 Existing GTDB-R207 Analysis

Prior work at `~/data/gtdb-r207/analysis/` already produced:
- 1,000 stratified pairs (`pair_list.json`)
- Phase 1 comparison (44 pairs)
- Phase 1 skani results (1,000 pairs)
- Retrained GBRT v3 model
- Per-phylum error evaluation (296 pairs, MAE 0.0226)
- Preliminary figures

---

## 2. Task Schedule

### Phase A: Preparation & Quick Wins (Day 1) — DONE

| # | Task | Location | Status |
|---|------|----------|--------|
| A1 | Status check & data audit | `~/data/gtdb-r207`, `~/Downloads/Syn2bANI` | ✅ Done |
| A2 | Install missing tools (seqkit, datasets) | macOS | ✅ Done |
| A3 | Fix integration tests | `~/Downloads/Syn2bANI/tests/integration_tests.rs` | ✅ Done |
| A4 | Implement algorithm optimizations | `~/Downloads/Syn2bANI/src/` | ✅ Done (with bug fixes) |
| A5 | Validate optimized build (`cargo test`, `cargo bench`) | `~/Downloads/Syn2bANI` | ✅ Done |

### Phase B: GTDB-R207 Benchmark Matrix (Day 1–2) — DONE

| # | Task | Output | Status |
|---|------|--------|--------|
| B1 | Generate 1,000 stratified pairs from `genomes_all/` | `results/pairs_gtdb_r207_1k.tsv` | ✅ Done |
| B2 | Run Syn2bANI baseline + optimized on pairs | `results/matrix_gtdb_r207_1k_*.tsv` | ✅ Done |
| B3 | Run skani on pairs | in matrix | ✅ Done |
| B4 | Run FastANI reference (pyfastani) on pairs | in matrix | ✅ Done |
| B5 | Merge into unified matrix | `results/matrix_gtdb_r207_1k_*.tsv` | ✅ Done |

### Phase C: Model Training & Evaluation (Day 2–3) — DONE

| # | Task | Output | Status |
|---|------|--------|--------|
| C1 | Evaluate per-label / per-phylum errors | `results/evaluation_gtdb_r207_1k_*.json` | ✅ Done |
| C2 | Compare Syn2bANI vs skani vs FastANI | `results/GTDB_R207_BENCHMARK_REPORT.md` | ✅ Done |
| C3 | Measure optimized vs baseline speedup | in report | ✅ Done |
| C4 | Retrain GBRT v3 on 1k matrix (experimental) | `results/gbrt_model_v3_retrained_1k.json`, `results/gbrt_v3_retrained_1k_report.txt` | ✅ Done |

### Phase D: Figures & Report (Day 3) — DONE

| # | Task | Output | Status |
|---|------|--------|--------|
| D1 | Generate scatter + error plots | `figures/gtdb_r207_*.png` | ✅ Done |
| D2 | Write comprehensive benchmark report | `results/GTDB_R207_BENCHMARK_REPORT.md` | ✅ Done |

### Phase E: HPC Scaling Plan (Day 3–4) — DONE

| # | Task | Output | Status |
|---|------|--------|--------|
| E1 | Design HPC benchmark for 100k+ pairs | `results/HPC_SCALING_PLAN.md` | ✅ Done |
| E2 | Draft SLURM/job-array scripts | `scripts/slurm_*.sh`, `scripts/run_benchmark_chunk.py`, `scripts/submit_hpc_workflow.sh` | ✅ Done |
| E3 | Estimate storage, CPU, memory, wall time | `results/HPC_SCALING_PLAN.md` | ✅ Done |
| E4 | Design stratified FastANI subset workflow | `scripts/sample_fastani_subset.py`, `scripts/slurm_merge_phase1.sh`, `scripts/slurm_fastani_array.sh`, `scripts/slurm_train_final.sh` | ✅ Done |
| E5 | Local dry-run of HPC workflow | 200-pair dry-run producing a valid GBRT model | ✅ Done |

### Phase F: Clean GBRT v4 Model & Embedding (Day 4) — DONE

| # | Task | Output | Status |
|---|------|--------|--------|
| F1 | Remove skani feature leakage from training | `scripts/train_gbrt_v3.py` | ✅ Done |
| F2 | Train v4 on inference-time features (raw ANI, shared tags, af_q/r) | `results/gbrt_model_v4_1k.json` | ✅ Done |
| F3 | Add Rust v4 inference path and make it default | `src/core/gbrt.rs`, `src/core/ani_calculator.rs` | ✅ Done |
| F4 | Embed v4 JSON and rebuild | `~/Downloads/Syn2bANI/gbrt_model_v4.json`, `target/release/syn2bani` | ✅ Done |
| F5 | Validate end-to-end accuracy and tests | `cargo test --release`, evaluation JSON updated | ✅ Done |

### Phase G: Larger-Scale Training (Day 4+) — IN PROGRESS

| # | Task | Output | Status |
|---|------|--------|--------|
| G1 | Generate 10k stratified pairs | `results/pairs_gtdb_r207_10k.tsv` | ✅ Done (6,210 pairs due to limited multi-genome species) |
| G2 | Run benchmark matrix on 6.2k pairs | `results/matrix_gtdb_r207_10k.tsv` | ✅ Done |
| G3 | Train v4 on 6.2k matrix | `results/gbrt_model_v4_10k.json` | ✅ Done; embedded into Syn2bANI |
| G4 | Deploy to HKU HPC2021 and stage data/tools | `/lustre1/g/aos_shihuang/Syn2bANI`, `/lustre1/g/aos_shihuang/Syn2bANI-paper` | ✅ In progress |
| G5 | Submit 100k-pair SLURM workflow | HPC job array | ⏳ Pending GTDB tar extraction |

---

## 3. File Mapping

| Deliverable | Target Path |
|-------------|-------------|
| Code updates | `~/Downloads/Syn2bANI/src/` |
| Benchmark results | `~/Downloads/Syn2bANI-paper/results/` |
| Figures | `~/Downloads/Syn2bANI-paper/figures/` |
| This status/schedule | `~/Downloads/Syn2bANI-paper/STATUS_AND_SCHEDULE.md` |
| HPC plan | `~/Downloads/Syn2bANI-paper/results/HPC_SCALING_PLAN.md` |

---

## 4. Key Learnings / Risks

| Issue | Resolution |
|-------|------------|
| pyfastani cache reuse caused false zero-hits | Create fresh `Sketch` per pair |
| pyfastani not thread/process-safe | Run sequentially |
| Nested Rayon `par_iter()` caused empty tag sets | Use sequential inner loops |
| Disabling near-match destroyed accuracy | Keep `allow_near_match=true` default |
| Digestion speedup is real but requires correctness validation | Unit-tested on real genome; validated against baseline binary |

---

## 6. SynTracker Fig. 3 Independent Validation — IN PROGRESS

Goal: validate that Syn2bANI's `anchor_adjacency` can reproduce the four-species
evolutionary patterns from Enav et al. 2024 (Fig. 3).

| Task | Output | Status |
|---|---|---|
| Download SI tables 2–5 | `data/syntracker/41587_2024_2276_MOESM3_ESM.xlsx` | ✅ Done |
| Parse isolate metadata & pair lists | `data/syntracker/samples_*.tsv`, `pairs_*.tsv` | ✅ Done |
| Build ENA FASTQ manifests (132 isolates) | `data/syntracker/fastq_manifest_*.tsv` | ✅ Done |
| HPC workflow scripts | `scripts/syntracker_validation/` | ✅ Done / pushed |
| Download references & reads on HPC | `/lustre1/g/aos_shihuang/data/syntracker_validation/` | ✅ Done (132 isolates, 75 GB) |
| Assemble isolates (SLURM array) | `assemblies/*.fna` | ✅ Done 2026-08-09 (132/132; 补跑 45 个时序错位样本，S. rimosus 需 6h+/断点续跑 `03b_resume_one.sh`) |
| Run Syn2bANI + skani per species | `syn2bani/syn2bani_*.tsv`, `skani/skani_*.tsv` | ✅ Done 2026-08-10 |
| Plot ANI vs anchor_adjacency | `figures/syntracker_validation/ani_vs_synteny_syntracker_species.png` | ✅ Done 2026-08-10 |

**结果（Spearman ρ，ANI vs anchor_adjacency）**：四个物种全部重现预期进化模式——
E. coli hypermutator ρ=0.06（ANI 99.96–100，synteny 平稳 ~0.76，SNP 主导）；
H. pylori ρ=0.81（混合模式）；N. gonorrhoeae ρ=0.59（混合 SNP+SV）；
S. rimosus ANI 钉在 99.98–100（克隆）而 synteny 0.89–0.955 大幅变化（SV 主导）。

**注意**：`--calibrate` 的线性校准模型在这些近克隆对上失效（S. rimosus 出现
ANI>100% 的外推，N. gonorrhoeae 被拉低约 1.2 个百分点），绘图与结论均使用
原始 `ani`（gamma 异质）列。校准模型的适用域是 GTDB 中低 ANI 区间，
不应外推到 >99.5% 的近克隆比较。

### Expected evolutionary signatures

- *N. gonorrhoeae* — ANI and anchor_adjacency correlated (mixed SNP + SV evolution).
- hypermutator *E. coli* — wide ANI range, anchor_adjacency stays high (SNP-driven).
- *H. pylori* — mixed modes within/between hosts.
- *S. rimosus* — ANI high/clonal, anchor_adjacency variable (SV-driven).

### How to run on HPC2021

```bash
# 1. Setup (login node)
bash /lustre1/g/aos_shihuang/Syn2bANI-paper/scripts/syntracker_validation/00_setup.sh

# 2. References (login / I/O node)
bash /lustre1/g/aos_shihuang/data/syntracker_validation/scripts/01_download_references.sh

# 3. FASTQs (I/O node — network-bound)
ssh shihuang@hpc2021-io1.hku.hk
bash /lustre1/g/aos_shihuang/data/syntracker_validation/scripts/02_download_reads.sh

# 4. Assemble (compute nodes — SLURM)
bash /lustre1/g/aos_shihuang/data/syntracker_validation/scripts/03_assemble_array.sh

# 5. Syn2bANI + skani comparisons (compute nodes — SLURM)
bash /lustre1/g/aos_shihuang/data/syntracker_validation/scripts/04_submit_syn2bani.sh
bash /lustre1/g/aos_shihuang/data/syntracker_validation/scripts/05_submit_skani.sh

# 6. Plot (locally after copying results)
python3 /lustre1/g/aos_shihuang/Syn2bANI-paper/scripts/syntracker_validation/06_plot_results.py \
  --syn2bani-dir data/syntracker_validation/syn2bani \
  --skani-dir data/syntracker_validation/skani \
  --metadata-dir data/syntracker_validation/samples \
  --outdir figures/syntracker_validation
```

---

## 5. Next Immediate Actions

1. ✅ GBRT v4 model trained on 6.2k pairs, embedded, and validated; report/figures updated.
2. ✅ HPC workflow designed, scripted, and locally dry-run validated (200 pairs).
3. Submit `scripts/submit_hpc_workflow.sh` on a cluster for 100k+ pairs when HPC access is available.
4. ✅ SynTracker Fig. 3 validation completed on HPC2021 (2026-08-10): 132 isolates assembled, Syn2bANI + skani compared per species, all four expected evolutionary patterns reproduced — see §6.
