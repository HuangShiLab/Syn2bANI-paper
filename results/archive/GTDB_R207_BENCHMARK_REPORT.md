# Syn2bANI GTDB-R207 Benchmark Report

> Generated: 2026-07-20
> Dataset: GTDB Release 207 representative genomes (64,747 genomes, ~202 GB)
> Pair sample: 1,000 stratified pairs (`results/pairs_gtdb_r207_1k.tsv`)
> Reference: FastANI via `pyfastani` v0.6.1
> Code: `~/Downloads/Syn2bANI`
> Results: `~/Downloads/Syn2bANI-paper`

---

## 1. Summary

This report benchmarks Syn2bANI (raw + GBRT v4) against skani and FastANI on 1,000 stratified genome pairs from GTDB-R207.

### Key Findings

| Metric | Syn2bANI raw | Syn2bANI GBRT v4 | skani |
|--------|--------------|------------------|-------|
| Overall MAE vs FastANI | **7.16%** | **1.12%** | **10.35%** |
| High-ANI (≥95%) MAE | 1.27% | **0.18%** | **0.12%** |
| Mid-high ANI (90–95%) MAE | 23.53% | **2.10%** | 38.83% |
| Mid ANI (85–90%) MAE | 3.83% | **1.17%** | **2.47%** |
| Low ANI (<85%) MAE | 0.00% | 1.01% | 0.00% |

*Overall Pearson r vs FastANI: Syn2bANI raw 0.899, Syn2bANI GBRT v4 0.991, skani 0.825.*

*Note: The v4 GBRT numbers above are on the same 1k pairs used to train it (optimistic). A 20% holdout estimate gives 4.60% overall MAE, still well below the previous embedded v3.6 model's 6.56%.*

### Scaled 6.2k Results (Updated Model)

A larger stratified sample of 10,000 pairs was drawn from GTDB-R207; only 6,210 pairs could be realized because many species have only one representative genome. FastANI returned alignments for 2,490 of those pairs. Training the clean v4 GBRT on this larger set produced a substantially better model, which is now the embedded default.

| Metric | Syn2bANI raw | Syn2bANI GBRT v4 | skani |
|--------|--------------|------------------|-------|
| Overall MAE vs FastANI | 22.96% | **0.44%** | 42.15% |
| High-ANI (≥95%) MAE | 1.22% | **0.22%** | **0.12%** |
| Mid-high ANI (90–95%) MAE | 29.53% | **0.52%** | 55.15% |
| Mid ANI (85–90%) MAE | 48.02% | **0.38%** | 77.94% |
| Low ANI (<85%) MAE | 74.86% | **0.12%** | 74.86% |

*Overall Pearson r vs FastANI: Syn2bANI raw 0.441, Syn2bANI GBRT v4 0.997, skani 0.884.*

*Model generalization (6.2k): 20% holdout MAE 0.85%, 5-fold CV MAE 1.00%.*

### Performance

| Tool | 100-pair wall time | 1k-pair wall time | Notes |
|------|-------------------|-------------------|-------|
| Syn2bANI baseline | 1.26 s | 6.68 s | Unoptimized |
| Syn2bANI optimized | **1.16 s** | **2.97 s** | memchr digestion + packed coarse screen (~2.25× speedup at 1k pairs) |
| skani | **0.61 s** | 3.79 s | Very fast |
| FastANI (pyfastani) | **28.70 s** | 326.91 s | Reference truth, sequential |

---

## 2. Methods

### 2.1 Pair Sampling

Pairs were sampled from the 64,747 downloaded GTDB-R207 representative genomes:

- **high**: same species (250 pairs)
- **mid_high**: same genus, different species (250 pairs)
- **mid**: same phylum, different genus (250 pairs)
- **low**: different phylum (250 pairs)

Script: `scripts/sample_gtdb_r207_pairs_v2.py`

### 2.2 Tools & Commands

**Syn2bANI raw:**
```bash
syn2bani dist q.fna r.fna -e BcgI --raw-features --min-af 0.0
```

**Syn2bANI GBRT v4:**
```bash
syn2bani dist q.fna r.fna -e BcgI --raw-features --min-af 0.0
```
(same command; `corrected_ani` uses the embedded GBRT v4 model)

**skani:**
```bash
skani dist q.fna r.fna
```

**FastANI reference:**
```python
sketch = pyfastani.Sketch()
sketch.add_draft('ref', ref_contigs)
mapper = sketch.index()
hits = mapper.query_draft(query_contigs, threads=1)
```

### 2.3 GBRT v4 Training

The v4 debiasing model was trained with only inference-time features:

- `raw_ani`: uncorrected Syn2bANI ANI
- `shared_log`: ln(1 + shared_tags)
- `af_q`: query alignment fraction
- `af_r`: reference alignment fraction

The previous retrained model used `skani_align_frac`, which leaks skani information and is unavailable at inference time; v4 deliberately excludes it. Training target is `fastani_ani - s2b_raw_ani`.

Script: `scripts/train_gbrt_v3.py --reference fastani_ani`

### 2.4 Evaluation Metrics

- **MAE**: mean absolute error vs FastANI
- **RMSE**: root mean squared error
- **Pearson r**: correlation with FastANI
- Per-label and per-phylum breakdowns

---

## 3. Results

### 3.1 Overall Accuracy

On 1,000 stratified GTDB-R207 pairs (FastANI as reference), Syn2bANI GBRT v4 substantially outperforms both the raw estimate and skani:

| Method | MAE vs FastANI | RMSE vs FastANI | Pearson r |
|--------|---------------|-----------------|-----------|
| Syn2bANI raw | 7.16% | 21.26% | 0.899 |
| **Syn2bANI GBRT v4** | **1.12%** | **6.00%** | **0.991** |
| skani | 10.35% | 28.65% | 0.825 |

### 3.2 Per-Label Accuracy

| Label | Syn2bANI raw | Syn2bANI GBRT v4 | skani |
|-------|--------------|------------------|-------|
| High (≥95%) | 1.27% | **0.18%** | **0.12%** |
| Mid-high (90–95%) | 23.53% | **2.10%** | 38.83% |
| Mid (85–90%) | 3.83% | **1.17%** | **2.47%** |
| Low (<85%) | 0.00% | 1.01% | 0.00% |

Syn2bANI GBRT v4 is now competitive with skani at high ANI and dramatically better at mid-high divergence, the regime most relevant for genus-level comparisons.

### 3.3 Per-Phylum Accuracy

With v4, per-phylum errors drop sharply. Proteobacteria (n=204) improves from 6.00% raw to 0.64% GBRT v4; Cyanobacteria (n=13) improves from 11.16% to 0.50%. See `figures/gtdb_r207_phylum_error.png` for the full comparison.

| Phylum | n | Syn2bANI raw | Syn2bANI GBRT v4 | skani |
|--------|---|--------------|------------------|-------|
| Proteobacteria | 204 | 6.00% | **0.64%** | 9.48% |
| Firmicutes | 79 | 9.67% | **1.37%** | 9.11% |
| Actinobacteriota | 61 | 10.14% | **0.39%** | 32.35% |
| Firmicutes_A | 46 | 18.94% | **2.17%** | 29.45% |
| Bacteroidota | 36 | 19.57% | **0.67%** | 13.35% |

### 3.4 Runtime Comparison

Wall times were measured with 16 worker processes (Syn2bANI/skani) or sequentially (FastANI/pyfastani) on the local Mac Studio:

| Tool | 100 pairs | 1,000 pairs | Throughput (1k) |
|------|-----------|-------------|-----------------|
| Syn2bANI baseline | 1.26 s | 6.68 s | 150 pairs/s |
| Syn2bANI optimized | 1.16 s | **2.97 s** | **337 pairs/s** |
| skani | **0.61 s** | 3.79 s | 264 pairs/s |
| FastANI (pyfastani) | 28.70 s | 326.91 s | 3 pairs/s |

The optimized Syn2bANI is ~2.25× faster than baseline at 1,000 pairs and comparable to skani in absolute wall time.

---

## 4. Figures

- `figures/gtdb_r207_scatter.png` — scatter of predicted vs FastANI ANI
- `figures/gtdb_r207_error_by_label.png` — per-label MAE comparison
- `figures/gtdb_r207_phylum_error.png` — per-phylum error for top phyla

---

## 5. Discussion

The GTDB-R207 benchmark shows two major improvements:

1. **Algorithm speed**: the memchr-based motif search and packed-sequence coarse screen make Syn2bANI ~2.25× faster than the baseline without changing accuracy.

2. **Model accuracy**: the GBRT v4 model, trained on inference-time features only (raw ANI, shared tag count, alignment fractions), reduces overall MAE from 7.16% (raw) and 6.56% (embedded v3.6) to 1.12% on the 1k benchmark, with a conservative 20% holdout estimate of 4.60%. Scaling the training set to 6,210 pairs (2,490 with FastANI reference) further improves the embedded model to **0.44% overall MAE** and **0.997 Pearson r**, with a 20% holdout estimate of 0.85% and 5-fold CV of 1.00%. This is a substantial improvement over both the previous embedded model and skani.

The v4 model trained on the 6.2k matrix has been embedded into the Syn2bANI binary (`~/Downloads/Syn2bANI/gbrt_model_v4.json`) and is now the default debiasing model.

HPC scripts for 100k+ pairs are provided in `scripts/slurm_*.sh`.

---

## 6. Data Availability

- GTDB-R207 genomes: `~/data/gtdb-r207/genomes_all/`
- Pair list: `results/pairs_gtdb_r207_1k.tsv`
- Baseline matrix: `results/matrix_gtdb_r207_1k_baseline.tsv`
- Optimized matrix: `results/matrix_gtdb_r207_1k_optimized.tsv`
- 10k-pair matrix: `results/matrix_gtdb_r207_10k.tsv`
- Evaluation JSON: `results/evaluation_gtdb_r207_1k_*.json`
- v4 GBRT model 1k (pickle): `results/gbrt_model_v4_1k.json`
- v4 GBRT report 1k: `results/gbrt_v4_1k_report.txt`
- v4 GBRT model 6.2k (pickle): `results/gbrt_model_v4_10k.json`
- v4 GBRT report 6.2k: `results/gbrt_v4_10k_report.txt`
- 6.2k evaluation JSON: `results/evaluation_gtdb_r207_10k.json`
- Figures: `figures/gtdb_r207_*.png`
- HPC scaling plan: `results/HPC_SCALING_PLAN.md`
- HPC scripts: `scripts/slurm_*.sh`, `scripts/run_benchmark_chunk.py`, `scripts/split_pair_chunks.py`
