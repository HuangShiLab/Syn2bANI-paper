# Syn2bANI-paper

> Paper repository for **"Syn2b-ANI: Strain-level ANI estimation via fixed restriction-site anchors for fragmented metagenome-assembled genomes"**

This repository contains the manuscript, analysis code, benchmark data, and figures for the Syn2bANI paper.

---

## Repository Structure

```
Syn2bANI-paper/
├── paper/              # Manuscript drafts, outlines, and methods documentation
│   ├── MANUSCRIPT_OUTLINE.md
│   ├── GBRT_DEBIAS_EXPLANATION.md
│   ├── GBRT_METHODS_FOR_MANUSCRIPT.md
│   └── GBRT_V2_TRAINING_REPORT.md
├── analysis/           # Python scripts for analysis and figure generation
│   ├── benchmark_pipeline.py
│   ├── performance_benchmark.py
│   ├── plot_comparison.py
│   ├── plot_final.py
│   ├── plot_multispecies.py
│   ├── plot_performance.py
│   ├── train_gbrt_v2.py
│   ├── validate_gbrt_v2.py
│   └── ...
├── data/               # Benchmark data and intermediate results
│   ├── benchmarks/     # CSV tables: accuracy, performance, multi-species, SV validation
│   ├── genomes/        # Key reference genomes (E. coli K-12, simulated MAGs)
│   └── models/         # GBRT model files (JSON decision trees, Python pickles)
├── figures/            # All publication-ready figures
│   ├── benchmark_*.png         # Enzyme extraction performance
│   ├── debiasing_*.png         # GBRT debiasing validation
│   ├── performance_*.png       # Runtime vs skani/FastANI
│   ├── multispecies_validation.png
│   ├── realistic_mag_results.png
│   └── skani_figure2_style.png
└── results/            # Detailed benchmark reports and summaries
    ├── BENCHMARK_REPORT.md
    ├── COMPREHENSIVE_BENCHMARK.md
    ├── FINAL_BENCHMARK_REPORT.md
    ├── FINAL_REPORT_v0.1.1.md
    └── HEAD_TO_HEAD_REPORT.md
```

---

## Key Figures

| Figure | Description | Script |
|--------|-------------|--------|
| `benchmark_16enzymes.png` | 16-enzyme panel extraction speed | `benches/benchmark.rs` |
| `benchmark_new_vs_old.png` | Fast2bRAD-M vs legacy margin-based | `benches/benchmark.rs` |
| `debiasing_comparison.png` | GBRT vs simple debias across divergence | `plot_final.py` |
| `performance_comparison_corrected.png` | Runtime vs skani/FastANI (log scale) | `performance_benchmark.py` |
| `multispecies_validation.png` | Cross-species GBRT generalization | `plot_multispecies.py` |
| `realistic_mag_results.png` | Chimera, contamination, duplication robustness | `task3_realistic_mag.py` |

---

## Data Files

### Benchmark CSVs

| File | Description | Rows |
|------|-------------|------|
| `comparison_results.csv` | Accuracy vs FastANI (completeness + N50) | ~20 |
| `performance_results.csv` | Runtime benchmark (Syn2bANI vs skani) | ~50 |
| `multispecies_results.csv` | Cross-species validation (5 species) | 5 |
| `sv_validation_results.csv` | Structural variation detection accuracy | 4 |
| `multienzyme_consensus.csv` | Multi-enzyme ANI consensus | 10 |
| `training_data_v2.csv` | GBRT v2 training data (49 species) | 1,260 |
| `genome_metadata.json` | 49-species metadata for GBRT training | 49 entries |

### Models

| File | Description | Size |
|------|-------------|------|
| `gbrt_model_v2.json` | Embedded GBRT model (300 trees, depth 5) | 1.08 MB |
| `gbrt_model_runtime.json` | Runtime-optimized GBRT (200 trees, depth 4) | 0.64 MB |
| `syn2bani_gbrt_debias_model.pkl` | Python sklearn model (for retraining) | 0.36 MB |

---

## Reproducing Figures

Most figures can be regenerated from the provided CSV data using the scripts in `analysis/`:

```bash
cd analysis
python3 plot_final.py          # Figure 1: Accuracy comparison
python3 plot_performance.py    # Figure 2: Performance benchmark
python3 plot_multispecies.py   # Figure 3: Cross-species validation
python3 task3_realistic_mag.py # Figure 4: Realistic MAG robustness
python3 sv_simulation.py       # Figure 5: SV detection validation
```

> **Note**: The Rust benchmark in `benches/benchmark.rs` (from the [Syn2bANI](https://github.com/HuangShiLab/Syn2bANI) code repo) generates the enzyme extraction speed figures. Run `cargo bench` in the code repo to regenerate.

---

## Manuscript Status

| Section | Status | Notes |
|---------|--------|-------|
| Abstract | 📝 Draft | Needs final polish |
| Introduction | 📝 Draft | 2bRAD-M + skani context |
| Methods | ✅ Complete | All algorithms documented |
| Results — Accuracy | ✅ Complete | Benchmarked vs FastANI |
| Results — Performance | ✅ Complete | vs skani & FastANI |
| Results — Fragmentation | ✅ Complete | N50 500 bp–100 kb |
| Results — SV Detection | ✅ Complete | 4 SV types validated |
| Results — Multi-species | ✅ Complete | 5 species, GBRT v2 |
| Discussion | 📝 Draft | Needs final synthesis |
| Figures | ✅ Complete | 13 figures generated |

---

## Citation

```bibtex
@article{syn2bani2025,
  title={Syn2b-ANI: Strain-level ANI estimation via fixed restriction-site anchors for fragmented metagenome-assembled genomes},
  author={HuangShiLab},
  journal={In preparation},
  year={2025}
}
```

---

## Related Repositories

- **[Syn2bANI](https://github.com/HuangShiLab/Syn2bANI)** — Main code repository (Rust)
- **[Syn2b](https://github.com/HuangShiLab/Syn2b)** — Synteny analysis tool (upstream)
- **[Fast2bRAD-M](https://github.com/HuangShiLab/Fast2bRAD-M)** — Fast tag extraction (methodology)

---

## License

Analysis scripts and data are released under the MIT License. The manuscript text is © 2025 HuangShiLab.
# Syn2bANI-paper
