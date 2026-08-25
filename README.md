# Syn2bANI-paper

> Paper repository for **"Syn2b-ANI: Strain-level ANI estimation and structural comparison via fixed restriction-site anchors"**

This repository contains the manuscript, analysis code, benchmark data, and figures for the Syn2bANI paper.

---

## Repository Structure

```
Syn2bANI-paper/
├── paper/              # Manuscript, supplementary, and project indices
│   ├── manuscript/     # Current manuscript and supplementary files
│   │   ├── manuscript.md
│   │   ├── manuscript.docx
│   │   ├── SUPPLEMENTARY.md
│   │   └── SUPPLEMENTARY.docx
│   ├── figures/        # Figure index (source files in figures/)
│   │   └── FIGURES.md
│   ├── data/           # Data index and simulation/performance report
│   │   ├── DATA.md
│   │   └── SIMULATION_AND_PERFORMANCE_REPORT.md
│   ├── scripts/        # Script index (source files in scripts/ and analysis/)
│   │   └── SCRIPTS.md
│   └── others/         # Legacy/exploratory drafts and notes
│       ├── MANUSCRIPT_OUTLINE.md
│       ├── GBRT_DEBIAS_EXPLANATION.md
│       ├── GBRT_METHODS_FOR_MANUSCRIPT.md
│       └── GBRT_V2_TRAINING_REPORT.md
├── analysis/           # Python scripts for analysis and figure generation
│   ├── analyze_gtdb_quality_vs_mae.py
│   ├── plot_supplementary_simulations.py
│   ├── plot_gtdb_quality_combined.py
│   └── ...
├── data/               # Benchmark metadata and intermediate inputs
│   ├── benchmarks/
│   ├── gtdb_metadata/
│   └── ...
├── figures/            # Publication-ready figures
│   ├── report/         # Main-text and supplementary figures
│   ├── gtdb50k/
│   └── syntracker_validation/
├── results/            # Benchmark reports and per-pair data
│   ├── gtdb50k/
│   ├── mag_validation/
│   ├── sv_validation/
│   ├── synteny_bench/
│   └── ...
└── scripts/            # Analysis pipelines and HPC submission scripts
    ├── gtdb50k/
    ├── mag_validation/
    ├── syntracker_validation/
    └── ...
```

---

## Key Figures

Main-text figures are in `figures/report/` and supplementary figures in `figures/report/` and selected subdirectories.

| Figure | Description | Script / Source |
|--------|-------------|-----------------|
| Fig. 1 | Syn2bANI algorithm schematic | `figures/report/fig1_algorithm_schematic.png` |
| Fig. 2 | Exact-truth ANI ladder accuracy | `figures/report/fig1_simulation_ladder.png` |
| Fig. 3 | Robustness under indels, fragmentation, GC, accessory content | `figures/report/fig2_robustness.png` |
| Fig. 4 | Enzyme-panel optimization | `figures/report/fig3_enzyme_panel.png` |
| Fig. 5 | Mid-ANI validation against ANIm | `figures/report/fig4_midani_anim_validation.png` |
| Fig. 6 | Large-scale comparison against FastANI (45,000 GTDB pairs) | `figures/report/fig5_gtdb_r207_benchmark.png` |
| Fig. 7 | Computational efficiency | `figures/report/fig6_efficiency.png` |
| Fig. 8 | ANIm-truth benchmark by ANI band | `figures/report/fig7_anim_by_band.png` |
| Fig. 9 | Structural-variant detection on real genomes | `figures/report/fig8_sv_detection.png` |
| Fig. 10 | Accuracy on binned CAMI2 MAGs | `figures/report/mag_validation.png` |
| Fig. 11 | Near-clonal ANI masks extensive rearrangements | `figures/syntracker_validation/syntracker_high_ani_low_synteny.png` |
| Fig. 12 | Database-scale structurally divergent top hits | `figures/gtdb50k/gtdb_discordant_high_ani.png` |

Supplementary figures S1–S9 include the inversion-ladder truth benchmark, GTDB held-out and unified benchmarks, genome-quality robustness, and the exact-truth simulation families (indel, GC, fragmentation, accessory, mosaic).

---

## Data and Results

Per-pair benchmark data, ground-truth files, and summary reports are in `results/`:

| Dataset | Location | Ground truth | Notes |
|---------|----------|--------------|-------|
| GTDB-R207 calibration/training set (2,520 pairs) | `results/panel_by_band/` | dnadiff/ANIm | 2,074 band-stratified + 467 targeted 95–99.5% pairs |
| GTDB-R207 43,334 held-out pairs | `results/gtdb50k/` | dnadiff/ANIm | Strict genome-level holdout from calibration set |
| GTDB-R207 high-ANI test set (727 pairs) | `results/gtdb50k/high_ani_results.tsv` | dnadiff/ANIm | Non-representative genomes, 95–100% |
| Unified 80–100% benchmark | `results/gtdb50k/` | dnadiff/ANIm | 43,334 held-out + 727 high-ANI test |
| Mid-ANI validation | `results/validation/` | dnadiff/ANIm | 15 pairs, low alignment coverage |
| Oral/gut validation | `results/validation/` | dnadiff/ANIm / FastANI/skani | 50 isolates, 1,225 pairs |
| CAMI2 MAG benchmark | `results/mag_validation/` | dnadiff/ANIm + CAMI2 assignment | 695 bins |
| SV validation | `results/sv_validation/` | dnadiff structural | Enterobacteriaceae pairs |
| Synteny benchmark | `results/synteny_bench/` | Exact-truth inversion ladder | 0–32 inversions |
| Syntracker validation | `results/syntracker_validation/` | Structural re-analysis | *E. coli*, *H. pylori*, *N. gonorrhoeae*, *S. rimosus* |

---

## Reproducing Figures and Analyses

Most figures can be regenerated from the provided data using scripts in `analysis/`:

```bash
cd analysis
python3 analyze_gtdb_quality_vs_mae.py      # Fig. S4: genome quality vs ANI accuracy
python3 plot_supplementary_simulations.py   # Figs. S5–S9: simulation families
python3 plot_gtdb_quality_combined.py       # combined quality figure (Fig. S4)
```

The main Syn2bANI tool (Rust) lives in the [Syn2bANI code repository](https://github.com/HuangShiLab/Syn2bANI); the simulation harness is in its `prototype/` directory.

---

## Manuscript Status

The manuscript draft is `paper/manuscript.md` (and `.docx`). Key results are frozen for the current submission version:
- Default enzyme panel: BcgI, AlfI, AloI, FalI
- Calibration model: v5 (ridge regression on internal features)
- Main accuracy claim: MAE 0.619 on 39,903 held-out GTDB-R207 pairs
- Synteny/structural outputs validated against dnadiff and an exact-truth inversion ladder

---

## Related Repositories

- **[Syn2bANI](https://github.com/HuangShiLab/Syn2bANI)** — Main code repository (Rust)
- **[Syn2b](https://github.com/HuangShiLab/Syn2b)** — Synteny analysis tool (upstream)
- **[Fast2bRAD-M](https://github.com/HuangShiLab/Fast2bRAD-M)** — Fast tag extraction (methodology)

---

## License

Analysis scripts and data are released under the MIT License. The manuscript text is © 2025 HuangShiLab.
