# Scripts Index

This directory documents the analysis pipelines and scripts used in the current manuscript. The actual scripts live under `analysis/` and `scripts/` in the repository root.

## Analysis Scripts (`analysis/`)

Scripts for figure generation, model training, and result summaries:

- `analyze_gtdb_quality_vs_mae.py` — Genome quality vs. ANI accuracy (Fig. S4)
- `plot_supplementary_simulations.py` — Supplementary simulation figures (Figs. S5–S9)
- `plot_gtdb_quality_combined.py` — Combined quality figure
- `plot_final.py` — Main-text figure assembly
- `plot_comparison.py` — Tool comparison plots
- `export_gbrt.py` — GBRT model export utilities
- `benchmark_pipeline.py` — Benchmark orchestration
- `multienzyme_benchmark.py` — Multi-enzyme panel benchmarks
- `performance_benchmark.py` — Efficiency benchmarking

## Pipeline Scripts (`scripts/`)

End-to-end benchmark and HPC submission scripts:

- `gtdb50k/` — GTDB-R207 43,334-pair held-out benchmark pipeline
- `mag_validation/` — CAMI2 MAG validation pipeline
- `syntracker_validation/` — Syntracker discordant-pair analysis
- `synteny_bench/` — Inversion-ladder synteny benchmark
- `hi95_anim/` — High-ANI ANIm truth generation
- `anim_mid_ani/` — Mid-ANI oral/gut validation
- `evaluate_vs_anim.py` — Evaluate Syn2bANI against ANIm truth
- `train_gbrt_v5.py` — Train calibration v5 model
- `evaluate_gbrt_v5.py` — Evaluate calibration v5 model
- `slurm_train_final.sh` — Final HPC training submission
- `slurm_matrix_v8_array.sh` — Feature-matrix extraction on GTDB-R207
- `slurm_fastani_v8_array.sh` — FastANI reference runs on GTDB-R207

## Legacy / Exploratory Scripts

Earlier GBRT-era and prototype scripts remain in `scripts/` for reproducibility but are not part of the current manuscript's main pipeline. They are superseded by the v8 chain-restricted MLE estimator and calibration v5 workflow.
