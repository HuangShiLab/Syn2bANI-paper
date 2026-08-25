# Results Archive

This directory contains superseded intermediate files from earlier stages of the project. They are kept for reproducibility but are not part of the current manuscript's main results.

## Contents

- **GBRT v3/v4 models and reports** — Older gradient-boosted calibration attempts (`gbrt_model_v3_*`, `gbrt_model_v4_*`, `gbrt_v3_*`, `gbrt_v4_*`). Superseded by the ridge calibration v5/v6.
- **Old evaluation and matrix files** — Small-scale (100/1k/10k) and v7 evaluation outputs (`evaluation_gtdb_r207_*`, `matrix_gtdb_r207_*`, `gtdb_r207_100k_*`). Superseded by the 43,334-pair held-out benchmark and v8 final matrix.
- **Old pair sampling files** — Early pair lists (`pairs_gtdb_r207_*`, `pairs_test*`). Superseded by `results/gtdb50k/pairs_50k.tsv` and related files.
- **Old benchmark logs** — `benchmark_*.log` from early timing experiments.
- **Old benchmark reports** — `BENCHMARK*.md`, `COMPREHENSIVE_BENCHMARK.md`, `FINAL_BENCHMARK_REPORT.md`, `FINAL_REPORT_v0.1.1.md`, `HEAD_TO_HEAD_REPORT.md`, `GTDB_R207_BENCHMARK_REPORT.md`, `GTDB_R207_V8_BENCHMARK_REPORT.md`, `PERFORMANCE_BENCHMARK.md`, `HPC_SCALING_PLAN.md`.
- **Old calibration files** — `panel_by_band/` subdirectory containing v2/v3/v4 ridge and linear calibration files. Superseded by v5/v6 calibration files in `results/panel_by_band/`.
- **Old truth files** — Gated/v8 truth versions for the 2,074-pair and hi95 sets. Superseded by the v9-rescue truth files in `results/`.
- **Other test/temporary files** — `test_*.tsv/json`, `sample_anim_truth.tsv`, `oral_gut_1225_v8current.tsv`, `oral_gut_validation_v8_report.md`, `COMPREHENSIVE_BENCHMARK.md`.

## Current results

The main results used in the manuscript remain in `results/` and its subdirectories:

- `results/panel_by_band/` — v5/v6 calibration and main ANIm-truth table
- `results/gtdb50k/` — 43,334 held-out pairs, 727 high-ANI test pairs, unified benchmark
- `results/validation/` — mid-ANI validation
- `results/mag_validation/` — CAMI2 MAG benchmark
- `results/syntracker_validation/` — Syntracker discordant-pair analysis
- `results/synteny_bench/` — exact-truth inversion-ladder synteny benchmark
- `results/sv_validation/` — structural-variation validation on Enterobacteriaceae
- `results/efficiency_v8/` — computational-efficiency benchmarks
- `results/db_scale/` — database-scale dist/search/triangle scaling
