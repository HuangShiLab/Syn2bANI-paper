# Task 0b Simulation Rerun Report

**Date:** 2026-08-25  
**Binary:** `/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani` (commit `fe0f36c`, calibration v5, default 4-enzyme panel BcgI/AlfI/AloI/FalI, rescue pass)  
**Harness:** `/Users/macstudio/Downloads/Syn2bANI/prototype/`  
**Outputs:** `/Users/macstudio/Downloads/Syn2bANI-paper/figures/report/fig_s*_simulation_*.png`

## 1. What was already up to date

The simulated genome FASTAs in `prototype/simindel/`, `prototype/simindel_sweep/`, `prototype/simfrag/`, and `prototype/simacc/` were present, but the corresponding `_4e.tsv` result tables were dated **2026-08-11**, while the release binary was built **2026-08-18**. Therefore all result TSVs were considered outdated and were regenerated with the current binary.

The GC-sweep source genomes (*F. nucleatum*, *S. mutans*, *B. longum*, *S. coelicolor*) are **not present** in the repository, so the GC sweep could not be rerun from scratch (see §4).

## 2. What was rerun

A single runner script was created at `prototype/run_simulations_4e.py` to regenerate the simulation genomes (where deterministic) and run `syn2bani ani` with the default 4-enzyme panel, 8 threads, and `--verbose`.

| Family | Script | Result TSV | Notes |
|--------|--------|------------|-------|
| Indel ladder | `simulate.py mg1655.fasta simindel 1.0` | `simindel_results_4e.tsv` | 12 ANI levels (0.85–0.999), inversion + 1 deletion/100 kb |
| Indel sweep | `simulate_indel_sweep.py mg1655.fasta simindel_sweep` | `simindel_sweep_results_4e.tsv` | Fixed 95% ANI, indel rates 0–4 per 100 kb |
| Fragmentation | `simulate.py mg1655.fasta simfrag_tmp_src 0.0` → `fragment.py q95.fasta simfrag` | `simfrag_results_4e.tsv` | 95% ANI draft fragmented to 20/50/100/200 contigs, flipped/shuffled |
| Accessory | `simulate_accessory.py mg1655.fasta simacc 0.95` | `simacc_results_4e.tsv` | Fixed 95% core ANI, accessory fractions 0–50% |
| Mosaic | `simulate_mosaic.py mg1655.fasta simmosaic 5` | `simmosaic_results_4e.tsv` | Gamma + bimodal rate heterogeneity, 5 kb blocks |

All commands were:

```bash
syn2bani ani --ql <queries.txt> --rl <refs.txt> -p -t 8 --verbose
```

The runner did **not** use `--calibrate`, so the reported `ani` / `ani_uniform` columns are the raw chain-restricted MLE estimates (gated by the new rescue logic). Exact-truth labels are preserved in the per-simulation `manifest.tsv` files.

## 3. Key MAE / accuracy numbers

| Family | n | MAE gamma (ANI points) | MAE uniform (ANI points) | Max |error| gamma |
|--------|---|------------------------|--------------------------|---------------------|
| Indel ladder | 12 | 0.0732 | 0.0732 | 0.2068 |
| Indel sweep | 5 | 0.0807 | 0.0807 | — |
| Fragmentation | 4 | 0.0695 | 0.0695 | — |
| Accessory | 6 | 0.1121 | 0.1121 | 0.2588 |
| Mosaic (gamma) | 6 | 1.3541 | 2.7507 | — |
| Mosaic (bimodal) | 3 | 2.4252 | 3.4247 | — |

Observations:

* The indel ladder, indel sweep, fragmentation, and accessory families all stay well below 0.12 MAE, confirming the current 4-enzyme panel is robust to these confounds.
* The mosaic family is harder: even the gamma estimator is biased high (mean +0.5 to +2.4 ANI points depending on regime), though it consistently outperforms the uniform estimator. The bimodal misspecification case shows the largest deviation.

## 4. Figures produced

| Figure | File(s) | Data source | Notes |
|--------|---------|-------------|-------|
| S5 Indel | `figures/report/fig_s5_simulation_indel.{png,pdf}` | `simindel_results_4e.tsv` + `simindel_sweep_results_4e.tsv` | Two panels: ANI ladder scatter + indel-sweep error |
| S6 GC | `figures/report/fig_s6_simulation_gc.{png,pdf}` | Historical GC sweep from `ALGORITHM_MLE.md` §4.8 | **Not rerun**: source genomes missing; values are the previously reported 4-enzyme vs 5-enzyme sweep |
| S7 Fragmentation | `figures/report/fig_s7_simulation_fragment.{png,pdf}` | `simfrag_results_4e.tsv` | Error vs contig count (log x-axis) |
| S8 Accessory | `figures/report/fig_s8_simulation_accessory.{png,pdf}` | `simacc_results_4e.tsv` | ANI error + aligned fraction vs accessory fraction |
| S9 Mosaic | `figures/report/fig_s9_simulation_mosaic.{png,pdf}` | `simmosaic_results_4e.tsv` | Error vs true ANI for gamma and bimodal regimes |

**Figure numbering.** Renumbered to `S5`–`S9` to avoid conflict with the existing supplementary figures S1–S4 (inversion ladder, GTDB held-out, unified benchmark, quality vs MAE).

## 5. Errors / blockers

* **GC sweep source genomes missing.** Only *E. coli* K-12 MG1655 is present in the repo. The other four GC-sweep species (*F. nucleatum*, *S. mutans*, *B. longum*, *S. coelicolor*) are not available, so the GC panel was plotted from the historical hard-coded table rather than rerun. If the source FASTAs can be provided, `prototype/gc_bench.py` can be rerun to update the figure.
* No other errors. All simulations completed in <1 minute total.

## 6. Files touched

In `/Users/macstudio/Downloads/Syn2bANI`:

* `prototype/run_simulations_4e.py` (new runner script)
* `prototype/simindel_results_4e.tsv` (regenerated)
* `prototype/simindel_sweep_results_4e.tsv` (regenerated)
* `prototype/simfrag_results_4e.tsv` (regenerated)
* `prototype/simacc_results_4e.tsv` (regenerated)
* `prototype/simmosaic_results_4e.tsv` (new)
* `prototype/simmosaic/` (new directory with simulated genomes)
* `prototype/simfrag_tmp_src/` (new temporary source directory)

In `/Users/macstudio/Downloads/Syn2bANI-paper`:

* `analysis/plot_supplementary_simulations.py` (new plotting script)
* `figures/report/fig_s5_simulation_indel.{png,pdf}`
* `figures/report/fig_s6_simulation_gc.{png,pdf}`
* `figures/report/fig_s7_simulation_fragment.{png,pdf}`
* `figures/report/fig_s8_simulation_accessory.{png,pdf}`
* `figures/report/fig_s9_simulation_mosaic.{png,pdf}`
* `results/SIMULATION_RERUN_REPORT.md` (this file)

No commits or pushes were performed.
