# Figures Index

This directory lists the figures used in the current manuscript. The actual figure files live under `figures/` in the repository root.

## Main-text Figures

| Figure | Description | Source file |
|--------|-------------|-------------|
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

## Supplementary Figures

Supplementary figures S1–S11 include the inversion-ladder truth benchmark, GTDB held-out and unified benchmarks, genome-quality robustness, exact-truth simulation families (indel, GC, fragmentation, accessory, mosaic), and the Syntracker isolate-collection extensions. They are located in `figures/report/` and `figures/syntracker_validation/`.

| Figure | Description | Source file |
|--------|-------------|-------------|
| Fig. S1 | Exact-truth inversion-ladder synteny benchmark | `figures/report/fig_synteny_ladder.png` |
| Fig. S2 | Held-out GTDB-R207 benchmark | `figures/report/fig_gtdb50k_heldout.png` |
| Fig. S3 | Unified GTDB-R207 80–100% benchmark | `figures/gtdb_r207_unified_benchmark.png` |
| Fig. S4 | Genome quality vs. ANI accuracy | `figures/report/gtdb_quality_vs_mae_combined.png` |
| Fig. S5 | Exact-truth indel ladder and indel sweep | `figures/report/fig_s5_simulation_indel.png` |
| Fig. S6 | GC coverage ladder | `figures/report/fig_s6_simulation_gc.png` |
| Fig. S7 | Simulated fragmentation | `figures/report/fig_s7_simulation_fragment.png` |
| Fig. S8 | Accessory-content confound | `figures/report/fig_s8_simulation_accessory.png` |
| Fig. S9 | Mosaic/rate-heterogeneity family | `figures/report/fig_s9_simulation_mosaic.png` |
| Fig. S10 | Breakpoints vs. ANI in Syntracker isolate collections | `figures/syntracker_validation/syntracker_breakpoints_vs_ani.png` |
| Fig. S11 | *N. gonorrhoeae* and *S. rimosus* ANI–synteny decoupling | `figures/syntracker_validation/syntracker_supp_ngonorrhoeae_srimosus.png` |

## Regenerating Figures

Most figures can be regenerated from data using scripts in `analysis/`:

```bash
cd analysis
python3 analyze_gtdb_quality_vs_mae.py      # Fig. S4: genome quality vs ANI accuracy
python3 plot_supplementary_simulations.py   # Figs. S5–S9: simulation families
python3 plot_gtdb_quality_combined.py       # combined quality figure (Fig. S4)
```
