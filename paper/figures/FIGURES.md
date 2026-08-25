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

Supplementary figures S1–S9 include the inversion-ladder truth benchmark, GTDB held-out and unified benchmarks, genome-quality robustness, and the exact-truth simulation families (indel, GC, fragmentation, accessory, mosaic). They are located in `figures/report/` and selected subdirectories.

## Regenerating Figures

Most figures can be regenerated from data using scripts in `analysis/`:

```bash
cd analysis
python3 analyze_gtdb_quality_vs_mae.py      # Fig. S4: genome quality vs ANI accuracy
python3 plot_supplementary_simulations.py   # Figs. S5–S9: simulation families
python3 plot_gtdb_quality_combined.py       # combined quality figure (Fig. S4)
```
