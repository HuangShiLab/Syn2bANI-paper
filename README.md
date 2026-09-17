# Syn2bANI-paper

> Paper repository for **"Syn2b-ANI: Strain-level ANI estimation and structural comparison via fixed restriction-site anchors"**

This repository contains the manuscript, analysis code, benchmark data, and figures for the Syn2bANI paper.

---

## Repository Structure

```
Syn2bANI-paper/
├── paper/
│   ├── manuscript/
│   │   ├── manuscript_nature_methods.md   # CURRENT manuscript
│   │   ├── SUPPLEMENTARY.md               # supplementary figures and tables
│   │   ├── manuscript.md                  # superseded long-form draft (do not cite)
│   │   └── *.docx
│   ├── figures/
│   │   ├── main/                          # Figs. 1-8
│   │   ├── supplementary/                 # Figs. S1-S12
│   │   ├── archive/                       # withdrawn figures
│   │   └── FIGURES.md                     # index (stale numbering; main/ is authoritative)
│   ├── REVIEW_2026-09-14.md               # internal pre-submission review
│   ├── data/, scripts/, others/           # indices and legacy notes
├── analysis/           # analysis and figure scripts
│   ├── breakpoint_ladder.py               # breakpoint_count regression check
│   ├── calibrated_comparators.py          # skani/FastANI under the same calibration
│   └── ...
├── case_studies/
│   ├── rerun_breakpoints.py               # all-vs-all re-run with a current binary
│   ├── fetch_ncbi_fna.py                  # genome fetch via the Datasets v2 API
│   ├── ecoli_o157_fitzgerald_2021/, fda_argos_s_aureus/, h_pylori_cagpai/
├── data/               # benchmark metadata and intermediate inputs
├── figures/            # generated figure output (gtdb50k/, syntracker_validation/, report/)
├── results/            # per-pair data and benchmark reports
└── scripts/            # pipelines and HPC submission scripts
```

---

## Key Figures

Main-text figures are in `paper/figures/main/`, supplementary figures in
`paper/figures/supplementary/`. The numbering below is the one used by
`paper/manuscript/manuscript_nature_methods.md`.

| Figure | Description | File |
|--------|-------------|------|
| Fig. 1 | The Syn2bANI estimator | `paper/figures/main/fig1_algorithm.png` |
| Fig. 2 | Accuracy and robustness under exact truth | `paper/figures/main/fig2_simulations.png` |
| Fig. 3 | Unified GTDB-R207 80–100% benchmark vs ANIm | `paper/figures/main/fig3_gtdb_r207_benchmark.png` |
| Fig. 4 | Computational efficiency | `paper/figures/main/fig4_efficiency.png` |
| Fig. 5 | Structural outputs vs alignment-based truth | `paper/figures/main/fig5_sv_validation.png` |
| Fig. 6 | Structural divergence at near-clonal ANI, four isolate collections | `paper/figures/main/fig6_syntracker_structure.png` |
| Fig. 7 | Locus-targeted chain coverage (*B. longum* abfA) | `paper/figures/main/fig7_b_longum_abfa.png` |
| Fig. 8 | cagPAI architecture in 528 *H. pylori* isolates | `paper/figures/main/fig8_cagpai_h_pylori.png` |

Withdrawn: `fig6_ani_sv_discordance` (GTDB-R207 discordant pairs, 61–989 SV
calls) and `fig_s8_syntracker_breakpoints` were produced with a pre-c974f5f
`breakpoint_count` and have been removed; see **Manuscript Status** below.

Supplementary figures S1–S12 cover the inversion-ladder truth benchmark, the GTDB held-out and per-band benchmarks, genome-quality robustness, the exact-truth simulation families (GC, mosaic), the CAMI2 MAG benchmark, the *H. pylori* participant resolution (S8), the *E. coli* O157:H7 (S9) and FDA-ARGOS *S. aureus* (S10) collections, and the cagPAI country/population breakdown (S11) and circular-origin filtering (S12). Legends are in `paper/manuscript/SUPPLEMENTARY.md`.

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
# Supplementary figures from the shipped per-pair data
python3 analysis/analyze_gtdb_quality_vs_mae.py     # Fig. S4
python3 analysis/plot_supplementary_simulations.py  # simulation families
python3 analysis/plot_gtdb_quality_combined.py      # combined quality figure

# Comparators under the same band-holdout calibration Syn2bANI gets
python3 analysis/calibrated_comparators.py          # results/gtdb50k/CALIBRATED_COMPARATORS.md

# breakpoint_count regression check (needs the syn2bani binary and a MG1655 FASTA)
python3 analysis/breakpoint_ladder.py \
    --genome ../Syn2bANI/prototype/mg1655.fasta \
    --syn2bani ../Syn2bANI/target/release/syn2bani

# Fig. 6 / S8: SynTracker cohorts from the shipped per-pair structural table
python3 scripts/syntracker_validation/09_analyze_structural_vs_ani.py \
    --structural data/syntracker_validation/syn2b_structural_raw/syn2b_structural_pairs_raw.tsv \
    --skani-dir data/syntracker_validation/skani \
    --metadata-dir data/syntracker_validation/samples \
    --outdir figures/syntracker_validation

# Figs. S9 / S10: case-study re-runs (downloads genomes first; *.fna is gitignored)
python3 case_studies/fetch_ncbi_fna.py \
    --accessions case_studies/ecoli_o157_fitzgerald_2021/scripts/accession_map.tsv \
    --outdir case_studies/ecoli_o157_fitzgerald_2021/genomes
python3 case_studies/rerun_breakpoints.py --study ecoli_o157_fitzgerald_2021 \
    --metadata results/metadata_with_lineage.tsv \
    --group-cols assigned_lineage,host_category
python3 case_studies/rerun_breakpoints.py --study fda_argos_s_aureus \
    --metadata results/assembly_metadata.tsv --group-cols country --dedup-stem
```

The main Syn2bANI tool (Rust) lives in the [Syn2bANI code repository](https://github.com/HuangShiLab/Syn2bANI); the simulation harness is in its `prototype/` directory.

---

## Manuscript Status

The current manuscript is `paper/manuscript/manuscript_nature_methods.md`.
`paper/manuscript/manuscript.md` is the superseded long-form draft and carries a
banner to that effect. Frozen choices for this submission version:

- Default enzyme panel: BcgI, AlfI, AloI, FalI
- Calibration model: v5 (ridge regression on internal features)
- Main accuracy claim: MAE 0.619 on the 39,903 of 43,334 held-out GTDB-R207
  pairs for which the calibrated estimator returns a value; with the same
  band-holdout linear calibration applied to the comparators, skani reaches
  0.75 and FastANI 0.55, so the gain over uncalibrated tools is largely the
  removal of a constant bias (`results/gtdb50k/CALIBRATED_COMPARATORS.md`)
- `breakpoint_count` was redefined in Syn2bANI v0.1.1 (block-level,
  positive-contradiction; see that repository's README). Every structural
  number in the current manuscript was recomputed with it, and
  `analysis/breakpoint_ladder.py` is the regression check. Three classes of
  result produced with earlier builds are withdrawn: the GTDB-R207
  ANI-vs-SV discordance figure, the pre-fix SynTracker summary
  (`results/syntracker_validation/syntracker_summary_pre_fix.tsv`, retained for
  the record), and the `breakpoint_count`-vs-dnadiff partial correlations on
  the 43,334-pair set

Open items: the simulation-ladder and Enterobacteriaceae skani comparisons
still use 0.1.0 (caveat stated in the manuscript; the cohort ANI axis has been
re-run with 0.3.2 and reproduces exactly), and the three-class cagPAI locality
variable has not been run through the stratified disease-stage tests because
that needs the cohort metadata. **Closed 2026-09-18:** the 43,334-pair
structural re-computation and the `breakpoint_count`-vs-dnadiff correlations
were completed on the HPC with the v0.1.1 binary (partial r = 0.411 vs the
withdrawn 0.414; large-indels 0.327 vs 0.453; ANIm ≥ 95% ρ = 0.633 vs 0.674;
plus a new cross-implementation concordance with Syn2b junctions, ρ = 0.86 on
43,312 pairs — outputs in the Syn2b-paper repository,
`results/gtdb50k/rerun_v011/` and `results/gtdb50k/sv_reanalysis_v011.md`),
and the four SynTracker cohorts were re-run through Syn2bANI's own structural
channel (medians 0/1/5/3; within-*H. pylori* 0 vs between 6; per-pair data in
`results/syntracker_validation/rerun_v032/` of the Syn2b-paper repository).
The withdrawn GTDB discordance figure is superseded by the Syn2b-paper
within-species census (in progress).

---

## Related Repositories

- **[Syn2bANI](https://github.com/HuangShiLab/Syn2bANI)** — Main code repository (Rust)
- **[Syn2b](https://github.com/HuangShiLab/Syn2b)** — Synteny analysis tool (upstream)
- **[Fast2bRAD-M](https://github.com/HuangShiLab/Fast2bRAD-M)** — Fast tag extraction (methodology)

---

## License

Analysis scripts and data are released under the MIT License. The manuscript text is © 2025 HuangShiLab.
