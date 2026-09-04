# Syn2bANI Manuscript & Figure Reorganization Outline

**Goal:** One clean logical flow built around a single authoritative dataset (GTDB-R207 50k), with minimal intermediate/exploratory figures and each case study collapsed to one composite figure.

**Core storyline**
1. What Syn2bANI does and why (algorithm schematic).
2. ANI accuracy: exact-truth simulations → real-genome GTDB-R207 benchmark → hybrid estimator.
3. Computational efficiency: same GTDB-R207 dataset, pairwise and database scale.
4. SV outputs: validation against dnadiff/minimap2 on GTDB-R207, then ANI+SV joint view.
5. Case studies: each is one composite figure showing ANI vs SV + locus/strain-level detail.

---

## 1. Main Text Figures

### Figure 1 | Algorithm overview and outputs
- **Content:** in-silico Type IIB digestion → shared anchors → chaining → MLE → ANI + AF + synteny/SV outputs.
- **Source:** merge current `figures/report/fig1_algorithm_schematic.png` with the output-panel concept; keep simple.
- **Action:** keep / lightly relabel.

### Figure 2 | ANI accuracy under exact truth and main robustness checks
- **Content:** single composite with four sub-panels:
  - (a) ANI ladder (85–99.9%) vs true ANI — Syn2bANI raw, skani, FastANI.
  - (b) Indel sweep at 95% ANI.
  - (c) Fragmentation ladder (20–201 contigs).
  - (d) Accessory-content sweep (0–50% shuffled blocks).
- **Rationale:** establishes that raw estimator is unbiased and robust before adding calibration.
- **Action:** merge `fig1_simulation_ladder.png` + `fig2_robustness.png` into one 4-panel figure; drop redundant GC ladder from main text (move to supplementary if kept).

### Figure 3 | GTDB-R207 ANI benchmark: the central result
- **Content:** single composite:
  - (a) 43,334 held-out same-genus pairs: estimated vs ANIm truth (Syn2bANI hybrid, skani, FastANI).
  - (b) Per-band MAE (80–85, 85–90, 90–95, 95–97, 97–100).
  - (c) Signed error distributions.
  - (d) High-ANI zoom (95–100%) showing raw gated vs v6-calibrated vs hybrid.
- **Rationale:** this is the paper's headline accuracy figure; all other ANI accuracy plots fold into it.
- **Action:** use `figures/gtdb_r207_unified_benchmark.png` as base; if needed, regenerate as 4-panel composite. Drop separate `fig7_anim_by_band.png` and `fig_gtdb50k_heldout.png` from main text (can become Fig. S2/S3).

### Figure 4 | Computational efficiency on GTDB-R207
- **Content:** single composite:
  - (a) Wall time vs number of genomes for pairwise `dist`, `search`, and `triangle` on GTDB-R207 subsamples.
  - (b) Peak memory per tool.
  - (c) Throughput (comparisons per second) at database scale.
  - (d) Time breakdown: sketch vs compare vs I/O.
- **Rationale:** efficiency claims must be anchored to the same GTDB-R207 dataset as accuracy.
- **Action:** regenerate from `results/gtdb50k/efficiency_v8/` or `results/efficiency_v8/` data; current `figures/report/fig6_efficiency.png` may be outdated.

### Figure 5 | SV outputs agree with alignment-based truth
- **Content:** single composite:
  - (a) `raw_inverted_fraction` vs dnadiff inverted fraction (all 43,334 pairs; Pearson r = 0.936).
  - (b) Strain-level zoom ANIm ≥ 97% (r = 0.996).
  - (c) `breakpoint_count` vs dnadiff relocations+inversions (partial correlation controlling ANIm + contig count).
  - (d) `af_query` vs dnadiff/minimap2 aligned fraction.
- **Rationale:** validates that SV outputs are meaningful; puts `raw_inverted_fraction` as the primary metric.
- **Action:** use `figures/report/fig_inverted_fraction_comparison_high_ani.png` and `results/gtdb50k/sv_reanalysis_metrics.tsv`; regenerate as unified 4-panel figure.

### Figure 6 | ANI and SV together: discordant high-ANI pairs in GTDB-R207
- **Content:** single composite:
  - (a) ANIm vs `raw_inverted_fraction` for GTDB-R207 pairs, color-coded by ANI band.
  - (b) Examples of pairs with ANI > 99% but high inverted fraction; call out species.
  - (c) ANI vs breakpoint count for the same pairs.
  - (d) A representative pair shown as chain/SV calls.
- **Rationale:** demonstrates the biological value of the dual-axis view.
- **Action:** generate from `results/gtdb50k/discordant_*.tsv` and `closed_inversion_pairs_annotated.tsv`.

### Figure 7 | B. longum abfA case study
- **Content:** single composite figure:
  - (a) abfA locus diagram / chain coverage across locus.
  - (b) Periphery coverage boxplot by cluster status (complete vs deleted/truncated).
  - (c) ANI distribution between groups (showing no separation).
  - (d) Selected strain phylogeny or pairwise ANI + SV summary.
- **Rationale:** one figure tells the entire case-study story.
- **Action:** merge current `figures/report/fig_b_longum_abfa.png` panels + add ANI/SV summary.

### Figure 8 | H. pylori cagPAI case study
- **Content:** single composite figure:
  - (a) Extended cagPAI state classification workflow.
  - (b) Stacked-bar distribution by FastBAPS lineage.
  - (c) Stacked-bar distribution by disease stage (for comparison, with CMH note).
  - (d) Engineered control validation (WT, Δ, inv, transloc).
- **Rationale:** collapse current Fig. 8 + Fig. S12 + S19–S20 into one main figure.
- **Action:** regenerate from filtered association data; supplementary only needs per-country/population breakdown if space permits.

---

## 2. Supplementary Figures

### Supplementary Figure S1 | Exact-truth inversion ladder
- `figures/report/fig_synteny_ladder.png` — keep, rename from `fig_synteny_ladder` to `fig_s1_synteny_ladder` if missing file needs to be restored.

### Supplementary Figure S2 | Held-out GTDB-R207 43k benchmark (v5 detail)
- Move current main-text `figures/report/fig_gtdb50k_heldout.png` here; shows the v5-only benchmark.

### Supplementary Figure S3 | ANI by band and tool comparison
- Move current main-text `figures/report/fig7_anim_by_band.png` here.

### Supplementary Figure S4 | Genome quality does not drive ANI error
- `figures/report/gtdb_quality_vs_mae_combined.png` — keep.

### Supplementary Figure S5 | GC coverage ladder
- `figures/report/fig_s6_simulation_gc.png` — keep as supplementary.

### Supplementary Figure S6 | Mosaic/rate-heterogeneity family
- `figures/report/fig_s9_simulation_mosaic.png` — keep.

### Supplementary Figure S7 | MAG validation
- `figures/report/mag_validation.png` — keep; currently main-text Fig. 5 candidate but better as supplementary unless MAG accuracy is a central claim.

### Supplementary Figure S8 | Syntracker isolate collections: ANI vs breakpoint count
- `figures/syntracker_validation/syntracker_ani_vs_breakpoints.png` — keep; currently main-text Fig. 6 but may move to supplementary if GTDB-R207 discordance (Fig. 6 above) becomes the main structural decoupling figure.

### Supplementary Figure S9 | E. coli O157:H7 case-study composite
- Merge lineage + host source breakpoint plots into one figure.
- Source: `figures/report/fig_s15_ecoli_o157_breakpoints_*.png`.

### Supplementary Figure S10 | FDA-ARGOS S. aureus case-study composite
- Merge country + source breakpoint plots into one figure.
- Source: `figures/report/fig_s16_saureus_breakpoints_*.png`.

### Supplementary Figure S11 | H. pylori cagPAI extended state by country and phylogenetic population
- Merge `fig_s17_cagpai_by_country.png` + `fig_s18_cagpai_by_population.png` into one figure.

### Supplementary Figure S12 | Breakdown of circular-origin artifact filtering in H. pylori
- New figure: before/after cagPAI state counts + example genome-spanning call.

---

## 3. Tables

### Table 1 | GTDB-R207 ANI benchmark MAE by band
- Keep current table in manuscript; update numbers after hybrid estimator is finalized.

### Table 2 | Computational efficiency summary
- New table: tool, n genomes, n pairs, wall time, peak memory, throughput, notes (skani excludes filtered pairs).

### Table 3 | SV output correlations with dnadiff/minimap2 truth
- Consolidate `sv_reanalysis_metrics.tsv` into a clean table: metric, vs dnadiff/minimap2, raw r, partial r (controlling ANIm + contig), note.

### Table 4 | Case-study summary
- One summary table: dataset, n genomes / pairs, ANI range, key SV metric, biological association, figure.

### Supplementary Table S1 | MAG accuracy by quality tier
- Keep.

### Supplementary Table S2 | Lineage-stratified cagPAI CMH results
- Keep current Supplementary Table S8; renumber to S2 if tables are consolidated.

### Supplementary Table S3 | Hybrid threshold sweep
- Keep current table in Supplementary Note 6; extract as formal table if desired.

---

## 4. File Actions

### Keep / rename in place
- `figures/report/fig1_algorithm_schematic.png` → Fig. 1
- `figures/report/fig_b_longum_abfa.png` → base for Fig. 7
- `figures/report/fig_inverted_fraction_comparison_high_ani.png` → base for Fig. 5
- `figures/gtdb_r207_unified_benchmark.png` → base for Fig. 3
- `figures/report/gtdb_quality_vs_mae_combined.png` → Fig. S4

### Merge into composite main figures
- `fig1_simulation_ladder.png` + `fig2_robustness.png` → Fig. 2
- cagPAI panels + engineered controls → Fig. 8

### Move from main to supplementary
- `fig_gtdb50k_heldout.png` → Fig. S2
- `fig7_anim_by_band.png` → Fig. S3
- `mag_validation.png` → Fig. S7 (unless kept as main)
- `syntracker_ani_vs_breakpoints.png` → Fig. S8 (unless kept as main)

### Merge case-study supplementary figures
- `fig_s15_ecoli_o157_breakpoints_lineage.png` + `fig_s15_ecoli_o157_breakpoints_host.png` → Fig. S9
- `fig_s16_saureus_breakpoints_country.png` + `fig_s16_saureus_breakpoints_source.png` → Fig. S10
- `fig_s17_cagpai_by_country.png` + `fig_s18_cagpai_by_population.png` → Fig. S11

### Delete or archive to `paper/others/`
- Old gamma/uniform comparison figures no longer referenced.
- Intermediate calibration figures (pre-v5, pre-v6) unless historically needed.
- Duplicative GTDB benchmark panels that are folded into Fig. 3.
- `figures/report/fig_diagnostic_af_query_vs_anchor_adjacency.*` — diagnostic, not for paper.
- `figures/report/fig_inverted_fraction_by_band.*` if redundant with Fig. 5/6.
- `figures/report/fig_inverted_fraction_comparison_high_ani_raw.*` / `_min.*` — variants; keep only the chosen one.

### Needs to be generated
- Fig. 4 (efficiency composite).
- Fig. 5 (SV validation composite).
- Fig. 6 (GTDB-R207 ANI–SV discordance composite).
- Fig. 7 composite (B. longum).
- Fig. 8 composite (cagPAI).
- Fig. S12 (circular-origin filtering demonstration).

---

## 5. Results Section Order (proposed)

1. **The Syn2bANI estimator** — Fig. 1.
2. **Accurate ANI estimation under exact truth** — Fig. 2.
3. **Real-genome accuracy on GTDB-R207** — Fig. 3 + Table 1.
4. **Computational efficiency** — Fig. 4 + Table 2.
5. **Structural outputs and their validation** — Fig. 5 + Table 3.
6. **ANI and structure are decoupled in GTDB-R207** — Fig. 6.
7. **Case study: B. longum abfA** — Fig. 7.
8. **Case study: H. pylori cagPAI** — Fig. 8 + Table 4.
9. **Discussion**
10. **Online Methods**

---

## 6. Repository Organization (proposed)

```
paper/
  manuscript/
    manuscript_nature_methods.md
    manuscript_nature_methods.docx
    SUPPLEMENTARY.md
    OUTLINE_REORGANIZED.md   (this file)
  figures/
    main/                    # final main-text figures, numbered
      fig1_algorithm.png
      fig2_simulations.png
      ...
    supplementary/           # final supplementary figures, numbered
      fig_s1_...
      ...
    archive/                 # old/exploratory figures (keep but not linked)
  data/
    tables/                  # final TSV/CSV tables referenced in paper
  scripts/
    generate_figures.py      # one script to regenerate all final figures
```
