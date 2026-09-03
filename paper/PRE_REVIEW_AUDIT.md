# Pre-review audit — Syn2bANI manuscript

## Manuscript status

- **Format**: Nature Methods markdown draft (`paper/manuscript/manuscript_nature_methods.md`)
- **Word count**: Abstract ~153 words; main text ~4,500 words (within target).
- **Sections**: Abstract, Introduction, Results (8 subsections), Discussion, Online Methods, Data/Code availability, References, Figure Legends, Tables.
- **Placeholders remaining**: Acknowledgements, Author contributions, Competing interests.

## Structural completeness

| Section | Status | Notes |
|---|---|---|
| Title / authors / affiliations | OK | Correspondence uses shihuang@hku.hk; line 267 still has old placeholder `huangshi@njau.edu.cn`. |
| Abstract | OK | Could mention structural outputs and case studies more explicitly. |
| Introduction | OK | Clear motivation and differentiation from k-mer tools. |
| Results | OK | 8 subsections; logic flows from estimator → simulations → real benchmarks → structural outputs → case studies → efficiency. |
| Discussion | OK | Three boundaries are well stated; could add a dedicated "Limitations" subheading. |
| Methods | OK | Very detailed; may be too long for Nature Methods (consider moving some derivations to Supplementary). |
| Data/Code availability | OK | Repo links present; commit `fe0f36c` should be verified as current. |
| References | CHECK | Ref 10 (Syntracker) year/format needs verification. Ref 16 (dRep) is cited in text? Verify all numbered refs are cited. |

## Critical issues

### 1. Figure file naming vs. figure numbers

The PNG files in `figures/report/` do **not** follow main-text figure numbers; captions instead reference source files by descriptive name. This is acceptable for assembly, but risks confusion. Recommended: rename or create a `figures/submission/` folder with `fig1.png`–`fig8.png` copies mapped to captions.

### 2. Missing supplementary figures

The following supplementary figures are referenced but not present in the repository:

| Figure | Caption location | Needed content | Current status |
|---|---|---|---|
| S10 | SUPPLEMENTARY.md:114 | `syntracker_breakpoints_vs_ani.png` — E. coli + H. pylori ANI vs breakpoints | **MISSING** |
| S11 | SUPPLEMENTARY.md:125 | `syntracker_supp_ngonorrhoeae_srimosus.png` — N. gonorrhoeae + S. rimosus decoupling | **MISSING** |
| S14 | manuscript:89 | inverted_fraction standard-error model / comparison (referenced in main text) | **MISSING** |
| S15 | manuscript:289 | E. coli O157:H7 ANI vs breakpoints by lineage/host | **MISSING** (source files exist in `case_studies/ecoli_o157_fitzgerald_2021/figures/`) |
| S16 | manuscript:291 | FDA-ARGOS S. aureus breakpoint distribution | **MISSING** (source files exist in `case_studies/fda_argos_s_aureus/figures/`) |
| S17–S20 | manuscript:287 | cagPAI extended state by country / population / FastBAPS / disease-stage stacked bars | **PARTIAL** (4 PNGs exist in `case_studies/h_pylori_cagpai/results/struct_extended/` but not linked as S17–S20) |

**Action**: generate/copy the missing panels and assign explicit S10–S20 filenames.

### 3. Table 3 inconsistencies

- Row 10: "GTDB scale set" says used in "Fig. S3", but S3 is the unified 80–100% benchmark. The 45,000 FastANI comparison may not be shown in any main/supplementary figure; clarify or remove reference.
- Row 11: "Unified high-ANI set" says 2,342 pairs, but text says 727 high-ANI test pairs. Reconcile.
- Row 15: "Real drafts" says "Fig. 5", but Fig. 5 is CAMI2 MAGs. Should likely be Fig. 2.
- Row 17/18: GTDB 50k held-out/SV benchmark — 43,334 pairs in text vs 45,000 in row 10; make consistent.

### 4. Correspondence email inconsistency

Line 11: `shihuang@hku.hk`  
Line 267: `huangshi@njau.edu.cn` (old placeholder)  
**Action**: update line 267.

### 5. Reference 10 citation

Ref 10: Enav, Paz & Ley, *Nat. Biotechnol.* (2024). Verify exact volume/page/year; the paper may have been published in 2025 or as online early. Also confirm the in-text citation at line 79 is the intended Syntracker reference.

## Technical / claim-level risks for peer review

1. **skani is faster and more accurate near-clonal** — stated openly, but reviewers may ask why one would use Syn2bANI if skani is better in the most common use case (strain-level). The hybrid rule addresses this; make sure it is emphasized in the Abstract.

2. **Calibration is regime-specific** — this is disclosed, but reviewers may worry about users applying `--calibrate` to MAGs. Consider adding a CLI warning or default behavior discussion.

3. **Structural outputs vs. base-level SV tools** — the Discussion contrasts with minigraph-Cactus/SyRI/Mauve. Reviewers may ask for direct comparison on SV recall/precision. Current validation is correlation-based; coordinate-level validation on real pairs is weak (see `results/closed_inversions/JUNCTION_COORDINATE_REPORT.md`).

4. **anchor_adjacency interpretation** — the manuscript carefully says it is not a proxy for base-pair synteny. Reviewers may still challenge its utility. The case studies (cagPAI, abfA) help, but a clearer biological interpretation would strengthen the paper.

5. **Inversion fraction correlation** — Pearson r = 0.936 is strong, but the low overall correlation of `breakpoint_count` (partial r ~0.4) may be seen as a limitation. Explain that breakpoint_count is a count confounded by fragmentation, while inverted_fraction is the fragmentation-invariant ratio.

6. **Enzyme panel choice** — default BcgI/AlfI/AloI/FalI is stated but the optimization process is not detailed in the main text. The paper earlier explored many enzymes; reviewers may ask how the panel was chosen. Add a panel-selection paragraph or supplementary note.

7. **Method length** — Online Methods is very long. Nature Methods may ask to move detailed likelihood derivations, calibration protocol, and database-screen details to Supplementary Information.

## Recommended revision plan

### Phase A — fix submission mechanics (1–2 days)
- [ ] Remove/update placeholder lines (affiliation email, acknowledgements, author contributions, competing interests).
- [ ] Resolve Table 3 figure/number inconsistencies.
- [ ] Generate/copy missing supplementary figures S10, S11, S14, S15, S16 and link S17–S20 explicitly.
- [ ] Create `figures/submission/` or rename final figure files to match main-text figure numbers.
- [ ] Verify all references are cited and formatted correctly.

### Phase B — strengthen scientific narrative (2–3 days)
- [ ] Add a "Limitations" subsection in Discussion or expand the existing three-boundary paragraph.
- [ ] Add a paragraph on enzyme-panel selection and optimization to Results or Methods.
- [ ] Consider moving lengthy Methods derivations to Supplementary Notes.
- [ ] Add a sentence in Abstract highlighting the structural/SV output as a second axis beyond ANI.

### Phase C — technical validation (ongoing, depends on HPC)
- [ ] Complete FracMinHash scale sweep analysis when HPC jobs finish.
- [ ] Decide whether the closed-genome coordinate-level validation result should be included, downplayed, or replaced with count-level correlation.

## Files to check

- `paper/manuscript/manuscript_nature_methods.md`
- `paper/manuscript/SUPPLEMENTARY.md`
- `figures/report/` and `figures/syntracker_validation/`
- `case_studies/*/figures/`
- `results/closed_inversions/JUNCTION_COORDINATE_REPORT.md`
