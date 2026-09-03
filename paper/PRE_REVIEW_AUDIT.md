# Pre-review audit — Syn2bANI manuscript

## Manuscript status

- **Format**: Nature Methods markdown draft (`paper/manuscript/manuscript_nature_methods.md`)
- **Word count**: Abstract ~153 words; main text ~4,500 words (within target).
- **Sections**: Abstract, Introduction, Results (8 subsections), Discussion, Online Methods, Data/Code availability, References, Figure Legends, Tables.
- **Placeholders remaining**: Acknowledgements, Author contributions, Competing interests.

## Structural completeness

| Section | Status | Notes |
|---|---|---|
| Title / authors / affiliations | OK | Correspondence email updated to shihuang@hku.hk throughout. |
| Abstract | OK | Could mention structural outputs and case studies more explicitly. |
| Introduction | OK | Clear motivation and differentiation from k-mer tools. |
| Results | OK | 8 subsections; logic flows from estimator → simulations → real benchmarks → structural outputs → case studies → efficiency. |
| Discussion | OK | Three boundaries are well stated; could add a dedicated "Limitations" subheading. |
| Methods | OK | Very detailed; may be too long for Nature Methods (consider moving some derivations to Supplementary). |
| Data/Code availability | OK | Repo links present; commit `fe0f36c` should be verified as current. |
| References | OK | Ref 10 updated to 2025. Ref 16 now cited in Introduction. |

## Critical issues

### 1. Figure file naming vs. figure numbers — FIXED

Created `figures/submission/` with source PNGs mapped to main-text figure numbers and a `README.md` assembly guide.

### 2. Missing supplementary figures — FIXED

| Figure | Status | Notes |
|---|---|---|
| S10/S11 | REMOVED | Main-text Fig. 6 now uses the existing four-panel `syntracker_ani_vs_breakpoints.png`; orphaned S10/S11 sections were removed from Supplementary. |
| S14 | OK | Linked to existing `figures/report/fig_inverted_fraction_comparison_high_ani.png`. |
| S15 | OK | Copied from `case_studies/ecoli_o157_fitzgerald_2021/figures/`. |
| S16 | OK | Copied from `case_studies/fda_argos_s_aureus/figures/`. |
| S17–S20 | OK | Copied from `case_studies/h_pylori_cagpai/results/struct_extended/` and assigned explicit numbers. |

### 3. Table 3 inconsistencies — FIXED

- Row 10: changed "Fig. S3" → "Supplementary analysis".
- Row 11: changed 2,342 → 727 pairs to match text; added "non-representative" and exclusion note.
- Row 15: removed erroneous "Fig. 5" reference.
- Fig. 6 caption and main-text references updated to use existing breakpoints figure; 45,000 FastANI comparison no longer claims to be in Fig. 6.

### 4. Correspondence email inconsistency — FIXED

Updated line 267 from `huangshi@njau.edu.cn` to `shihuang@hku.hk`.

### 5. Reference 10 citation — FIXED

Updated to: Enav, H., Paz, I. & Ley, R. E. *Nat. Biotechnol.* **43**, 773–783 (2025).

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
- [x] Remove/update placeholder lines (affiliation email).
- [x] Resolve Table 3 figure/number inconsistencies.
- [x] Generate/copy missing supplementary figures S14–S20; remove orphaned S10/S11.
- [x] Create `figures/submission/` with main-text figure files and assembly README.
- [x] Update placeholders: Author contributions, Competing interests (Acknowledgements left for funding details).
- [x] Verify references: Ref 10 updated to 2025 Nat. Biotechnol. 43:773–783; Ref 16 (dRep) now cited in Introduction.
- [x] FracMinHash scale sweep completed and reported in `results/gtdb50k/FRACMINHASH_VALIDATION.md`.

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
