# Revision plan — Syn2bANI Nature Methods submission

Response to the mock peer review in `paper/MOCK_REVIEW.md`. This plan incorporates the post-review finding that the *H. pylori* cagPAI `complete_rearranged` state is largely a circular-origin coordinate artifact. Items are tagged by priority.

---

## Strategic repositioning of the four structural columns

Based on the GTDB-R207 SV reanalysis and the fragmentation theorem in the companion Syn2b repository, the four structural outputs should be repositioned as follows:

| Metric | New role | Evidence |
|---|---|---|
| `raw_inverted_fraction` | **Headline structural output** | Pearson r = 0.936 vs dnadiff across 43k pairs; r = 0.996 at ANIm ≥ 97%; fragmentation-invariant by construction. |
| `breakpoint_count` | Secondary count-based descriptor | Partial r = 0.414 vs `dnadiff breakpoints`, but rises to ~0.79 vs `dnadiff relocations + inversions`; must be reported with contig control. |
| `synteny_blocks` | Assembly-structure diagnostic | 62% of variance tracks contig count; not a rearrangement metric. |
| `af_query` / `af_reference` | Coverage / aligned-fraction metrics | Correlate with dnadiff/minimap2 coverage; not synteny. |
| `anchor_adjacency` | Local anchor-order conservation | Separate combinatorial metric; not a proxy for base-pair synteny. |

The manuscript's structural section should be rewritten around this hierarchy.

---

## Phase 1 — Blockers (must complete before manuscript rewrite)

### 1. Submit the three high-ANI SLURM arrays

**What:** `s7_s2b_high_ani_final.slurm`, `s8_skani_high_ani_final.slurm`, `s9_fastani_high_ani_final.slurm` in `results/gtdb50k/`.

**Why:** The 727 high-ANI test set is underpowered for the ≥97% regime. The unified `high_ani_pairs_final.tsv` has 2,348 pairs (95–97%: 322; 97–100%: 2,026), which is needed for v6 calibration and for validating `raw_inverted_fraction` in the strain-level range.

**Steps on HPC:**
```bash
cd /lustre1/g/aos_shihuang/Syn2bANI-paper
git pull origin main
cd results/gtdb50k
mkdir -p logs
sbatch ../../scripts/gtdb50k/s7_s2b_high_ani_final.slurm
sbatch ../../scripts/gtdb50k/s8_skani_high_ani_final.slurm
sbatch ../../scripts/gtdb50k/s9_fastani_high_ani_final.slurm
```
After completion:
```bash
python3 scripts/gtdb50k/merge_high_ani_results.py
python3 scripts/calibration_v6.py
```
Pull `high_ani_results.tsv`, `linear_cal_v6.json`, and `calibration_v6_*.tsv` back to the local repo.

**Priority:** Critical blocker.  
**Effort:** Minutes to submit; hours to run on HPC.  
**Owner:** HPC execution; manuscript updates after pull.

### 2. Fix circular-origin handling in Syn2bANI `struct`

**What:** Add per-contig circularity, circular start-coordinate normalization, and a genome-spanning SV filter.

**Why:** The *H. pylori* cagPAI case study is 95% artifact without this. The closed-genome cohort's 36% median inversion rate is likely the same issue.

**Implementation notes:**
- Syn2b already has per-contig circularity (`0c4c541`). Syn2bANI's `struct` does not.
- Before pairwise comparison, rotate circular contigs to a common start (e.g., smallest k-mer / dnaA / fixed coordinate).
- After calling, flag any SV whose span exceeds a threshold fraction of the contig length (default 50%) as `COORDINATE_ARTIFACT` and exclude it from rearrangement counts.
- The companion Syn2b repository's fix should be ported or shared via the same logic.

**Priority:** Critical blocker for cagPAI case study.  
**Effort:** Small to medium (Rust code change + recompile + re-run).  
**Owner:** Syn2bANI Rust code; `case_studies/h_pylori_cagpai/` reanalysis.

### 3. Diagnose closed-genome `raw_inverted_fraction` distribution and 0/371 seed bug

**What:** Two independent problems in the closed-genome all-vs-all output (`results/gtdb50k/syn2b_inverted_fraction_closed.tsv`):
- Median inversion fraction 0.36 is biologically implausible for same-species pairs.
- 371 expected seed pairs appear 0 times in 61,537 outputs.

**Why:** `raw_inverted_fraction` is the proposed headline metric; it cannot be presented until these are explained.

**Hypotheses to test:**
- Median 36% is circular-origin / global orientation artifact (same as *H. pylori* but on many species).
- 0/371 seed pairs is an identifier mismatch or a filtering bug in the runner.

**Priority:** Critical for headline metric credibility.  
**Effort:** Small (diagnostic scripts on existing tables).  
**Owner:** `results/gtdb50k/CLOSED_GENOME_INVERSION_REPORT.md` update; new diagnostic script.

---

## Phase 2 — Manuscript rewrite

### 4. Restructure the structural outputs section

**Proposed revision:**
- Lead with `raw_inverted_fraction` as the primary length-weighted structural metric.
- Present `breakpoint_count` as a count-based descriptor, always with contig-count control and validated against `dnadiff relocations + inversions` (not `dnadiff breakpoints`).
- Demote `synteny_blocks` to an assembly-structure diagnostic.
- Demote `af_query` to a coverage metric.
- Keep `anchor_adjacency` as a local anchor-order metric with the existing careful caveats.
- Add a table or flowchart mapping each metric to its biological question.

**Priority:** Essential.  
**Effort:** Medium (rewriting Results 2.8 and Discussion).  
**Owner:** Manuscript.

### 5. Recompute or remove the *H. pylori* cagPAI rearrangement associations

**Proposed revision:**
- Re-run `syn2bani struct` vs 26695 with circular-origin normalization.
- Apply the genome-spanning SV filter and reclassify `complete_rearranged`.
- Recompute χ² associations for disease stage **with FastBAPS lineage as a stratifying variable** (CMH or logistic regression) — the raw disease association may be entirely lineage-driven.
- If the rearranged state collapses to a small number of true local events, reframe the result as: "cagPAI presence/absence is associated with lineage and disease stage; large-scale rearrangement calls were dominated by coordinate artifacts."

**Priority:** Essential.  
**Effort:** Small to medium (depends on re-running struct).  
**Owner:** `case_studies/h_pylori_cagpai/`.

### 6. Sharpen the innovation claim

**Proposed revision:**
- Frame the advance as: coordinate-carrying anchors enable ANI + structural descriptors in one pass, whereas k-mer tools regress chain statistics into a scalar and discard geometry.
- Update Figure 1 with a contrast panel.
- Add a paragraph in the Discussion comparing Syn2bANI with skani's chain-statistic approach.

**Priority:** Essential.  
**Effort:** Writing + one figure panel.  
**Owner:** Manuscript (`Introduction`, `Discussion`, Figure 1).

### 7. Make output selection transparent

**Proposed revision:**
- Add the hybrid threshold sweep to Supplementary Note 6.
- Add a user-facing decision flowchart.
- Move the raw/calibrated/hybrid recommendation to the Abstract and Discussion.
- Update CLI help/README to state `--calibrate` is for complete genomes only.

**Priority:** Essential.  
**Effort:** Low.  
**Owner:** Manuscript, Supplementary Note 6, README.

---

## Phase 3 — Validation strengthening (if time permits)

### 8. Validate structural outputs against a dedicated SV tool

**Proposed revision:**
- Run SyRI or minigraph-Cactus on the Enterobacteriaceae completes and/or top Syntracker discordant pairs.
- Report per-event precision/recall for inversions, translocations, and indels >1 kb.
- If this is not feasible, reframe structural outputs as "fast synteny/structural descriptors" rather than "SV calls" and add a clear limitations sentence.

**Priority:** Recommended if "call" is retained; optional if reframed as descriptors.  
**Effort:** High (new HPC runs + parsing).  
**Owner:** New results directory.

### 9. Runtime decomposition and screen recall

**Proposed revision:**
- Add a Supplementary table decomposing runtime by stage.
- Report absolute numbers for Stage 1 screen recall on GTDB-R207.

**Priority:** Recommended.  
**Effort:** Low (re-aggregate logs).  
**Owner:** Supplementary Information.

---

## Phase 4 — Polish

### 10. Minor revisions

| # | Item | Priority | Action |
|---|---|---|---|
| 10.1 | Reserve "ANI" for default output | Recommended | Rename ambiguous table rows/text. |
| 10.2 | Figure clarity | Recommended | Update Figure 1 schematic; enlarge Figure 6 annotations. |
| 10.3 | Supplementary completeness | Essential | Fix missing S10–S13 captions; copy Syn2b math review as Supplementary Note. |
| 10.4 | Reproducibility | Recommended | Add `Dockerfile` or `environment.yml`. |
| 10.5 | Citations | Recommended | Spot-check Nature Methods reference style. |

---

## Execution order

1. Submit high-ANI arrays (Phase 1.1).
2. In parallel: fix circular-origin handling (Phase 1.2) and diagnose closed-genome bugs (Phase 1.3).
3. Once 1.1 returns: train v6 calibration and validate in the ≥97% range.
4. Once 1.2 is fixed: re-run *H. pylori* cagPAI and recompute associations (Phase 2.5).
5. Rewrite structural section (Phase 2.4), innovation claim (Phase 2.6), and output selection (Phase 2.7).
6. Optional: dedicated SV tool validation (Phase 3.8) and runtime decomposition (Phase 3.9).
7. Polish and final citation/style check (Phase 4).

---

## Fallbacks

- **If high-ANI arrays cannot be submitted soon:** The manuscript can still proceed with the existing 727-pair high-ANI set, but the ≥97% validation will remain weak.
- **If circular-origin fix is delayed:** Remove the cagPAI rearrangement associations and use only the presence/absence (`empty`/`partial`/`complete`) axis, which is unaffected.
- **If dedicated SV validation is infeasible:** Reframe structural outputs as "fast synteny descriptors" rather than "SV calls."
