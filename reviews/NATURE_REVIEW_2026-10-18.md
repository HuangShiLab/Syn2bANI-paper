# Nature-style review package — 18 October 2026

## Review process

Three independent mock reviews of
`paper/manuscript/manuscript_nature_methods.md` used the installed
`nature-reviewer` workflow. The reviews separately examined claims, methods,
statistics and reproducibility. This file records the consensus major items and
the revision response; minor wording items were handled during the same pass.

## Consensus blocking issues

### R1-M1 / R2-M1 — census tool attribution

**Issue.** The large censuses use the companion Syn2b structural channel and
skani for ANI. They are not an end-to-end Syn2bANI run and must not be
described as such.

**Response.** Census claims were reframed as an anchor-derived hybrid census:
the Syn2b structural channel was joined to skani 0.3.2 ANI by canonical
unordered genome pairs. Syn2bANI-specific claims are restricted to its own
benchmark and case-study reruns.

### R1-M2 / R2-M2 / R3-M1 — v5/v6 calibration provenance

**Issue.** The unified benchmark used calibration v6, but manuscript prose and
Table 1 called it v5. The deployed software snapshot contains v5.

**Response.** The manuscript now distinguishes deployed v5 from benchmark v6.
V6 was trained on the 2,520 v5 pairs plus 1,614 usable high-ANI training pairs
and was evaluated only on genomes excluded from both training sources. Table 1
identifies the 2,074-pair CV row as v5 and the 43,334-pair and unified rows as
v6. Provenance is archived in
`results/gtdb50k/CALIBRATION_PROVENANCE.md`.

### R1-M3 / R2-M3 — structural-validation and Figure 5 provenance

**Issue.** The v0.1.1 breakpoint rerun and v032 SynTracker rerun were not fully
archived, and the Figure 5 wording did not state the regenerated statistics.

**Response.** `results/gtdb50k/rerun_v011/`,
`results/gtdb50k/sv_reanalysis_v011.*`,
`results/syntracker_validation/rerun_v032/`, the regeneration script and the
regenerated Figure 5 are archived. The caption now gives r = 0.935 (n = 43,312)
for panel a, r = 0.996 (n = 3,826; slope = 0.986) for panel b, and partial
Pearson r = 0.411 (n = 43,334) for panel c.

### R1-M4 / R2-M4 / R3-M2 — FDA-ARGOS AMR source

**Issue.** The FDA-ARGOS *E. coli* resistance genotypes were described as CARD,
but they were derived by AMRFinderPlus through NCBI Pathogen Detection.

**Response.** Results, Methods, Table 3 and Table 4 now identify AMRFinderPlus
and NCBI Pathogen Detection. The separate *H. pylori* collaborator metadata,
which does use curated CARD entries, remains explicitly distinct.

### R1-M5 / R3-M3 — GTDB breakpoint denominator

**Issue.** Figure 11 primarily showed reported-pair fractions even though
10.08% of reported pairs lacked a numeric breakpoint.

**Response.** Figure 11b now displays paired bars for both the reported and
breakpoint-informative denominators. Results gives both values, and the caption
defines the informative denominator. Full sensitivity remains in
`results/census/gtdb/effective_quality_sensitivity.tsv`.

### R1-M6 / R2-M5 / R3-M4 — dependence-aware uncertainty

**Issue.** Benchmark and census pairs are not independent because genomes occur
in many pairs.

**Response.** New genome-level bootstrap intervals were added for the GTDB
43k benchmark. Species-level sensitivity intervals were added for OAPGC and
HROM. The report states explicitly that species-level resampling does not
remove within-species pair dependence, shared ancestry or quality correlation.
Outputs are in `results/revision/`.

### R1-M7 / R2-M6 / R3-M5 — HROM quality confounding

**Issue.** Structural calls can be more complete for high-quality assemblies,
and species labels should not be equated with measured mutation or
recombination rates.

**Response.** The manuscript retains the quality-stratified HROM analysis,
uses “candidate sequence-level” and “candidate architecture-level” labels, and
now includes dependence-aware species intervals. It explicitly avoids direct
hypermutator or hyper-recombinator claims for census-derived lists.

### R1-M8 / R2-M7 / R3-M6 — biological overclaiming

**Issue.** Treatment-outcome language for *B. longum* abfA, outbreak language
for FDA-ARGOS, “SNP” as a synonym for ANI-derived distance, and “within-site”
without defining its two habitats risk overinterpretation.

**Response.** Wording was narrowed: the abfA result supports a curated
locus-status association rather than outcome prediction; FDA-ARGOS *E. coli*
are very-high-ANI reference-genome pairs rather than inferred outbreak clones;
OAPGC uses ANI-derived sequence distance rather than base-level SNP counts; and
within-site is defined as within-habitat (oral-oral or airway-airway).

## Additional reconciliation

- HROM all-pair mode counts were reconciled to 96 structural-enriched, 160
  ANI-enriched, 36 both and 229 not enriched. The stale quality report and TSV
  were corrected.
- The GTDB ANI merge `header_ok: false` is explained as C-locale header
  placement in a sorted stream; all substantive merge checks pass
  (`results/census/gtdb/GTDB_ANI_MERGE_QC_NOTES.md`).
- The *E. coli* GTDB shortfall (1,622,961 missing pairs; 99.5210% coverage) and
  the contiguity sensitivity at N50 ≥100 kb remain explicit limitations.

## Post-revision verification

- Reference 41 is now complete as Song, X. et al., *Proc. Natl. Acad. Sci.
  USA* **123**, e2603926123 (2026); reference 44 is complete. The reference
  list contains 46 consecutive entries.
- Figure 3 and the regenerated Supplementary Figures S2 and S4 identify the
  benchmark model as v6; the separate 2,074-pair CV result remains explicitly
  v5.
- Figure 10 was regenerated with ANI-derived sequence-distance and
  within-habitat terminology. Figure 11 caption now explicitly describes the
  paired reported/informative denominators.
- A regenerated Supplementary Figure S15 was added from the FDA-ARGOS
  pair-level PAF, AMR and validation tables; its source counts are 41/16/46
  structural classes and 60/11/3 AMR-discordant pairs at the three ANI
  thresholds.
- Supplementary Table S7 was added for FDA-ARGOS structural and AMR results.
  The submission package now contains Supplementary Figures S1–S15, Tables
  S1–S7 and Notes 1–7. Both Markdown and DOCX copies were regenerated.
- Data Availability now links the calibration-provenance, dependence-aware
  uncertainty and GTDB merge-QC outputs, while continuing to state that the
  pair-level census tables are planned for Zenodo.
- Final checks passed: `git diff --check`; Python syntax checks for edited
  plotting/analysis scripts; cross-reference audit for all Supplementary
  figures, tables and notes; a 46-reference audit; and numerical spot checks
  of the GTDB, HROM and OAPGC headline results against their source tables.
