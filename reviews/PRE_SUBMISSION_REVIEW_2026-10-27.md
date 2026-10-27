# Pre-submission review: Syn2bANI (Nature Methods format audit)

**Date:** 2026-10-27  
**Audited files:** `paper/manuscript/manuscript_nature_methods.md`, `paper/manuscript/SUPPLEMENTARY.md`, packaged figures/tables  
**Scope:** final editorial and submission-readiness audit after the 2026-10-21 scientific review. Frozen reviewer reports were not modified.

## Executive verdict

**Recommendation: ready for submission after this pass.** The manuscript now fits a compact Nature Methods Article format, uses one consistent display-item scheme, removes remaining editorial blockers, and preserves the scientific caveats introduced in earlier review rounds.

## Before/after severity summary

| Severity | Before | After |
|---|---:|---:|
| CRITICAL | 2 | 0 |
| MAJOR | 6 | 0 |
| MINOR | 4 | 0 |

The two pre-pass critical issues were a main body of ~6,874 words versus the ~3,000-word Article target and 11 main figures plus 4 main tables. The revised main body is ~2,240 words excluding Abstract, with 6 main figures and 2 main tables. Five former main figures are now Extended Data Figures 1–5, and former Tables 3–4 are Extended Data Tables 1–2.

## CRITICAL issues found before the pass

1. **Length exceeded the target Article format.**  
   *Evidence:* Abstract 203 words; Introduction 540; Results 5,891; Discussion 443; main body excluding Methods ~6,874.  
   *Fix:* condensed prose while retaining key effect sizes, complete-case denominators, calibration boundaries, quality stratification and census caveats; moved full result-level detail to Supplementary Note 8.  
   *Status:* fixed.

2. **Too many main display items.**  
   *Evidence:* Figures 1–11 and Tables 1–4 were all labeled as main displays.  
   *Fix:* retained main Figures 1–6 and Tables 1–2; reassigned computational, collection-structure and catalog-census displays to Extended Data Figures 1–5; reassigned dataset and case-study tables to Extended Data Tables 1–2.  
   *Status:* fixed.

## MAJOR issues found before the pass

1. **Tool-name inconsistency.** The title used `Syn2b-ANI`; the method is `Syn2bANI`. Fixed in the main title and supplement title.
2. **Prose em dashes.** Manuscript had 69 and supplement 51. All are now removed without changing data claims.
3. **Internal drafting note.** The obsolete status line containing manuscript word counts appeared before References. Removed.
4. **Incomplete source paths.** Shorthand paths for AMR, cagPAI locality, HROM enrichment, OAPGC effects, FDA-ARGOS validation, and case-study summaries were expanded to repository paths. The repository-path audit reports no missing manuscript paths.
5. **Figure 6 vector/label issue.** The source existed only as an unlabeled PNG. A labeled RGB-backed PDF is now packaged as Extended Data Figure 2 with panels a–d.
6. **Overlapping panel labels in cagPAI figure.** Regenerated the cagPAI figure, moved panel labels clear of titles, refreshed the PDF and packaged it as main Figure 6.

## Final checks

- Main body: Abstract 169 words; Introduction 266; Results 1,623; Discussion 301; total excluding Abstract ~2,190 words.
- Main displays: 6 figures and 2 tables.
- Extended Data: 5 figures and 2 tables.
- References: 46 consecutive entries, 1–46, with no duplicates.
- Cited main figures: 1–6; cited Extended Data figures: 1–5; no stale main references to old Figures 7–11.
- Repository paths cited in the manuscript/supplement: all audited paths exist.
- Manuscript and supplement em dashes: 0.
- Python files changed or used for regeneration compile successfully.
- `git diff --check`: passed at the time of review.

## Remaining submission-stage action

The aggregate tables, scripts, checksums and reduced review verification set are available now. The complete pair-level census tables and immutable Zenodo accession/checksum manifest should be deposited before acceptance, as already stated in Data availability. This is a data-deposit step, not a manuscript-text blocker.

## Recommendation

Proceed with submission. Use `paper/manuscript/manuscript_nature_methods.docx`, `paper/manuscript/SUPPLEMENTARY.docx`, `paper/figures/main/final/figure1.pdf` through `figure6.pdf`, and `paper/figures/extended_data/extended_data_figure1.pdf` through `extended_data_figure5.pdf`.
