# Syn2bANI skill audit, 2026-11-12

## Scope

Applied the four installed skills derived from the Strain2bScan retrospective:

1. `benchmark-design-guard`;
2. `paper-architecture`;
3. `claim-evidence-check`;
4. `supervisor-review`.

The Strain2bScan-specific historical blockers were not treated as Syn2bANI
facts. Their generalized failure modes were mapped to this project: comparator
version drift, clean-checkout reproducibility, claim-to-source agreement,
figure regression direction, denominator transparency, display architecture and
pre-delivery review.

## Initial findings before repair

| Severity | Finding | Evidence |
|---|---|---|
| CRITICAL | Abstract, Results and Figure 2 caption quoted skani MAE 0.377, but Methods stated skani 0.3.2 and the rendered Figure 2 reported 0.477. | `results/sv_validation/SKANI_0.3.2_RERUN.md`; `paper/figures/main/final/figure2.pdf` |
| CRITICAL | A clean checkout could fail Figure 2 generation because the simulation script fell back to a host-local prototype directory for required TSVs. | `scripts/generate_fig2_simulations.py`; absence of `results/simulations/` |
| MAJOR | No persisted `PROJECT_FRAME.md`, `PAPER_ARCHITECTURE.md` or `CLAIM_EVIDENCE_MATRIX.md` existed for the new skills. | repository audit |
| MAJOR | Figure 4 quoted the reverse-regression slopes without saying that the plotted ordinary least-squares lines use the opposite x/y orientation. | rendered slopes 0.872 and 0.986; source-report reverse slopes 1.004 and 1.006 |
| MAJOR | The packaged S. aureus Supplement Fig. S10 was current, but its generator would have recreated the obsolete pre-v0.1.1 triangle and embedded an old binary commit in the title. | old generator used `triangle.tsv` and `/Volumes/MoneyCat`; packaged PDF used v0.1.1 rerun |
| MINOR | Older internal reports contain historical skani 0.1.0 or duplicated fragmentation-row values; they were not labeled as superseded at the exact claim-source level. | `paper/data/SIMULATION_AND_PERFORMANCE_REPORT.md`; `results/SIMULATION_RERUN_REPORT.md` |

## Repairs made

1. Copied the self-contained simulation outputs and manifests to
   `results/simulations/`.
2. Regenerated Figure 2 from `results/simulations/` and copied the PDF into
   `paper/figures/main/final/figure2.pdf`.
3. Updated Abstract, Results, Figure 2 caption and Data availability to the
   skani 0.3.2 MAE of 0.477 and the self-contained input paths.
4. Clarified Figure 4 regression orientation in Results, caption and
   Supplement. The displayed fits are dnadiff-on-Syn2bANI; the quoted 1.004 and
   1.006 slopes are Syn2bANI-on-dnadiff reverse regressions.
5. Added `PROJECT_FRAME.md`, `PAPER_ARCHITECTURE.md`,
   `CLAIM_EVIDENCE_MATRIX.md` and `results/simulations/CLAIM_VALUE_AUDIT.md`.
6. Replaced the stale S. aureus supplement-figure generator with the
   self-contained v0.1.1 rerun and regenerated Fig. S10.
7. Packaged mosaic simulation outputs and removed host-local fallbacks from
   Fig. 2 and Fig. S6 generators.
8. Added the simulation artifacts, checksums and result mapping to
   `REPRODUCIBILITY_MANIFEST.md`.

## Post-repair checks

- Main figures: 6; main tables: 2.
- Extended Data figures: 5; Extended Data tables: 2.
- References: 46 consecutive entries.
- Manuscript/supplement repo-path audit: no missing paths.
- Extracted Figure 2 labels: Syn2bANI 0.073, skani 0.477, FastANI 0.742.
- Exact simulation source recalculation: 0.0732, 0.4775, 0.7423.
- Em dashes: 0 in manuscript and Supplement.
- Relevant Fig. 2, Fig. S6 and Fig. S10 plotting scripts compile; current display generators no longer require host-local fallbacks.
- Regenerated DOCX files are current.
- Repository `git diff --check`: clean after regeneration.

## Supervisor-review verdict

**Ready for submission after the 2026-11-12 repair.** The original 2026-10-27
verdict is superseded by this audit because the earlier pass did not catch the
skani-version claim mismatch or the host-local figure dependency. The remaining
non-blocking release action is unchanged: deposit complete pair-level census
tables and checksum/accession manifests on Zenodo before acceptance.

## Severity accounting

| Severity | Before skill audit | After repair |
|---|---:|---:|
| CRITICAL | 2 | 0 |
| MAJOR | 3 | 0 |
| MINOR | 1 | 1 retained historical-note issue only |

The remaining MINOR is not a manuscript blocker: older internal technical
reports remain in the repository as historical records. The current manuscript,
`results/simulations/CLAIM_VALUE_AUDIT.md`, packaged figures and
`REPRODUCIBILITY_MANIFEST.md` identify the current source tier and values.

## Final package

- Manuscript DOCX regenerated at 2026-11-12.
- Supplement DOCXs regenerated and synchronized.
- Main display package unchanged except corrected Figure 2.
- Supplementary Figures S6 and S10 regenerated from self-contained or repository data.
- No manuscript path is missing.
