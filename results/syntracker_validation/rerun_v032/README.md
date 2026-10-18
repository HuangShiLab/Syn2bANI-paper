# skani 0.3.2 rerun of the SynTracker cohorts — drift check

Rerun 2026-09-17 (HPC job 4071698; `skani 0.3.2`, note: it reports a
"Learned ANI mode" regression). Files: `skani_<cohort>.tsv` in this
directory (unique pairs, C(n,2) per cohort). The 0.1.0 files in
`data/syntracker_validation/skani/` are retained for figure reproducibility;
their raw triangles also contain cross-contig rows (min values down to 1.00)
that the figure pipeline filters.

| cohort | pairs | 0.3.2 min / median / max | manuscript claim | verdict |
|---|---:|---|---|---|
| E. coli hypermutator | 253 | 99.93 / 99.99 / 100.00 | 99.93–100.00, median 99.99 | identical |
| H. pylori | 2,926 | 94.62 / 95.68 / 100.00 | 94.6–100% | identical |
| N. gonorrhoeae | 66 | 99.61 / 99.93 / 100.00 | 99.6–100% | identical |
| S. rimosus | 190 | 99.91 / 99.97 / 100.00 | >99.9% | identical |

Every ANI value quoted in Results section 7 / Figure 4 reproduces exactly
under skani 0.3.2 in the >=94.6% regime, as the Methods note anticipated.
No figure regeneration required.

# syn2bani v0.1.1 structural channel (T4) — dual-implementation check

Rerun 2026-09-17 (HPC job 4074367; `syn2bani triangle --edge-list --verbose`,
per-pair `breakpoint_count`). Files: `syn2bani_el_<cohort>.tsv`.

| cohort | pairs | syn2bani median (max) | Syn2b-channel median | % zero |
|---|---:|---|---:|---:|
| E. coli hypermutator | 253 | 0 (2) | 0 | 74% |
| H. pylori | 2,926 | 5 (12) | 7 | 28% |
| N. gonorrhoeae | 66 | 1 (4) | 3 | 32% |
| S. rimosus | 190 | 3 (9) | 10 | 8% |

Within-host H. pylori (n = 476): median 0 (max 2); between-host (n = 2,450):
median 6 — the Syn2b channel reports 0 vs 8. The epidemiological signal
reproduces through the second, independent implementation; the chain-based
counter reads systematically lower than Syn2b's adjacency counter (consistent
with the 43,312-pair GTDB concordance, Spearman 0.86, median ratio 0.30),
which is a scale difference in strict event evidence, not a ranking failure.
