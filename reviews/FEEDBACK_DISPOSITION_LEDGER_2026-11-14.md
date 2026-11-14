# Feedback disposition ledger: Syn2bANI, 2026-11-14

This ledger consolidates the synthesized 2026-10-21 review round and the two
subsequent skill audits. The 2026-10-21 reviewer files remain frozen; this file
is the current disposition record.

| ID | Source | Location | Request / failure mode | Disposition | Verification |
|---|---|---|---|---|---|
| C01 | Reviewer synthesis R1-M1 | Results, Table 1, hybrid analysis | Do not validate the benchmark-selected 98% threshold as independent | fixed | Results/Table 1 label hybrid as operational sensitivity; threshold selection stated in Methods and `results/gtdb50k/CALIBRATION_PROVENANCE.md` |
| C02 | Reviewer synthesis R2-M1 | OAPGC interpretation | Do not equate `100 − skani ANI` with SNP load or directional migration | fixed | Methods/Results use ANI-derived distance; recent-strain power limitation retained; alias table is machine readable |
| C03 | Reviewer synthesis R2-M2 | Structural and HROM interpretation | Do not treat breakpoint counts as a complete SV census or labels as measured hypermutator states | fixed | Results call breakpoint count a detection statistic; HROM labels are candidate anchor-statistic phenotypes |
| C04 | Reviewer synthesis R2-M3 | GTDB census scope | Resolve cluster scope, missing-breakpoint denominator and dominance | fixed | `results/revision/gtdb_cluster_scope_sensitivity.tsv`; Results distinguishes 1,311 complete species and 1,316 all eligible clusters; denominators shown |
| C05 | Reviewer synthesis R3-M1 | Provenance | Pin software, models, scripts and artifacts | fixed | `REPRODUCIBILITY_MANIFEST.md`; `scripts/audit_release.py` verifies checksums |
| C06 | Reviewer synthesis R3-M4 | Data availability | Separate immediate review set from full archive plan | fixed | Data availability lists review outputs and immutable Zenodo release condition |
| C07 | Reviewer synthesis R3-M5 | Structural validation | Validate structural metrics separately | fixed | Results/Fig. 4 state inversion-fraction, breakpoint and coverage-specific evidence and limits |
| C08 | Skill audit 2026-11-12 | Fig. 2 comparator | Use skani 0.3.2 consistently | fixed | Source and rendered MAE 0.477; manuscript updated; `scripts/audit_release.py` passes |
| C09 | Skill audit 2026-11-12 | Fig. 4 regression | Avoid mixing plotted and quoted regression directions | fixed | Caption states displayed OLS direction and reverse-regression slopes |
| C10 | Skill audit 2026-11-14 | Simulation manifest | Remove host-local absolute paths | fixed | `results/simulations/simindel/manifest.tsv` uses source-repository-relative paths; checksum updated |
| C11 | Skill audit 2026-11-14 | Fragmentation denominator | Correct 20–201 to current 20–200 contig levels | fixed | Results, Methods, caption and Table 3 now say 20–200; audit gate checks absence of 20–201 |
| C12 | Skill audit 2026-11-14 | Calibration wording | Avoid implying all complete-genome calibration is unavailable | fixed | Abstract says benchmark-only v6 is not embedded; v5 remains deployed in Methods |
| C13 | Skill audit 2026-11-14 | abfA denominator and claim | Give locus-analysis n and avoid stronger deletion claim than curated labels | fixed | Fig. 5 states 183 non-reference genomes, 127/56; title uses cluster difference, not loss |
| C14 | Release audit 2026-11-14 | Code/census archives | Provide durable archive or venue-compatible release condition | fixed for submission | Code Zenodo DOI and immutable census archive are explicitly promised before acceptance; GitHub commits pin review snapshot |

No actionable feedback row is `deferred` or `unresolved`. The acceptance-stage
Zenodo actions are release conditions, not unresolved manuscript-text questions.

## Individual minor-comment dispositions

The following rows are keyed to the frozen 2026-10-21 reviewer files. The
corresponding substantive revisions are in the current manuscript, Supplement,
source tables or audit script; repository history contains the frozen reports.

| ID | Source | Request / issue | Disposition | Verification |
|---|---|---|---|---|
| R1-m1 | Reviewer 1 | Use finite-row v5 accounting and state the 21 omitted source pairs | fixed | Supplementary Note 6 says 2,053 finite + 467 = 2,520 and states 21 lacked finite features; `results/gtdb50k/CALIBRATION_PROVENANCE.md` agrees |
| R1-m2 | Reviewer 1 | Clarify v5/v6 roles in Discussion | fixed | Discussion states deployed versus benchmark-only calibration and raw gated recommendation for MAGs/flagged pairs |
| R1-m3 | Reviewer 1 | State Figure 3 source and complete-case denominators | fixed | Fig. 3 legend gives 43,334+727 source pairs and 40,629 plotted complete-case pairs; non-call handling stated |
| R1-m4 | Reviewer 1 | Explain FastANI 363-pair subset and non-call policy | fixed | Table 1 caption states 363 available-call pairs and non-calls were not imputed |
| R2-m1 | Reviewer 2 | Distinguish estimator hybrid rule from two-layer census | fixed | Results/Methods use hybrid estimator and two-layer census terminology; census is explicitly not end-to-end Syn2bANI |
| R2-m2 | Reviewer 2 | Make GTDB Fig. 11a/11b scope difference prominent | fixed | Results and ED Fig. 5 distinguish 1,311 complete species from 1,316 eligible clusters and both breakpoint denominators |
| R2-m3 | Reviewer 2 | Avoid calling HROM labels divergence modes | fixed | Results/ED Fig. 3 use candidate anchor-statistic phenotypes, not measured mutation/recombination modes |
| R2-m4 | Reviewer 2 | Rename `SNP_auc` or document ANI-derived meaning | fixed | Methods documents `100 − skani ANI`; released alias is `species_effects_all_pairs_ani_distance.tsv` |
| R2-m5 | Reviewer 2 | Keep locality/resistance analyses exploratory and lineage aware | fixed | Methods label locality analyses exploratory post hoc; Results state macrolide association did not survive lineage stratification |
| R2-m6 | Reviewer 2 | Clarify AMR genotype rather than phenotype | fixed | Results, Methods and Supplement state genotype not phenotype was tested |
| R2-m7 | Reviewer 2 | Give companion Syn2b provenance/parameters/format | fixed for submission | Methods pins Syn2b snapshot `64717ab`, four enzymes, `.tgt` anchors and census block/pair construction; exact command now added in Methods |
| R3-m1 | Reviewer 3 | Use ANI points consistently | fixed | Abstract, Results, Table 1 and captions use ANI points/ANI-derived distance |
| R3-m2 | Reviewer 3 | Audit dataset-to-figure mapping | fixed | Extended Data Table 1 mappings checked by release audit and updated after display renumbering |
| R3-m3 | Reviewer 3 | Explain 484 runtime denominator | fixed | Table 2 caption states ordered n² comparisons of 22 genomes including self |
| R3-m4 | Reviewer 3 | Replace vague S. aureus country statement | fixed | Methods/captions state 351/2,211 pairs have no country on either side |
| R3-m5 | Reviewer 3 | Refine S4 quality title/caption | fixed | Supplement title says little association but extreme fragmentation has more error |
| R3-m6 | Reviewer 3 | List all eight cagPAI pilot genomes | fixed | Methods lists wild type, 26695, mut1, three structural events, and mut1+event backgrounds |
| R3-m7 | Reviewer 3 | Separate prespecified versus exploratory locality tests | fixed | Fig. 6 caption and Methods state primary associations and exploratory locality analysis separately |
| R3-m8 | Reviewer 3 | Clarify hardware/run configuration | fixed for submission | Table 2 caption reports hardware/core context and Methods states 32-thread database run; run-size row is explicit |
| R3-m9 | Reviewer 3 | Avoid misleading 50k shorthand | fixed | Results uses 43,334; the historical path name remains only for artifact compatibility |
| R3-m10 | Reviewer 3 | Explain blank unified high-ANI cells | fixed | Table 1 caption says source 95–100% pairs are split at 97% and summarized in all column |

No individual minor row remains deferred or unresolved.

## Current release-audit dispositions

| ID | Source | Request / issue | Disposition | Verification |
|---|---|---|---|---|
| V01 | Supervisor review 2026-11-14 | Refresh external ledger and delivery gate | fixed | This ledger and `reviews/SUPERVISOR_LOOP_2026-11-14.md` use Supervisor-Skills HEAD `2edb32733e4042a0408824441f0e832e2476b07b` and the current 55-check gate |
| V02/V03 | Supervisor review 2026-11-14 | GC sweep lacked reproducible inputs and source access | fixed by scope change | Unsupported GC ladder removed from Methods, Supplement, ED Table 1 and release scripts |
| V04 | Supervisor review 2026-11-14 | ED Fig. 2 lacked generation provenance | fixed | `scripts/generate_extended_data_figure2.py` regenerates packaged PDF/PNG from `paper/figures/main/fig6_syntracker_structure.png` |
| V06 | Supervisor review 2026-11-14 | Superseded artifacts remained in release repository | fixed | Old manuscript/archive/report files removed from current release branch; history remains in Git |
| V07/V08 | Supervisor review 2026-11-14 | Remove process history and unverified test-suite claim | fixed | Main text/caption process language removed; unsupported 120-test statement removed |
