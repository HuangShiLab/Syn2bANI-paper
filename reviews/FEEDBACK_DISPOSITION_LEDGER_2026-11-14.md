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
