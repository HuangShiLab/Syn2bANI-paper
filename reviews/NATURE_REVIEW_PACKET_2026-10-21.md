# Nature Methods review packet, 21 October 2026

## Scope

This packet is for three mutually blind Nature-style reviews of the current
Syn2bANI manuscript. Each reviewer receives the same source packet and
emphasis brief. Reviews must be grounded only in the files listed here.

## Primary manuscript

- `paper/manuscript/manuscript_nature_methods.md`
- `paper/manuscript/SUPPLEMENTARY.md`

## Source anchors

Use these files to verify numerical and reproducibility claims. Do not infer
additional experiments or data.

- Benchmark provenance: `results/gtdb50k/CALIBRATION_PROVENANCE.md`
- Benchmark reports: `results/gtdb50k/GTDB50K_HELDOUT_REPORT.md`, `results/gtdb50k/CALIBRATED_COMPARATORS.md`
- Uncertainty audit: `results/revision/DEPENDENCE_AWARE_UNCERTAINTY_REPORT.md`, `results/revision/dependence_aware_uncertainty.tsv`
- GTDB census: `results/census/gtdb/gtdb_species_ani_summary.tsv`, `results/census/gtdb/gtdb_species_ani_sv_joint.tsv`, `results/census/gtdb/GTDB_ANI_MERGE_QC_NOTES.md`
- GTDB sensitivity: `results/census/gtdb/effective_quality_sensitivity.tsv`, `results/census/gtdb/effective_species_summary.tsv`
- HROM census and quality analysis: `results/census/hrom/quality_stratified/QUALITY_STRATIFIED_REPORT.md`, `results/census/hrom/quality_stratified/quality_stratified_summary.tsv`
- OAPGC census and site tests: `results/census/oapgc/site_effects/species_effects_all_pairs.tsv`, `results/census/oapgc/site_effects/overall_all_pairs.tsv`
- Figure scripts: `scripts/generate_fig3_gtdb_benchmark.py`, `scripts/gtdb_census/plot_gtdb_species_ani_sv.py`, `scripts/oapgc_census/plot_oapgc_site_axes.py`, `scripts/revision/dependence_aware_uncertainty.py`

## Journal criteria

Assess the manuscript for a Nature Methods audience against originality,
scientific importance, interdisciplinary readability, technical soundness, and
methodological reproducibility. Identify only concerns supported by the packet.

## Reviewer emphasis briefs

- Reviewer 1: benchmark design, calibration provenance, leakage control,
  comparator fairness, uncertainty and dependence structure.
- Reviewer 2: census interpretation, biological claims, quality and ancestry
  confounding, cross-habitat inference, and limits of ANI-derived distance.
- Reviewer 3: clarity for nonspecialists, methods reproducibility, data and
  code availability, figure/table consistency, and claim-to-evidence mapping.

## Output contract

Each reviewer writes exactly one report to the assigned path. Reports must
follow the `nature-reviewer` skill output structure. Do not read other review
files before freezing the assigned report.

- Reviewer 1 output: `reviews/2026-10-21_R1.md`
- Reviewer 2 output: `reviews/2026-10-21_R2.md`
- Reviewer 3 output: `reviews/2026-10-21_R3.md`
