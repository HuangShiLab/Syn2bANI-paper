# H. pylori cagPAI cohort analysis — Syn2bANI case study

Goal: show that ANI-based genome search misses clinically relevant structural
differences, using cagPAI (36 kb pathogenicity island, complete/partial/empty)
in the Song et al. 2026 PNAS cohort (528 genomes, Correa's cascade stages)
as the functional phenotype.

Data sources:
- New genomes: NCBI BioProject PRJNA1419121
- Public genomes: ~65 BioProjects listed in the paper's Data Availability
- Per-strain metadata (accession × disease stage NAG/AG/IM/GC): SI Table S1
  (request from Prof. Liang Wang; place at `metadata.tsv` next to this README)

Requires syn2bani >= 4d679cf (struct shadow-chain fix).

## Workflow

1. `s1_prepare_genomes.sh accessions.txt outdir/` — download assemblies
   (one accession per line, GCA_/GCF_ accessions from SI Table S1).
2. `s2_struct_vs_26695.slurm` — `syn2bani struct` of every genome vs
   H. pylori 26695 (cagPAI+) and vs the in-silico empty-site construct;
   classifies each strain's cagPAI state (complete / partial / empty) from
   event sizes at the cagPAI locus (chr 547,328–583,481 in 26695).
3. `s3_triangle_dist.slurm` — all-vs-all `syn2bani triangle --edge-list`
   for ANI + af asymmetry screen.
4. `analyze_cagpai_cohort.py` — merge cagPAI states + pairwise ANI/af +
   disease-stage metadata; stratified statistics and figures:
   - rate of cagPAI-discordant pairs within ANI ≥ 95 / 97 / 99 % bins
     (the pairs an ANI search would return as "same strain");
   - ANI vs |Δaf| scatter, discordant pairs highlighted;
   - disease-stage crossing rate of discordant pairs
     (NAG/AG vs IM/GC: functional consequence).
