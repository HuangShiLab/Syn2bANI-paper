# HROM within-species structural census

**Completed:** 2026-09-26.
**Tool:** Syn2b structural channel, enzyme panel `BcgI,AlfI,AloI,FalI`; ANI from skani 0.3.2.
**Scope:** all HROM conspecific genomes in multi-genome clusters.

## Census size

| Item | Value |
|---|---:|
| HROM species clusters | 5,113 |
| Total HROM genomes | 145,149 |
| Multi-genome clusters | 2,241 |
| Genomes included | 142,277 |
| Unique unordered pairs | 62,872,901 |
| Structural tasks | 3,525 |
| Structural task failures | 0 |
| Missing clusters/pairs | 0 |
| Pair rows with matched skani ANI | 62,872,901 / 62,872,901 |

## QC

- All 3,525 block-pair tasks completed and produced checkpoints.
- The merged unique-pair count is exactly `sum(C(n,2))` over included clusters:
  **62,872,901**.
- No cluster has fewer pairs than its genome-level expectation.
- Every structural pair received an ANI value.
- The first block-pair implementation retained within-block rows in cross-block
  tasks; the merge step suppressed **10,644,750** duplicate rows. The final
  pair table is therefore complete and non-redundant. The worker filter was
  corrected for future runs (`bi != bj` now requires the second genome to lie
  in the partner block).
- HROM FASTA headers are contig IDs, so digestion used
  `--ensure-filename-genome-id` to key all contigs by their HROM genome
  accession.

## Overall structural prevalence

Using the exact species-level counts:

| Metric | Value |
|---|---:|
| Pairs with `breakpoints >= 2` | 17,535,013 |
| Fraction with `breakpoints >= 2` | 27.8896% |

Thus, over half of HROM within-species pairs have ANI >= 97%, and about 28% of
all within-species pairs carry at least two alignment-visible block-level
junctions. This is a structural layer that an ANI-only comparison does not
provide.

## Results by skani ANI bin

| ANI bin | pairs | % with breakpoints >= 2 | median breakpoints | mean breakpoints |
|---|---:|---:|---:|---:|
| <95% | 2,376,169 | 45.128 | 1 | 8.016 |
| 95–97% | 28,516,234 | 31.085 | 0 | 1.698 |
| 97–98% | 22,576,854 | 20.837 | 0 | 1.080 |
| 98–99% | 8,946,967 | 31.610 | 0 | 1.719 |
| 99–99.5% | 363,255 | 14.104 | 0 | 0.581 |
| >=99.5% | 93,418 | 15.751 | 0 | 0.628 |

The increase at 98–99% relative to 97–98% shows that ANI rank alone is an
imperfect proxy for structural divergence. The `raw_inverted_fraction` column
should not be used as a stand-alone biological signal for these mostly draft
assemblies: its median is near 0.5 in every bin, as expected when arbitrary
contig orientation dominates that ratio. The primary interpretable count here
is block-level `breakpoints`.

## Top clusters by number of pairs with >=2 breakpoints

The complete per-cluster table is in `species_summary.tsv`. The largest absolute
contributions include:

| Cluster | genomes | pairs | pairs with >=2 breakpoints | % |
|---|---:|---:|---:|---:|
| HROM_Genome_1099 | 3,123 | 4,875,003 | 1,839,165 | 37.726 |
| HROM_Genome_0976 | 2,627 | 3,449,251 | 1,693,482 | 49.097 |
| HROM_Genome_0999 | 2,151 | 2,312,325 | 1,297,490 | 56.112 |
| HROM_Genome_4823 | 2,316 | 2,680,770 | 1,091,165 | 40.703 |
| HROM_Genome_4892 | 2,196 | 2,410,110 | 667,294 | 27.687 |
| HROM_Genome_0302 | 2,472 | 3,054,156 | 584,115 | 19.125 |
| HROM_Genome_4825 | 1,186 | 702,705 | 556,331 | 79.170 |

## Files

- `hrom_within_species_sv.tsv.gz`: structural census table (1.9 GB compressed).
- `hrom_within_species_sv_ani.tsv.gz`: structural table plus `skani_ani`
  (2.0 GB compressed).
- `species_summary.tsv`: per-cluster structural summary.
- `ani_bin_summary.tsv`: ANI-bin summary.
- `census_qc.md`, `merge_ani_qc.md`: QC reports.

The two large compressed tables are intended for Zenodo/figshare rather than
GitHub; this directory currently carries summaries and QC in Git.
