# B. longum abfA case study — analysis report

Date: 2026-08-26
Data: 185 JNU B. longum genomes (HPC `Strain2b/data/JNU_genomes/genome2023/fna`),
all-vs-all Syn2bANI (17,020 pairs, `b_longum_s2b_matrix.tsv`),
metadata `metadata.tsv` (abfA = cluster2, hypba = cluster1; phenotype unknown
except 5 representative strains not yet provided).

## Headline numbers

- All 17,020 pairs have ANI ≥ 98.67% (16,313 pairs ≥ 99%): ANI cannot
  distinguish any of these strains.
- af_query spans 0.598–0.991 across pairs, but whole-genome AF does **not**
  separate abfA-complete from abfA-deleted strains (median 0.861 vs 0.867
  against the abfA+ reference FSHHK16M1; MWU p = 0.66): a ~20 kb island is
  swamped by accessory-genome and assembly variance at whole-genome scale.

## Locus-targeted analysis (vs abfA+ reference FSHHK16M1, 20 contigs)

abfA locus mapped by blastn of the curated cluster genes
(`gene_list_abfa.txt`, 1,601 homologs) to FSHHK16M1: contig10,
8,546–37,075 (28.5 kb). The six canonical AHWH21M4 cluster genes map to
the core 8,840–17,163.

Chain-coverage profiles (`struct --paf`, 1 kb bins, `abfA_locus_figure.pdf`):

| Region | complete median | deleted median | MWU p |
|---|---|---|---|
| core (8.8–17.2 kb) | 1.000 | 0.913 | 6.6e-04 |
| periphery (17.2–37.1 kb) | 1.000 | 0.160 | **1.4e-03** |
| downstream (37.1–113 kb) | 0.784 | 0.787 | 0.40 (n.s.) |

- The coverage difference is **specific to the cluster periphery**; the
  rest of the contig is indistinguishable — a genuine locus-scale signal.
- As a classifier (periphery coverage < 0.5): sensitivity 52%,
  specificity 75%. Not a clean discriminator.

## Why the discrimination is imperfect (verified facts)

1. The six canonical cluster genes are BLAST-present full-length (≥96%
   identity) in curated "deleted" strains as well — the curated binary
   call reflects peripheral/partial gene content, not physical absence of
   the core cluster. The structural reality is a continuum of partial
   haplotypes, so a binary metadata label cannot be recovered cleanly.
2. Chain spans bridge deleted intervals by design (the gap-penalised DP is
   what makes ANI robust to indels), so chain-span coverage overestimates
   base-level presence; struct Deletion events at the locus fire in only
   14/56 curated-deleted genomes on these draft assemblies.

## Conclusions for the manuscript

- ANI-blindness confirmed on a real, functionally annotated strain
  collection (phenotype-linked locus).
- Syn2bANI chain coverage localises the functional difference to the abfA
  cluster periphery with locus specificity (p = 1.4e-3), something no
  whole-genome ANI/AF metric shows (p = 0.66).
- But the effect is a moderate association, not a clean classifier —
  present it as a demonstration of locus-level structural signal, not as
  "Syn2bANI predicts abfA status". The clean mechanistic demonstrations
  remain the cagPAI pilot (exact 36,154 bp calls) and the GTDB discordant
  pairs (Mb-scale rearrangements at ANI ≥ 95%).

## Files

- `abfA_region_coverage.tsv` — per-strain core/periphery/downstream coverage
- `abfA_locus_coverage.tsv`, `abfA_locus_coverage_multiref.tsv` — wide-locus
  and 3-reference variants (superseded by the region-resolved analysis)
- `abfA_locus_figure.pdf/.png` — coverage profile + periphery boxplot
- `abfA_pair_metrics.tsv`, `ABfA_ANALYSIS_REPORT.md` — pair-level base stats
- HPC: `struct_vs_ref/`, `paf_vs_ref/`, `paf_vs_FGDLZ73M1_ctg/`,
  `paf_vs_FHNBA1M1_ctg/`, `blast/` under `results/b_longum_abfA/`
