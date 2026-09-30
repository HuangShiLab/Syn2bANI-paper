# OAPGC local within-species census feasibility

## Data selection

Zenodo metadata used:

- record 19212794, genome sets 1–2 plus `cluster.annotation.xlsx`;
- record 19220405, genome sets 3–5 plus annotations;
- CC-BY-4.0.

The annotation contains 149,921 HQMAG rows. A non-null
`Represent(STRAINs)` field identifies **99,215 strain-level non-redundant
representatives**. Grouping these by SGB gives **1,488 multi-genome clusters**;
98,229 representatives enter pairwise analysis and yield **30,560,247 unique
unordered pairs**.

The raw genome tarballs are 67.3 GB compressed. Only representative FASTAs will
be extracted; raw files remain outside Git.

## Local resource estimate

| Metric | Value |
|---|---:|
| Apple Silicon local cores used | 12–16 |
| structural tasks (batch 250) | 2,120 |
| planning Syn2b cost | 405.6 core-h |
| digestion cost | <0.01 core-h |
| elapsed estimate at 16 cores | 25.4 h |
| realistic wall time including I/O | 1–2 d |
| gzipped structural output | 0.92 GB |
| likely extracted FASTA footprint | 60–100 GB |
| local free space at estimate | 1.2 TB |

Conclusion: **yes, a local within-species structural census is feasible**, but
should proceed as a largest-20-cluster smoke test before the full run. If the
smoke test shows materially >13 core-h/million pairs because OAPGC MAGs are more
fragmented than HROM, switch the full run to HPC and retain only local analysis.

## Current status

- Zenodo records inspected and metadata table downloaded.
- Five representative-only genome tarball downloads started.
- Census planning code added under `scripts/oapgc_census/`.
- Full planning files generated outside Git in `/Users/macstudio/Downloads/OAPGC/census_plan`.
