# cagPAI status summary for 528 H. pylori genomes

## Methods

- **Aligner:** minimap2 v2.31 (`-cx asm20`)
- **Markers:** 28 CDS loci from H. pylori 26695 cagPAI (HP0520–HP0547)
- **Presence rule:** coverage ≥ 0.8 AND identity ≥ 0.8
- **Classification:** complete ≥ 0.85, empty ≤ 0.15, otherwise partial

## Cohort classification counts

| status   | count | fraction |
|----------|-------|----------|
| complete |   432 | 0.818    |
| partial  |    11 | 0.021    |
| empty    |    85 | 0.161    |
| **total**|   528 | 1.000    |

## Fraction-present distribution

- mean:   0.773
- median: 0.929
- min:    0.000
- max:    1.000

## Pilot validation

| strain            | n_present | fraction_present | status   | expected        |
|-------------------|-----------|------------------|----------|-----------------|
| wt                | 28        |            1.000 | complete | complete        |
| hp26695           | 28        |            1.000 | complete | complete        |
| mut1              | 28        |            1.000 | complete | complete        |
| inv               | 28        |            1.000 | complete | complete        |
| transloc          | 28        |            1.000 | complete | complete        |
| mut1_inv          | 28        |            1.000 | complete | complete        |
| mut1_transloc     | 28        |            1.000 | complete | complete        |
| del               | 0         |            0.000 | empty    | empty           |
| mut1_del          | 0         |            0.000 | empty    | empty           |

## Threshold interpretation

Pilot controls classified correctly: 9/9. All intact-island strains (wt, hp26695, mut1, inv, transloc, mut1_inv, mut1_transloc) score fraction_present ≥ 0.857, while the two deletion mutants (del, mut1_del) score 0. The thresholds therefore cleanly separate engineered cagPAI-present from cagPAI-absent controls.
