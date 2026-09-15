# ecoli_o157_fitzgerald_2021: re-run with syn2bani 0.1.0 (639b31d)

- genomes: 74; non-self pairs: 2,701
- ANI range: 99.8858 – 99.9999
- breakpoint_count: median 7, IQR 4–9, max 19, pairs at 0: 7.1%
- synteny_blocks: median 59
- Spearman(ANI, breakpoint_count) = -0.449 (p = 1.9e-134)
- pairs with ANI ≥ 99.9% and breakpoint_count ≥ 10: 652

## breakpoint_count by assigned_lineage

- same-assigned_lineage pairs: n = 1,376, median 6
- different-assigned_lineage pairs: n = 1,325, median 9

| assigned_lineage (both members) | pairs | median | IQR | max |
|---|---:|---:|---:|---:|
| I/II | 55 | 0 | 0–6 | 13 |
| II | 1275 | 6 | 3–8 | 18 |
| Ia | 45 | 6 | 3–9 | 13 |
| Ic | 1 | 0 | 0–0 | 0 |

## breakpoint_count by host_category

- same-host_category pairs: n = 1,426, median 6
- different-host_category pairs: n = 1,275, median 8

| host_category (both members) | pairs | median | IQR | max |
|---|---:|---:|---:|---:|
| bovine | 1275 | 6 | 3–9 | 19 |
| human | 15 | 9 | 7–12 | 17 |
| other/unknown | 136 | 8 | 4–9 | 17 |

