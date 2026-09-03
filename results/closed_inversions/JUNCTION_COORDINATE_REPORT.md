# Closed-genome inversion cohort: junction coordinate agreement

## Cohort

- 701 near-closed genomes from GTDB r207 selected by `select_closed_genomes.py`:
  - *Streptococcus pneumoniae* 101
  - *Salmonella enterica* 200
  - *Bordetella pertussis* 200
  - *Pseudomonas aeruginosa* 200
- 64,750 all-vs-all pairs computed with Syn2b (four-enzyme panel).
- Top 100 pairs with `raw_inverted_fraction > 0.3` and `breakpoints > 0` selected for
  dnadiff + coordinate comparison (HPC job 3982403, completed 2026-09-03).

## Methods

- dnadiff run on each pair, keeping `dd.report` and `dd.1coords`.
- `collect_junction_coordinates.py` derives rearrangement boundaries from `dd.1coords`
  (INV + JMP) and compares the derived count to `dd.report`.
- `compare_junction_positions.py` performs greedy one-to-one matching between Syn2b
  and dnadiff reference-frame junction coordinates.

## Results (all 100 pairs, `--all-pairs`)

| metric | value |
|--------|-------|
| pairs | 100 |
| Syn2b boundaries | 483 |
| dnadiff boundaries | 3,548 |
| matched one-to-one | 483 |
| Syn2b boundaries with dnadiff partner | 100.0% |
| dnadiff boundaries with Syn2b partner | 13.6% |
| median matched distance | 16,396 bp |
| p75 matched distance | 105,003 bp |
| p90 matched distance | 314,810 bp |
| within 5 kb | 40.4% |
| within 50 kb | 62.1% |
| median landmark spacing (4-enzyme) | ~440 bp |

Only 3/100 pairs passed the strict derivation check (`count_diff == 0`);
`dd.report` counts and the `dd.1coords`-derived counts rarely agree exactly.

## Interpretation

1. **Sensitivity gap**: dnadiff reports ~7× more boundaries than Syn2b. Many of these
   are small collinearity breaks (e.g., low-identity blocks, prophage regions) that
   fall below the landmark-spacing detection floor of Syn2b.
2. **Coordinate agreement is worse than landmark spacing**: the median matched
   distance (16 kb) is ~37× the median landmark spacing (~440 bp). This is not the
   ~spacing-level agreement seen in controlled E. coli K-12 simulations.
3. **Every Syn2b boundary has a dnadiff partner**: in the greedy matching, all 483
   Syn2b boundaries find a nearest dnadiff boundary. This is expected when dnadiff
   is much denser; it does not mean the two methods call the same event.

## Why the simulations and real divergent pairs differ

In the E. coli K-12 shattering controls cited in `compare_junction_positions.py`,
the *same* large inversion is present in both genomes and there are no other
rearrangements. In real divergent closed genomes, multiple overlapping events
(fragmentation of homology by indels, prophage insertions, repeated inversions)
create many dnadiff boundaries, and greedy one-to-one matching pairs Syn2b
landmarks with whichever dnadiff boundary happens to be closest, often a different
biological event.

## Recommendation

Coordinate-level validation on highly rearranged real pairs is not the strongest
way to demonstrate Syn2b's accuracy. Better strategies:

1. **Simulated controls**: known inversions/relocations of defined size in an
   otherwise identical background (already done; spacing-level agreement holds).
2. **Count-level correlation**: regress Syn2b `breakpoint_count` on dnadiff
   `Breakpoints` / `Relocations + Inversions` across the full 64k closed pairs.
3. **Phenotype cohorts**: link high-synteny-discordant pairs to known phenotypes
   (antibiotic resistance, virulence, host range) where the biological signal is
   at the strain-cluster level, not the single-junction level.

## Files

- `results/closed_inversions/junction_coordinates.tsv` — per-pair boundary counts and
  raw coordinates from HPC.
- `results/closed_inversions/position_agreement_allpairs.tsv` — one-to-one matching
  summary for all 100 pairs.
- HPC dnadiff outputs: `/lustre1/g/aos_shihuang/data/closed_inversions/top100_dnadiff/`
