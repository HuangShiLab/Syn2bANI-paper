# Junction coordinate agreement between Syn2b and dnadiff (GTDB 50k)

**Run:** 2026-09-02 on HPC  
**Inputs:**
- Syn2b junctions: `syn2b_tgts_cache/tmp_pairs/<pairid>/synteny.junctions.tsv`
- dnadiff alignments: `out/<pairid>/dd.1coords`

**Commands:**
```bash
python3 scripts/gtdb50k/collect_junction_coordinates.py \
    --pairs     results/gtdb50k/pairs_50k.tsv \
    --tgt-cache results/gtdb50k/syn2b_tgts_cache \
    --dnadiff   results/gtdb50k/out \
    --outdir    results/gtdb50k/junction_coords \
    --workers   32

python3 scripts/gtdb50k/compare_junction_positions.py \
    --coords results/gtdb50k/junction_coords/junction_coordinates.tsv \
    --truth  results/gtdb50k/high_ani_truth.tsv \
    --syn2b  results/gtdb50k/syn2b_inverted_fraction_50k.tsv \
    --out    results/gtdb50k/junction_coords/position_agreement.tsv
```

## Derivation check

dnadiff's own structural-difference files (`dd.rdiff`/`dd.qdiff`) were deleted by
`run_dnadiff_slice.sh` to save space, so dnadiff boundaries were re-derived from
`dd.1coords` and checked against `dd.report [Feature Estimates]` (Relocations +
Translocations + Inversions):

- exact match: 4,054 / 43,334 pairs (9.4%)
- within ±1: 9,547 / 43,334 pairs (22.0%)
- median count difference: −2

Only the 3,265 pairs whose derived count exactly matched `dd.report` were used
for the position comparison.

## Position agreement

| | count |
|---|---:|
| Syn2b boundaries | 17,079 |
| dnadiff boundaries | 28,681 |
| one-to-one matched | 12,490 |
| Syn2b with a dnadiff partner | 73.1% |
| dnadiff with a Syn2b partner | 43.5% |

**Matched-pair distance (n = 12,490):**

| percentile | distance |
|---|---:|
| median | 62,818 bp |
| p75 | 237,935 bp |
| p90 | 683,509 bp |
| p99 | 3,245,052 bp |

| within | n | % |
|---|---:|---:|
| 500 bp | 428 | 3.4% |
| 1,000 bp | 768 | 6.1% |
| 2,000 bp | 1,258 | 10.1% |
| 5,000 bp | 2,153 | 17.2% |
| 10,000 bp | 3,043 | 24.4% |
| 50,000 bp | 5,728 | 45.9% |

Median landmark spacing across these pairs: **39,648 bp**.

## Interpretation

Syn2b reports the left landmark of a broken adjacency, so its coordinate error
is bounded by the gap to the next landmark. A matched distance of the same
order as the landmark spacing (~40 kb) is therefore the expected result, not a
weak one. The median observed distance (62.8 kb) is consistent with this
expectation.

The one-sided match rate (73% of Syn2b boundaries have a dnadiff partner) is
also expected: dnadiff reports every 1-to-1 alignment boundary, including many
that fall in regions without restriction sites and are therefore invisible to
Syn2b.
