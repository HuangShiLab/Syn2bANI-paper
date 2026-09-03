# Closed-genome inversion diagnostic report

**Date:** 2026-09-04  
**Inputs:**
- `results/gtdb50k/syn2b_inverted_fraction_closed.tsv` (61,537 pairs, 701 genomes)
- `results/gtdb50k/closed_inversion_pairs.tsv` (371 seed pairs)
- `results/gtdb50k/closed_inversions_genomes.tsv` (genome metadata)

**Script:** `scripts/diagnose_closed_inversions.py`

---

## 1. The 0/371 seed-pair overlap is an identifier mismatch, not a filtering bug

The 371 seed pairs share **zero** accessions with the 701 genomes used for the
closed-genome all-vs-all run:

| Set | N accessions | GCA | GCF | Present in `closed_inversions_genomes.tsv` |
|---|---:|---:|---:|---:|
| Seed pairs (`closed_inversion_pairs.tsv`) | 662 | 162 | 500 | **0** |
| Closed-genome outputs (`syn2b_inverted_fraction_closed.tsv`) | 701 | 12 | 689 | **701** |

Even stripping version suffixes gives zero overlap. The seed pairs and the
closed-genome cohort are **different genome sets**. The original report's
"Seed pairs present in all-vs-all output: 0" therefore does not indicate a
filtering bug in the runner; it indicates that the closed-genome cohort was
downloaded/selected independently of the seed-pair list.

**Consequence:** the closed-genome all-vs-all output cannot validate the 371
seed pairs. Either the correct seed genomes must be re-downloaded and re-run,
or the validation must use the 701-genome cohort directly (e.g., by picking the
pairs with highest/lowest inverted fraction and checking them against dnadiff).

---

## 2. The 36% median inverted fraction is largely an orientation artifact

The original report quoted `syn2b_raw_inverted_fraction`, which has median
0.3595 across all 61,537 same-species pairs. The corrected column
`syn2b_inverted_fraction` (defined as `min(raw, 1 - raw)`, capped at 0.5) has
median **0.1837**:

| Metric | mean | 50% | 90% | max |
|---|---:|---:|---:|---:|
| `syn2b_raw_inverted_fraction` | 0.398 | 0.359 | 0.929 | 1.000 |
| `syn2b_inverted_fraction` (corrected) | 0.187 | 0.184 | 0.349 | 0.500 |

23,908 / 61,537 pairs (38.9%) have `raw > 0.5`; for these the corrected value is
`1 - raw`. The high raw median is therefore driven by arbitrary global
orientation: when two assemblies are randomly oriented with respect to each
other, ~50% of the genome is classified as inverted. The corrected metric
removes this confounder.

The top 20 "most inverted" pairs by raw value are all *Pseudomonas aeruginosa*
and have `raw = 1.0000`. Under the corrected metric these pairs have
`inverted_fraction = 0.0000`, i.e. they are perfectly collinear but opposite
orientation. This strongly suggests that the high raw values record coordinate
conventions, not biological inversions.

---

## 3. Species-level patterns

Median `syn2b_raw_inverted_fraction` by species (top 4):

| Species | n pairs | median raw | median corrected |
|---|---:|---:|---:|
| *Salmonella enterica* | 19,900 | 0.500 | 0.168 |
| *Bordetella pertussis* | 19,900 | 0.427 | 0.143 |
| *Pseudomonas aeruginosa* | 19,900 | 0.110 | 0.055 |
| *Streptococcus pneumoniae* | 5,050 | 0.082 | 0.041 |

The species ranking reverses partially after correction. *S. enterica* goes from
median 0.50 (raw) to 0.17 (corrected), confirming that the raw median is
inflated by orientation convention. The residual corrected values still need
biological interpretation; they may reflect real rearrangement rates or remaining
coordinate artifacts (e.g., circular start positions not normalized).

---

## 4. Correlation with assembly fragmentation is negligible

Because the cohort was selected as near-complete genomes, contig count is not
the driver:

| | raw inverted fraction | q_contigs | r_contigs |
|---|---:|---:|---:|
| raw inverted fraction | 1.000 | 0.002 | -0.020 |
| q_contigs | 0.002 | 1.000 | 0.134 |
| r_contigs | -0.020 | 0.134 | 1.000 |

The main confounders are therefore **global orientation** and possibly
**circular start-coordinate convention**, not fragmentation.

---

## 5. Recommendations

1. **Use `syn2b_inverted_fraction` (corrected) for undirected all-vs-all
   comparisons.** The raw column is reference-oriented and only meaningful when
   one genome is a fixed biological reference.
2. **Resolve the seed-pair mismatch.** Either re-download the 662 seed accessions
   and re-run Syn2b, or abandon the seed-pair validation and instead validate the
   701-genome cohort against dnadiff on the top/bottom inverted-fraction pairs.
3. **Add circular-origin normalization to the closed-genome pipeline.** Even the
   corrected metric may be contaminated by arbitrary start coordinates on circular
   chromosomes. Rotating each circular genome to a common start (e.g., dnaA or
   smallest lexicographic k-mer) before comparison should reduce the residual
   median further.
4. **Re-generate `CLOSED_GENOME_INVERSION_REPORT.md`** using the corrected
   metric after the seed-pair and circular-origin issues are fixed.
