# ANI–synteny discordance analysis (first pass)

## Data sources

- `results/gtdb50k/s2b_50k.tsv` + `truth_50k.tsv`: 43,334 held-out pairs
- `results/gtdb50k/high_ani_results.tsv`: 727 high-ANI test pairs (95–100%)

## Discordance definition

Pairs with ANIm truth ≥ threshold but Syn2bANI synteny_score < 0.98 (i.e.
near-clonal ANI but detectable structural divergence).

## Counts

| Set | ANI ≥ 95% | ANI ≥ 97% | ANI ≥ 98% | ANI ≥ 99% |
|---|---|---|---|---|
| 43,334 held-out | 211 | — | — | 1 |
| 727 high-ANI test | 176 | 147 | 122 | 73 |
| **Combined** | **387** | — | — | **74** |

The 43,334 held-out set is sparse above 95% ANI because it was stratified
toward divergent bands; the high-ANI test set provides most of the
near-clonal discordant cases.

## Top discordant pairs (ANI ≥ 99%)

| pairid | ANIm | synteny | breakpoints |
|---|---|---|---|
| GCA_002725245.1__GCA_013330135.1 | 99.72 | 0.917 | 67 |
| GCA_002431505.1__GCA_003488365.1 | 99.89 | 0.921 | 393 |
| GCA_002344645.1__GCA_003501905.1 | 99.95 | 0.931 | 327 |
| GCA_002295625.1__GCA_003513645.1 | 99.82 | 0.935 | 166 |
| GCA_903824795.1__GCA_903958965.1 | 99.95 | 0.938 | 170 |
| GCA_903937715.1__GCA_903945645.1 | 99.58 | 0.939 | 190 |
| GCA_900761705.1__GCA_905197205.1 | 98.88 | 0.936 | 155 |
| GCA_008363445.1__GCA_013360885.1 | 99.30 | 0.942 | 148 |
| GCA_903820925.1__GCA_903860705.1 | 99.99 | 0.950 | 111 |
| GCA_903863355.1__GCA_903872275.1 | 99.72 | 0.951 | 410 |

## Taxonomic distribution

The most discordant high-ANI pairs are enriched in:
- *Proteobacteria* (especially *Escherichia*, *Salmonella*, *Dickeya*,
  *Serratia*, *Burkholderia*)
- *Actinobacteriota* (*Streptomyces*)
- *Firmicutes_A*
- *Cyanobacteria*

This matches known high-rearrangement-rate taxa.

## Output files

- `results/gtdb50k/discordant_ani95_syn98.tsv` — 211 held-out pairs
- `results/gtdb50k/discordant_high_ani_test_syn98.tsv` — 176 high-ANI test pairs

## Next steps

1. If 387 preliminary pairs are sufficient for the manuscript figure, proceed
   to run `syn2bani struct` on the top 50–100 cases and compare against
   dnadiff/minimap2 SV truth.
2. If more 95–100% pairs are needed, sample additional same-species or
   same-genus pairs from GTDB-R207 non-representative genomes and run the
   full pipeline (Syn2bANI + ANIm truth + dnadiff SV truth).
3. Integrate the B. longum abfA case as a functionally validated example.
