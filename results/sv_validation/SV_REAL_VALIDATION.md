# SV validation on real genomes: syn2bani `struct` vs dnadiff

Date: 2026-08-14. Code: Syn2bANI main @ 98177dc + two working-tree fixes made
during this validation (unchained-anchor breakpoint fix in `chain_ani.rs`,
relocation guard in `sv.rs`; both uncommitted, `cargo test --lib` = 73 green).
Binary: `target/release/syn2bani` (Mac Studio). Truth: dnadiff (MUMMER/3.23)
on hpc2021.hku.hk, `dnadiff -p out <ref> <query>` per pair; inputs uploaded
and outputs mirrored under `results/sv_validation/<pair>/`
(`out.1coords/qdiff/rdiff/report`, `time.log`, `struct_sv.tsv`, `chains.paf`).
Analysis: `analysis/sv_real_validation.py` (run to regenerate
`sv_summary.tsv`, `sv_indel_compare.tsv`, `sv_inversion_compare.tsv`,
`compare_output.txt`). Figure: `analysis/plot_report_figures.py::fig8_sv_detection`.

## Pairs

| pair | ANI (syn2bani het) | why |
|---|---|---|
| E. coli K-12 MG1655 vs W3110 | 99.99 | near-clonal control; carries the documented ~780 kb rrnD–rrnE inversion |
| K-12 MG1655 vs O157:H7 Sakai | 98.13 | collinear enteric backbones, many prophage indels |
| S. Typhi CT18 vs Typhimurium LT2 | 98.48 | heavily rearranged (multiple nested inversions + relocations) |

All `struct` runs used the default 4-enzyme panel (BcgI,AlfI,AloI,FalI),
`indel_min = 1000`. Truth model: inversions = clustered reverse-orientation
blocks in `out.1coords` (merge gap 50 kb); indels ≥1 kb = `GAP` (|diff|),
`DUP` (copy number), and `BRK` (unaligned stretch — how dnadiff surfaces
prophage/mobile regions) entries in qdiff/rdiff, GAPs deduplicated across the
two files. Matching is span-based and one-to-many, because dnadiff fragments
large accessory regions into many small events while syn2bani calls the whole
junction once.

## Headline numbers (`sv_summary.tsv`)

| pair | inversions (syn / truth clusters) | translocations | indel calls | TP | FP | truth ≥1 kb | covered | missed (≥5 kb) |
|---|---|---|---|---|---|---|---|---|
| MG1655 vs W3110 | 1 / 1 | 0 | 10 | 10 | 0 | 13 | 10 | 3 (0) |
| MG1655 vs Sakai | 0 / 0 large (7 tiny repeat flips ≤2 kb) | 1 | 121 | 112 | 9 | 197 | 174 | 23 (9) |
| CT18 vs LT2 | 8 / all inside 2 rearranged regions | 4 | 70 | 62 | 8 | 104 | 70 | 34 (9) |
| **total** | | | **201** | **184 (91.5%)** | **17** | **314** | **254 (80.9%)** | **60 (18)** |

Size accuracy on the 128 one-to-one matched indels: median size ratio
**1.000**, 121/128 (94.5%) within ±25%, 10–90th percentile 0.997–1.086.

## Per-pair findings

### MG1655 vs W3110 (control)

- **Inversion: exact.** syn2bani q 3,429,150–4,207,770 (778,620 bp); dnadiff
  reverse cluster 3,423,662–4,213,165 (789,503 bp). Endpoints sit 5.4–5.5 kb
  inside the truth span; the truth junctions lie in the rrnD/rrnE repeat arms
  (~3.2 kb near-identical copies) where no tag anchor can phase uniquely, so
  this offset is at the resolution limit of the method.
- **Indels: 10/10 correct, 0 false.** Includes the 6,790 bp MG1655-specific
  island (size match exact) and nine ~1.2–1.3 kb IS-element copy-number
  differences (sizes match dnadiff DUP lengths within 0.5%).
- **3 missed**, all IS-copy DUPs (1,203–1,339 bp) in tandem-repeat contexts
  where the extra copy's anchors multi-map and cannot be chained.

### MG1655 vs Sakai

- **Inversions: none called — correct.** dnadiff's 14 "Inversions" in
  `out.report` are all reverse-alignment fragments ≤2 kb at 84–97% identity
  (inverted-repeat matches), below the 4-anchor chain minimum. syn2bani
  reports none of them.
- **1 translocation call**, verified against `out.1coords`: a 12.5 kb query
  block at q 0.57 Mb whose homolog sits at r 1.62–1.63 Mb — a genuine small
  relocation.
- **Indels: 112 TP / 9 FP / 23 missed.** All 9 FPs overlap dnadiff `JMP`
  (relocation-jump) features: the reference segment exists in the query but
  *elsewhere*, so dnadiff says relocation while syn2bani reports the local
  junction as an indel. None is a fabricated event. The 23 missed are small
  (16 < 5 kb; largest 21 kb `BRK`), inside anchor-poor accessory/prophage
  regions where flanking chains drop below the 4-anchor minimum.

### Typhi CT18 vs Typhimurium LT2

- **8 inversions + 4 translocations, all true-region.** nucmer finds two
  heavily rearranged regions (reverse/relocation blocks covering ref
  ~1.3–1.9 Mb and ~3.6–4.4 Mb, 95 reverse blocks in 10 clusters). Every
  syn2bani inversion/translocation call falls inside these clusters; no
  rearrangement is called anywhere else on the genome. syn2bani subdivides
  the regions into nested events (e.g. q 3.42–4.26 Mb split into 5 inversions
  + junctions), consistent with the documented nested-inversion history of
  this pair; per-event endpoint comparison against dnadiff is not meaningful
  there because dnadiff itself fragments the regions into ~95 blocks.
- **Indels: 62 TP / 8 FP / 34 missed.** As with Sakai, all 8 FPs overlap
  dnadiff `JMP` relocations; missed events are small (median ~3 kb; largest a
  64 kb `BRK` insertion at q 4.44 Mb and a 28 kb `BRK` deletion), again in
  anchor-poor accessory sequence.

## Runtime (`sv_summary.tsv`)

| pair | syn2bani struct | dnadiff | speedup |
|---|---|---|---|
| MG1655 vs W3110 | 48 ms | 8.24 s | 172× |
| MG1655 vs Sakai | 54 ms | 9.59 s | 178× |
| CT18 vs LT2 | 50 ms | 8.52 s | 170× |

syn2bani: Mac Studio, min of 5 runs, end-to-end from FASTA (digest + chain +
call). dnadiff: HPC login node, `/usr/bin/time -v` elapsed, full pipeline
(nucmer, delta-filter, show-diff, show-snps). Machines differ; the ~170× gap
is far larger than any hardware effect.

## Bugs found and fixed during this validation

1. **`breakpoint_count` counted chaining-rejected anchors** (Part 1;
   `results/sv_validation/breakpoint_formula.md`).
2. **Relocation junctions were misreported as giant indels.** In `sv.rs`
   branch (b), a small relocated block sandwiched between collinear chains
   produced a "deletion" sized by the reference jump distance (observed: a
   12.5 kb relocated block on MG1655 vs Sakai → spurious 961 kb deletion;
   five such calls across the three pairs, sizes 42 kb–961 kb). Fix: the
   between-chain indel branch now requires that no third chain anchors inside
   the reference interval the junction spans (`interval_occupied` guard);
   the order-contradiction junction is still reported as a translocation.
   Unit test `relocation_junction_is_not_a_giant_indel`; `cargo test --lib`
   73 passed. After the fix the largest indel call on Sakai is 116 kb and
   corresponds to a real Sakai-specific prophage region.

## Caveats (for reviewers)

- Equal-length homologous (HGT) replacements keep tags phased and do not
  break chains: they are not reported. Correct by design, but such events are
  invisible to `struct`.
- SVs whose flanking chains have < 4 anchors are not called. This is the
  dominant source of missed small events inside accessory regions
  (all 60 missed events; 42 of them < 5 kb).
- dnadiff and syn2bani disagree on classification at relocation junctions:
  all 17 calls without a dnadiff indel counterpart overlap dnadiff `JMP`
  features. Counting them as errors is a worst-case reading; none is a
  fabricated event.
- Junction-level granularity differs: dnadiff fragments large accessory
  regions into many GAP/BRK/DUP pieces; syn2bani reports one call per
  junction. All matching above is one-to-many to account for this.
- Inversion endpoints are limited by tag spacing and unphaseable repeat
  copies at the junctions (observed ≤5.5 kb on the 779 kb W3110 inversion).
