# Mock peer review — Syn2b-ANI: strain-level ANI and structural comparison from fixed restriction-site anchors

*Reviewer perspective: Nature Methods, computational biology / microbiome methods.*

---

## Reviewer summary

Syn2bANI proposes a novel restriction-site-anchor framework that simultaneously estimates ANI and reports synteny/structural metrics from the same chaining pass. The idea is elegant, the exact-truth simulation suite is thorough, and the GTDB-R207 held-out benchmark is large. However, the manuscript currently mixes a strong ANI-estimation story with a less mature structural-comparison story. A post-review analysis has additionally revealed that the flagship *H. pylori* cagPAI case study is confounded by circular-chromosome coordinate-origin artifacts, which changes the revision priorities. Before publication in *Nature Methods*, I would ask the authors to (i) sharpen the claim of what is fundamentally new compared with sparse k-mer chaining, (ii) restructure the structural section around the metrics that are empirically robust (`raw_inverted_fraction`), (iii) fix circular-origin handling in the `struct` pipeline, and (iv) make the output-selection logic (raw vs calibrated vs hybrid) transparent. I assess the work as **major revision**.

---

## Major concerns

### 1. Innovation: what is the *fundamental* advance over sparse k-mer chaining?

The manuscript argues that Type IIB tags are "long enough that a match is almost certainly homologous" and that the *absence* of a tag is informative. These properties are shared by sufficiently long k-mers, and skani already uses sparse k-mer chaining. The real differentiator appears to be that anchors carry exact coordinates, enabling structural descriptors in the same pass. The authors should state this more crisply and directly compare coordinate-carrying anchors with skani's chain statistics: can skani's chains (or post-processing of skani output) report analogous breakpoint/inversion metrics? If not, why not? If yes, the advantage must be framed in accuracy or computational cost.

### 2. The structural story is using the wrong headline metric

The manuscript currently leads with `breakpoint_count`, `anchor_adjacency`, and `synteny_blocks`. A recent reanalysis shows that:

- `raw_inverted_fraction` correlates strongly with dnadiff (Pearson r ≈ 0.936 across GTDB-R207; r ≈ 0.996 in the strain range), is fragmentation-invariant by construction, and is the most reliable structural output.
- `breakpoint_count` is useful only after controlling for assembly fragmentation and should be validated against `dnadiff relocations + inversions`, not `dnadiff breakpoints` (the latter is itself contaminated by contig-count artifacts).
- `synteny_blocks` is dominated by assembly fragmentation and should be reported as an assembly-structure diagnostic, not a rearrangement metric.
- `af_query` is a coverage/aligned-fraction metric and should be framed as such.
- `anchor_adjacency` measures local anchor-order conservation and is not a proxy for base-pair synteny.

The manuscript should be rewritten so that `raw_inverted_fraction` is the headline structural output, `breakpoint_count` is a secondary count-based descriptor, and the other metrics are clearly mapped to their interpretive niches. The mathematical distinction (count-based vs ratio-based metrics under fragmentation) should be summarized in the Supplementary Information.

### 3. Circular-chromosome coordinate-origin artifacts invalidate the cagPAI rearrangement analysis

In the 528-genome *H. pylori* cohort, 95.2% of genomes classified as `complete_rearranged` carry a structural call spanning >50% of the chromosome, most commonly a translocation `1459–1666206` that covers 99.8% of the *H. pylori* 26695 genome. This is not a biological translocation; it is the consequence of comparing a circular chromosome to a reference with a different arbitrary start coordinate. Because the cagPAI window (547,327–583,481) is always overlapped by a genome-spanning call, the `complete_rearranged` state is largely an artifact.

This has two consequences for the manuscript:

- The disease-stage and FastBAPS associations driven by `complete_rearranged` must be recomputed after filtering out genome-spanning calls and/or normalizing circular start coordinates. Preliminary data suggest the presence/absence axis (`empty` / `partial` / `complete`) is unaffected and remains biologically interesting.
- The Syn2bANI `struct` pipeline needs a fix: per-contig circularity handling, circular-origin normalization before comparison, and a filter that flags any call spanning more than a fixed fraction of the genome (e.g., 50%) as a coordinate-system artifact rather than a rearrangement.

Until this is fixed, the cagPAI case study cannot stand as evidence for SV-based phenotype association.

### 4. The hybrid estimator threshold and the calibration black-box

The 98% hybrid threshold is justified by a threshold sweep, but the manuscript reports only MAE values at 97/98/99%. The sweep itself is not shown. Because the hybrid rule is the proposed default for mixed input, readers need to see the sweep and the robustness to training-set composition. Relatedly, the ridge calibration is trained on 2,520 complete GTDB-R207 pairs. While the authors are admirably transparent that it fails on MAGs and mid-ANI pairs, the calibration is still a post-hoc ML layer whose features are not independently interpretable. A reviewer will ask whether the calibration is simply memorizing GTDB-specific rate-heterogeneity patterns and whether it would generalize to future GTDB releases or other genome collections.

### 5. Structural outputs are under-validated against dedicated SV tools

Even after the above corrections, the structural validation relies almost entirely on dnadiff-derived breakpoints and coverage scores. dnadiff is a reasonable reference for ANI and aligned fraction, but it is not a gold standard for SV recall/precision. The engineered *H. pylori* cagPAI panel is excellent, but it tests only one deletion/insertion size. I would expect at least one comparison with a dedicated pangenome/structural tool such as SyRI, minigraph-Cactus, or progressiveMauve on the Enterobacteriaceae completes or a subset of Syntracker discordant pairs, reporting per-event precision/recall for inversions, translocations, and indels >1 kb. Without this, claims such as "inversion/translocation calls" are not fully supported.

### 6. Closed-genome all-vs-all inversion distribution needs explanation

A closed-genome cohort of 701 near-complete GTDB-R207 genomes (64,750 pairs) shows a median `raw_inverted_fraction` of 0.36, with the 90th percentile at 0.93 and the 95th at 1.0. A 36% median inversion rate among conspecific pairs is biologically implausible and likely reflects the same coordinate-origin arbitrariness seen in *H. pylori*, compounded by global strand/orientation conventions. Additionally, 371 expected seed pairs appear 0 times in 61,537 outputs, suggesting an identifier or filtering bug. These two issues must be diagnosed before `raw_inverted_fraction` can be presented as a trustworthy headline metric.

---

## Minor concerns

1. **Terminology.** The manuscript uses "ANI" interchangeably for raw gated, calibrated, and hybrid estimates. Consider reserving "ANI" for the hybrid/default output and labeling raw/calibrated columns explicitly.

2. **Figure clarity.** Figure 1 should visually contrast Syn2bANI with skani/FastANI: k-mer methods also chain seeds, so the difference is not obvious from the current schematic. Figure 6's annotation of top-discordant cases is hard to read at manuscript size.

3. **Supplementary completeness.** Some supplementary figures (e.g., S10, S11, S13) are referenced in the numbering but absent from the Supplementary document. The math review and detection-power documents are cited but hosted in a different repository; they should be included as Supplementary Notes or stable links.

4. **MAG recommendation.** The recommendation to use raw gated output for MAGs is buried in Results. It should be in the abstract/Discussion and enforced in the CLI help/default behavior.

5. **Code availability.** The repository is available, but a Dockerfile or conda environment file would strengthen reproducibility for the benchmarking pipeline.

---

## Overall recommendation

**Major revision.** The ANI-estimation core is solid and the structural-comparison concept is timely, but the manuscript currently promises more on SV calling than the validation supports, and a critical circular-origin artifact has undermined the flagship case study. The revision should: (a) fix circular-origin handling in `struct`, (b) restructure the structural section around `raw_inverted_fraction` with clear metric definitions, (c) recompute or remove the cagPAI rearrangement associations, and (d) add direct SV validation against a dedicated tool if the word "call" is retained.
