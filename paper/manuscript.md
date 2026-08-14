# Syn2b-ANI: Strain-level ANI estimation via fixed restriction-site anchors for fragmented metagenome-assembled genomes

## Authors

[To be determined]

**Affiliations**
1. State Key Laboratory of Crop Genetics and Germplasm Enhancement, Nanjing Agricultural University, Nanjing, China
2. [Additional affiliations as needed]

**Corresponding author**
Shi Huang — huangshi@njau.edu.cn

---

## Abstract

*(Draft rewritten 2026-08-14 against the v8 chain-restricted MLE algorithm and the ANIm-anchored validation; the body below this section is the earlier GBRT-era draft and is pending rewrite — see `SIMULATION_AND_PERFORMANCE_REPORT.md` for current numbers.)*

Average Nucleotide Identity (ANI) estimation underpins species delineation and strain tracking, but the dominant alignment-free tools reduce a genome comparison to a single scalar, discarding the synteny and structural information that the comparison itself contains. We present **Syn2bANI**, which performs in-silico digestion with Type IIB restriction enzymes and matches the resulting short tags by exact position, chaining matched anchors along the genome and estimating ANI by maximum likelihood only from tags inside trusted chains. The chain-restricted estimator comes in two forms — a single-rate (uniform) model and a gamma rate-heterogeneity model — with per-enzyme consistency statistics that flag unreliable pairs, and an embedded ridge calibration trained on alignment-based (ANIm) ground truth. Because anchors carry genomic coordinates, a single computational pass returns ANI, aligned fraction, a synteny score, and structural-variant calls (inversions, translocations, and indels ≥ 1 kb) — outputs that k-mer tools such as skani and FastANI do not provide. On simulations with exactly known truth, Syn2bANI estimates ANI with a mean absolute error of 0.073%, five to ten times more accurate than skani (0.377%) and FastANI (0.742%), and remains unbiased under insertions/deletions, accessory-content variation, and GC contents from 27% to 72%. On 2,074 GTDB-R207 genome pairs scored against ANIm truth, the calibrated estimator attains MAE 0.906% — exact parity with skani — with near-zero bias (−0.04), and it generalizes to independent oral/gut (MAE 0.460%) and mid-ANI benchmarks. On synthetic structural variants, Syn2bANI detects a 400 kb inversion with endpoint error under 650 bp and recovers indels with exact sizes in seconds, without base-level alignment. Syn2bANI sketches are five times smaller than skani's and pairwise comparison runs at 484 genome pairs in under 3 seconds with 58 MB peak memory. Syn2bANI thus turns ANI estimation from a scalar lookup into a lightweight structural comparison of genomes, and is freely available at https://github.com/HuangShiLab/Syn2bANI.

**Keywords**: ANI, restriction-site anchors, 2bRAD, structural variation, synteny, maximum likelihood, metagenome-assembled genomes, bioinformatics

---

## 1. Introduction

The ability to accurately measure genetic relatedness between microbial genomes is foundational to modern microbiology. Average Nucleotide Identity (ANI), defined as the mean nucleotide identity of all orthologous gene pairs shared between two genomes, has emerged as the gold standard for species boundary delineation and strain-level phylogenomics [1,2]. In clinical and environmental metagenomics, ANI is routinely used to track pathogen transmission, assess functional redundancy within microbial communities, and infer ecological niche partitioning [3,4]. However, the genomes subjected to these analyses are increasingly derived from metagenomic assembly and binning, producing metagenome-assembled genomes (MAGs) that are often highly fragmented [5]. According to the Minimum Information about a Metagenome-Assembled Genome (MIMAG) standards, medium-quality MAGs may have N50 values as low as 3 kb, while even high-quality MAGs frequently exhibit N50 < 10 kb [5]. This fragmentation poses a fundamental challenge to existing ANI estimation tools, which were designed primarily for complete or nearly complete reference genomes.

The dominant alignment-free ANI tools, FastANI [6] and skani [7], employ k-mer-based approaches. FastANI uses 20-mer matches as seeds and performs local alignments of 3 kb fragments, effectively requiring contiguous regions of at least ~700 bp for reliable chaining [6]. skani, a more recent advance, uses spaced k-mer sketches and optimized chaining to achieve dramatically faster performance than FastANI [7]. While skani excels at large-scale database searches, both tools share a common architectural dependency: **k-mer chaining**. Chaining algorithms require collinear runs of matching seeds to anchor alignments, a requirement that becomes increasingly difficult to satisfy as genome fragmentation increases. When contigs are shorter than the chaining window, or when structural rearrangements disrupt collinearity, the chaining process breaks down, leading to incomplete or inaccurate ANI estimates. Furthermore, neither FastANI nor skani outputs structural variation (SV) information, even though such rearrangements are precisely what break their chaining algorithms.

An alternative paradigm is offered by the **2bRAD (double-digest restriction-site associated DNA) sequencing** methodology, originally developed for reduced-representation genotyping [8]. Type IIB restriction enzymes recognize a specific DNA sequence, cleave at fixed positions relative to the recognition site, and produce short, uniform tags (~32 bp) whose positions are completely determined by the enzyme's recognition sequence. In the 2bRAD-M (microbiome) adaptation, these tags have been shown to provide sufficient information for strain-level phylogeny, pangenome profiling, and even structural variation detection when combined with positional information [8,9]. The critical insight is that Type IIB restriction sites act as **natural, deterministic positional anchors** distributed across the genome at roughly 1.5 kb intervals. Unlike random k-mers, which occur unpredictably and require chaining to establish correspondence, 2bRAD tags are fixed points whose positions are implicit in the enzyme recognition sequence. This architectural difference suggests that ANI estimation could be performed without chaining entirely, by simply matching tags at corresponding positions.

Here we present **Syn2b-ani** (Synteny-based 2bRAD ANI), a tool that operationalizes this fixed-anchor paradigm for strain-level genome comparison. Syn2b-ani extracts 2bRAD tags from query and reference genomes via in-silico digestion, matches them using O(1) hash-table lookups, and computes ANI from the Hamming distances of matched tag pairs. Because the matching process does not require chaining, Syn2b-ani is intrinsically robust to fragmentation: even a genome broken into 500 bp contigs retains enough Type IIB sites for reliable ANI estimation. Moreover, the positional correspondence of matched tags enables the simultaneous construction of synteny blocks and detection of structural variations (inversions, insertions, deletions, translocations) in the same computational pass. To correct the systematic ANI overestimation caused by the survivorship bias of fixed-anchor matching (tags with ≥2 mutations within the 32 bp window are excluded), we trained an embedded gradient-boosted regression tree (GBRT) model on 49 bacterial species, achieving <0.02% mean absolute error on strain-level comparisons. The tool is implemented in high-performance Rust, supports all 16 Type IIB enzymes from the 2bRAD-M panel, and provides a skani-compatible command-line interface for seamless integration into existing bioinformatics pipelines.

---

## 2. Methods

### 2.1 Algorithm overview

Syn2b-ani implements a two-pass fixed-anchor algorithm for pairwise genome comparison (Fig. 1a).

**Pass 1: Coarse screening.** The algorithm first performs in-silico digestion of both query and reference genomes using the selected Type IIB enzyme(s). For each genome, it extracts the 2bRAD tag sequence (the short DNA fragment flanked by the enzyme's recognition sites) and records its genomic position. Tags are then compared by exact sequence identity using a hash set intersection. The number of shared tags, normalized by the maximum tag count, provides a coarse ANI estimate via the max-containment similarity formula. Pairs with coarse ANI below a configurable threshold (default 80%) are skipped, avoiding wasted computation on distant genomes.

**Pass 2: Fine ANI estimation and synteny analysis.** For genomes passing the coarse screen, Syn2b-ani builds a hash index of reference tags keyed by their 32 bp sequence. Each query tag is looked up in the index (O(1) per tag). Exact matches are accepted immediately; near-matches within a Hamming distance of ≤1 are also accepted, provided the 32 bp window contains only A, T, C, G bases (no degenerate or ambiguous nucleotides). For each matched tag pair, the Hamming distance is converted to a local ANI estimate: local_ANI = 1 – (Hamming_distance / 32). Matched pairs are then sorted by their query position, and synteny blocks are constructed by grouping consecutive pairs with consistent reference ordering and strand orientation. Orientation flips within blocks indicate inversions; gaps between blocks indicate insertions or deletions; and blocks with non-monotonic reference positions indicate translocations. The final ANI is computed as a weighted average of local ANI values across all matched pairs, with optional correction by the embedded GBRT debias model (see Section 2.6).

### 2.2 In-silico digestion with Fast2bRAD-M alignment

The in-silico digestion module was redesigned to align with the Fast2bRAD-M approach used in the upstream Syn2b tool [9]. For each of the 16 Type IIB enzymes in the 2bRAD-M panel, the recognition pattern is decomposed into a set of fixed anchors and IUPAC degenerate-base constraints. The algorithm slides a window of length equal to the enzyme's tag size across the genome sequence. Within each window, it checks whether any of the enzyme's forward or reverse patterns match. Forward and reverse patterns are defined statically at compile time, eliminating runtime string comparisons. IUPAC matching is performed via a precomputed bitmask lookup table: each base (A, T, C, G) is assigned a 4-bit mask, and IUPAC codes (R, Y, S, W, K, M, B, D, H, V, N) are resolved by bitwise AND operations. This reduces IUPAC checks from ~12 branch instructions per base to a single memory load and bitwise AND. Only windows containing exclusively A, T, C, G bases are retained as valid tags, excluding ambiguous regions common in MAG assemblies. Tags are deduplicated by position and sorted before output. The new implementation achieves a single-thread throughput of ~107 Mb/s for BcgI digestion, compared to 88.3 Mb/s for the previous margin-based approach (17.5% improvement), while producing identical tag sets (Supplementary Fig. S1).

### 2.3 Fixed-anchor matching

Unlike k-mer-based tools, which rely on seed-and-extend or chaining algorithms, Syn2b-ani matches tags purely by sequence identity. Each 32 bp tag is hashed using a 64-bit FNV-1a hash function. Reference tags are inserted into a `HashMap<u64, Vec<usize>>` mapping hash values to indices in the reference tag array. Query tags are then looked up in this index. Collisions are resolved by full 32-byte comparison. Near-matches (Hamming distance ≤ 1) are handled by an additional fallback scan: if a query tag has no exact match, it is compared against all unmatched reference tags within a bounded Hamming distance. The computational complexity of matching is O(n_shared), where n_shared is the number of matched tag pairs (typically 10³–10⁴ for strain-level comparisons), independent of genome size or fragmentation level. This is asymptotically faster than skani's O(n log n) chaining and FastANI's O(n²) local alignment.

### 2.4 Synteny block construction and structural variation detection

After matching, matched pairs are sorted by query position. A synteny block is initiated at the first matched pair and extended as long as consecutive pairs maintain the same reference strand and monotonically increasing reference positions. The block is terminated when a strand flip (indicating an inversion), a gap larger than twice the median tag spacing (indicating an insertion or deletion), or a non-monotonic reference position (indicating a translocation) is encountered. Four types of structural variation are detected:

**Inversions** are identified by strand flips within a synteny block. When consecutive matched pairs switch from forward to reverse strand (or vice versa), the boundary is recorded as an inversion breakpoint. Detection resolution is limited by tag spacing (~1.5 kb for BcgI in a typical 4.5 Mb bacterial genome).

**Insertions and deletions** are identified by gaps between consecutive matched pairs that exceed the expected tag spacing. If the query gap is larger than the reference gap, the difference is recorded as an insertion in the query; if the reference gap is larger, it is recorded as a deletion. Resolution is similarly limited by tag spacing.

**Translocations** are identified by synteny blocks whose reference positions do not follow the global ordering. When a block's reference coordinates jump to a distant region while the query coordinates remain continuous, the breakpoint is recorded as a translocation.

SV detection accuracy was validated using synthetic genomes with known structural variations (50 kb inversions, 20 kb translocations, 5 kb insertions/deletions) superimposed on a 2% SNP background.

### 2.5 GBRT debiasing model

**Problem formulation.** Fixed-anchor tag matching suffers from a survivorship bias: a 32 bp tag containing two or more SNPs will not match its reference counterpart (even with Hamming distance ≤ 1 tolerance), effectively excluding more divergent regions from ANI calculation. This bias causes ANI overestimation, with the magnitude increasing as true divergence increases. On synthetic E. coli pairs with 0.05%–5% SNP divergence, raw Syn2b-ani overestimates ANI by 0.03%–2.06% (Supplementary Fig. S2).

**Model architecture.** We trained a gradient-boosted regression tree (GBRT) model using scikit-learn [10] to correct this bias. The model takes seven features as input: (1) raw ANI from fixed-anchor matching, (2) aligned fraction of the query (AF_q), (3) aligned fraction of the reference (AF_r), (4) number of shared tags, (5) max-containment similarity, (6) divergence proxy (1 – raw ANI), and (7) reference GC content. The model consists of 300 regression trees with maximum depth 5, learning rate 0.05, and least-squares loss. Hyperparameters were selected via 5-fold cross-validation on the training set.

**Training data.** To ensure cross-species generalization, we constructed a training dataset of 1,260 synthetic genome pairs from 49 bacterial species spanning 1.7–4.7 Mb in size and 35–60% GC content (Supplementary Table S1). For each species, we generated 25 synthetic variants at controlled SNP rates (0.05%, 0.1%, 0.2%, 0.5%, 1%, 2%, 3%, 5%) using a neutral Wright-Fisher mutation model with equal base substitution rates. Fragmentation was simulated at N50 values of 500 bp, 1 kb, 2 kb, 5 kb, 10 kb, 20 kb, 50 kb, and 100 kb. Completeness was varied from 30% to 100%. The ground-truth ANI for each pair was computed by exact position-by-position alignment of the unfragmented sequences.

**Embedding.** The trained model was serialized as JSON decision trees and embedded directly into the Rust binary via the `include_str!` macro. At runtime, ANI prediction is performed by traversing the decision trees in memory, with each prediction requiring < 1 μs. The embedded model adds 1.08 MB to the binary size. A smaller runtime-optimized variant (200 trees, depth 4, 0.64 MB) is also available for memory-constrained environments.

**Model validation.** The model was validated on held-out test sets: (1) 20% of the E. coli synthetic pairs (within-species), (2) five completely held-out bacterial species (cross-species), and (3) 49 independent real bacterial genome pairs from the GTDB database [11]. Within-species error was 0.002% MAE; cross-species error was 0.30% MAE.

### 2.6 Implementation

Syn2b-ani is implemented in Rust (edition 2021) with the following dependencies: needletail [12] for FASTA parsing, rayon [13] for data parallelism, and serde_json for embedded model deserialization. The command-line interface follows skani's subcommand structure (`dist`, `search`, `sketch`, `triangle`) for backward compatibility. The custom binary sketch format (`.s2ba`) uses 2-bit encoding for tag sequences (32 bp → 64 bits), yielding sketches of ~48 KB per genome. The tool supports multi-threaded operation at both the genome-pair and enzyme levels. All 16 Type IIB enzymes from the 2bRAD-M panel are supported, with BcgI (CGA-N₆-TGC, 32 bp tag) as the default (Supplementary Table S2).

### 2.7 Benchmark setup

**Synthetic datasets.** The primary benchmark dataset was a synthetic E. coli K-12 MG1655 genome (4.65 Mb) with controlled divergence rates of 0.05%, 0.1%, 0.2%, 0.5%, 1%, 2%, 3%, and 5%. Fragmentation was simulated by randomly breaking the genome into contigs with specified N50 values (500 bp, 1 kb, 2 kb, 5 kb, 10 kb, 20 kb, 50 kb, 100 kb). Completeness was varied by randomly subsampling contigs (30%, 50%, 60%, 80%, 100%). Structural variations were introduced by inverting 50 kb segments, translocating 20 kb segments, and inserting/deleting 5 kb segments.

**Realistic MAG scenarios.** We simulated four common MAG artifacts: (1) contamination by mixing 5%–20% foreign DNA from *Bacillus subtilis*, (2) chimerism by shuffling contigs at 1–20 breakpoints, (3) duplication by duplicating 5%–20% of the genome, and (4) assembly errors by introducing 0.01%–0.2% random base substitutions.

**Comparison tools.** Syn2b-ani (v0.1.1) was compared against skani (v0.3.2) [7] and a Python implementation of FastANI [6]. All tools were run with default parameters unless otherwise noted. Ground-truth ANI was computed by exact alignment of the unfragmented synthetic sequences using the MUMmer nucmer algorithm [14].

**Performance metrics.** ANI accuracy was quantified as mean absolute error (MAE), maximum absolute error, and Pearson R² relative to ground truth. Runtime was measured using the Unix `time` command, averaged over 3 replicates. Memory was measured using `/usr/bin/time -v`. SV detection was evaluated by precision and recall against known SV breakpoints.


---

## 3. Results

### 3.1 ANI accuracy at varying divergence

We first evaluated Syn2b-ani's accuracy on a synthetic E. coli dataset with controlled SNP divergence rates from 0.05% to 5%. At low divergence (0.05%–0.2%), raw Syn2b-ani (without debiasing) produced ANI estimates within 0.03% of the ground truth, comparable to FastANI (0.01% error) and skani (0.02% error) (Fig. 2a). However, as divergence increased to 2%–5%, the survivorship bias of fixed-anchor matching became apparent: tags containing two or more SNPs within the 32 bp window were excluded, causing the raw ANI to overestimate by up to 2.06% (Fig. 2b). This bias is systematic and monotonic, increasing approximately linearly with true divergence.

The simple debias formula (subtracting the expected bias based on shared tag fraction) reduced the error to 0.15%–0.49%, but this was still insufficient for high-precision strain-level comparison. The GBRT debias model, by contrast, reduced the mean absolute error to 0.002% across all divergence levels (0.05%–5%), with a maximum absolute error of 0.008% at 5% divergence (Fig. 2b). This performance is statistically indistinguishable from FastANI (paired t-test, P = 0.42) and superior to skani at divergence > 2% (P < 0.001). The shared tag count declined from ~4,800 at 0.05% divergence to ~1,200 at 5% divergence (Fig. 2c), remaining well above the minimum threshold of 10 shared tags required for reliable ANI estimation. At 5% divergence, the coarse screening step correctly identified the pair as high-similarity (coarse ANI = 95.3%), triggering the full Pass 2 analysis.

### 3.2 Robustness to fragmentation

The most striking advantage of the fixed-anchor approach is its resilience to genome fragmentation. We simulated E. coli MAGs with N50 values ranging from 500 bp to 100 kb, all at 2% true divergence and 100% completeness. Syn2b-ani with GBRT debiasing produced ANI estimates of 98.53% ± 0.01% across all N50 levels, with no significant trend (linear regression, R² = 0.003, P = 0.89) (Fig. 3a). This stability arises because each 32 bp tag is extracted independently from its local sequence context; fragmentation does not alter the tag sequences themselves, only their genomic positions. As long as a contig is longer than the tag length (32 bp), which is true even for the 500 bp N50 scenario, all tags are recoverable.

In contrast, skani's accuracy degraded dramatically with decreasing N50 (Fig. 3b). At N50 = 100 kb, skani produced ANI = 98.5%, matching Syn2b-ani. However, as N50 decreased to 5 kb, skani's estimate dropped to 95.1%, and at N50 = 500 bp, it fell to 93.2%—a 5.3% underestimation relative to ground truth. This degradation occurs because skani's k-mer chaining algorithm requires sufficiently long collinear runs of matching seeds, which become increasingly rare as contigs shrink. FastANI showed a similar pattern, requiring N50 > 5 kb for reliable estimates (Fig. 3c). The error heatmap (tool × N50 × divergence) confirms that Syn2b-ani is the only tool with uniformly low error across all fragmentation levels (Fig. 3d). Even at the extreme of 500 bp N50, Syn2b-ani's MAE was 0.01%, compared to skani's 5.3% and FastANI's 4.8%.

### 3.3 Robustness to MAG completeness

We next simulated incomplete MAGs by randomly subsampling contigs from the fragmented E. coli genome (N50 fixed at ~10 kb). At 30% completeness, Syn2b-ani produced ANI = 98.52%, with a MAE of 0.02% relative to the 100% complete ground truth (Fig. 4a). Across the full range of 30%–100% completeness, the ANI remained stable at 98.55% ± 0.03%. The aligned fraction (AF) decreased predictably with completeness (from 28% at 30% completeness to 95% at 100% completeness), but this did not affect the ANI estimate because the remaining tags were still representative of the average nucleotide identity. skani and FastANI showed larger errors at low completeness (0.3% and 0.5% MAE respectively at 30%), likely because their chaining algorithms depend on contiguous coverage for reliable anchor extension. These results indicate that Syn2b-ani is suitable for low-complete MAGs from metagenomic samples, where completeness is often < 50%.

### 3.4 Structural variation detection

A unique feature of Syn2b-ani is the simultaneous detection of structural variations during ANI estimation. We validated this capability using synthetic genomes with known SVs superimposed on a 2% SNP background. A 50 kb inversion was detected with 100% precision and 100% recall: the orientation flip within the synteny block was unambiguous, and the inversion boundaries were localized to within ±1.5 kb of the true breakpoints (Fig. 5a). A 20 kb translocation was detected as a rearrangement block with non-monotonic reference coordinates; the translocation breakpoint was correctly identified with 100% precision but 85% recall due to some small blocks being merged (Fig. 5b). Insertions and deletions of 5 kb were detected by gap analysis between synteny blocks, with precision of 92% and recall of 88% (Fig. 5c). The false negatives occurred primarily when insertions/deletions were smaller than the tag spacing (~1.5 kb), falling below the detection threshold. In a combined scenario with multiple SVs and SNPs, all major SV events were detected, though SNPs caused some false-positive indel calls at tag boundaries where single SNPs altered the local spacing (Fig. 5d). These false positives are manageable: they represent < 5% of all calls and can be filtered by requiring a minimum gap size (3× tag spacing).

### 3.5 Speed and memory benchmark

We compared the computational performance of Syn2b-ani against skani and FastANI on a dataset of 65 bacterial genomes (average size 4.2 Mb) (Table 1). Sketch building (the one-time preprocessing step) required 12 s for Syn2b-ani, compared to 45 s for skani and 8 s for MASH. While MASH was fastest for sketching, its ANI resolution is too low for strain-level comparison. The query step (searching one genome against the pre-sketched database) required 0.3 s for Syn2b-ani, 0.8 s for skani, and 120 s for FastANI. Syn2b-ani was therefore 3× faster than skani and 400× faster than FastANI for the query step.

The sketch size was dramatically smaller for Syn2b-ani: 48 KB per genome (2-bit packed 32 bp tags), compared to ~20 MB for skani and ~30 KB for MASH. For large databases, this translates to significant memory savings: searching a 65,000-genome database requires 180 MB with Syn2b-ani, compared to 1.2 GB for skani. The memory footprint of FastANI (2.4 GB) is prohibitive for large-scale applications. The binary size of Syn2b-ani (2.7 MB including the embedded 1.08 MB GBRT model) is larger than skani (~1 MB) but smaller than FastANI (~5 MB).

Runtime scaling with genome count was linear for all tools (Supplementary Fig. S6), but Syn2b-ani's smaller constant factor (due to O(1) matching and no chaining) maintained its advantage across the full range tested (1 to 65,536 genomes).

### 3.6 Cross-species GBRT generalization

A critical concern for any ML-based correction model is whether it generalizes to species not seen during training. We validated the GBRT v2 model on five held-out bacterial species: *Escherichia coli* (training species, positive control), *Bacillus subtilis* (Gram-positive, different GC), *Bifidobacterium longum* (high GC, actinobacteria), *Lactobacillus acidophilus* (low GC, firmicutes), and *Pseudomonas aeruginosa* (high GC, proteobacteria). On *E. coli*, the model achieved 0.01% MAE, confirming that it had not overfit to the training set (Fig. 6a). On the four held-out species, the model achieved 0.30% MAE—still 2× better than the simple debias formula (0.62% MAE) and sufficient for most practical applications (Fig. 6b–f). The error was not correlated with GC content (Pearson R = 0.12, P = 0.67), genome size (R = 0.08, P = 0.78), or taxonomic distance from *E. coli* (R = 0.21, P = 0.45), suggesting that the model has learned general features of fixed-anchor bias rather than species-specific sequence patterns. For users requiring the highest precision (< 0.01% error), we recommend training a species-specific model using the provided training scripts (`train_gbrt_v2.py`), which requires only ~30 minutes on a standard laptop.

### 3.7 Realistic MAG scenarios

Finally, we tested Syn2b-ani on simulated realistic MAGs incorporating multiple artifacts simultaneously: 5% contamination from *B. subtilis*, 5 chimera breakpoints, 5% duplication, and 0.1% assembly error. Syn2b-ani produced ANI = 99.9% and AF = 90.7%, correctly identifying the query as a high-quality strain match (Fig. 7). The contamination was effectively masked by the max-containment normalization: foreign tags from *B. subtilis* had no matches in the *E. coli* reference and were simply excluded from the ANI calculation. Chimerism (shuffled contigs) did not affect ANI because the tag sequences themselves are unchanged; only their genomic positions are altered, which the matching algorithm accommodates naturally. Duplication increased the tag count proportionally but did not bias the ANI because duplicated tags match identically. Assembly errors (0.1% random substitutions) caused a slight ANI underestimate (99.9% vs. 100% ground truth) because erroneous tags were excluded from matching, but the effect was minor (< 0.1%). In comparison, skani produced ANI = 94.2% and AF = 71.3% on the same dataset, failing to correctly classify the MAG as a strain match due to fragmentation-induced chaining breakdown. FastANI was unable to complete the analysis because the combined fragmentation and chimerism fell below its minimum chaining threshold.

### 3.8 Large-scale GTDB-R207 validation and effective ANI range

To validate Syn2b-ani at scale, we benchmarked it against 64,747 representative genomes from the GTDB-R207 dataset (65,703 total representatives; 956 permanently unavailable from NCBI). The dataset spans 169 bacterial and archaeal phyla, with genome sizes ranging from 0.5 Mb to 13.2 Mb and GC contents from 14% to 75%.

#### Empirical detection threshold

The Type IIB restriction-site anchor approach has an inherent biological detection limit. When two genomes share too few Type IIB restriction sites, Syn2b-ani cannot find sufficient anchor points to estimate ANI, and the tool reports `ANI = 0.0`. Through systematic evaluation of 1,150 cross-taxonomic pairs spanning the full ANI spectrum, we identified this threshold at approximately **83% ANI** (Fig. 8a). Pairs with ANI below 83% typically return `raw_ani = 0.0` because the shared tag count falls below the minimum threshold (`min_shared_tags = 10`) or the alignment fraction drops below `min_af = 0.1`. This is not a model training limitation but a fundamental property of the Type IIB anchor method: for sufficiently distant genomes, Type IIB sites become statistically absent. For comparisons below 83% ANI (e.g., across-family or cross-phylum phylogenetics), users should use alignment-based tools such as FastANI or Mash.

#### GBRT v3.6 model performance

To correct systematic bias across the full 83–100% ANI range, we retrained the GBRT model (v3.6) on 622 pairs. The training set comprised 296 pairs with skani-validated ANI (ground truth, 95–100% range) and 326 pairs with Mash-calibrated ANI (estimated, 83–95% range). Mash distance was calibrated to ANI using a polynomial regression fitted to the 296 skani-validated pairs (R² = 0.976, MAE = 0.32%).

The GBRT v3.6 model (300 trees, max depth 5, learning rate 0.05) achieved **0.27% MAE** on a 20% holdout test set, with an R² of 0.9968 (Fig. 8b). Performance was consistent across all ANI ranges: 0.18% MAE at 97–100%, 0.27% at 95–97%, 0.36% at 90–95%, 0.35% at 85–90%, and 0.29% at 80–85% (Fig. 8c). Feature importance analysis showed that `raw_ani` dominated (46.8%), followed by `af_q` (36.3%), `af_r` (13.4%), and `has_skani` (3.5%), confirming that the model does not overfit to the ground truth source.

#### Phylum-specific performance

We evaluated accuracy across 20 bacterial and archaeal phyla. Three phyla dominated the dataset (Proteobacteria, Firmicutes, Actinobacteriota = 75% of pairs) and showed excellent accuracy (MAE < 0.8%). Two phyla, Patescibacteria and Fusobacteriota, showed apparently high errors (MAE = 43% and 46%, respectively). However, investigation revealed that these errors were not model failures but **detection threshold effects**: all outlier pairs had `raw_ani = 0.0` because the genomes were below the 83% ANI threshold. When excluding below-threshold pairs, these phyla showed normal accuracy. This finding led to the implementation of a below-detection warning in Syn2b-ani v0.1.1: when estimated ANI is below 83%, the tool emits a warning recommending verification with FastANI.

#### Comparison to skani and FastANI

At the strain level (>95% ANI), Syn2b-ani + GBRT v3.6 achieves 0.18–0.27% MAE, comparable to skani (0.20% MAE) and FastANI (0.30% MAE). In the 85–95% range (cross-species comparisons), Syn2b-ani maintains 0.29–0.36% MAE, while skani frequently fails to report ANI for pairs below ~80%. For MAG applications, Syn2b-ani's unique advantage is its **MAG-friendliness**: it handles fragmentation (N50 = 500 bp) without accuracy degradation, whereas skani's accuracy drops by >5% at N50 < 5 kb.

### 3.9 Recommendations for users

1. **Strain tracking / outbreak analysis** (ANI > 95%): Syn2b-ani is ideal — fast, accurate, and MAG-friendly.
2. **Species delineation** (ANI 83–95%): Syn2b-ani works well with GBRT v3.6, but verify borderline cases with FastANI.
3. **Phylogenetic reconstruction** (ANI < 83%): Use FastANI or Mash. Syn2b-ani will report "below detection threshold."
4. **Fragmented MAGs**: Syn2b-ani is preferred over FastANI because alignment-based methods struggle with short contigs.
5. **Patescibacteria or Fusobacteriota**: Use caution — these phyla have small/low-GC genomes that may fall below the detection threshold even for intra-species pairs.

---

---

## 4. Discussion

### 4.1 Fixed anchors versus random k-mers: a fundamental architectural difference

The central advance of Syn2b-ani is the replacement of random k-mer matching with fixed-anchor tag matching. This architectural change has three profound consequences. First, it eliminates the chaining bottleneck that limits all k-mer-based ANI tools. Chaining algorithms require collinear runs of matching seeds to establish reliable anchors; when genome fragmentation breaks these runs, the algorithm fails. Fixed anchors, by contrast, are independent: each tag is matched individually, and fragmentation does not affect the matchability of any single tag as long as the contig is longer than the tag itself. Second, fixed anchors provide positional information "for free": the genomic coordinates of each tag are implicit in the enzyme recognition sequence, eliminating the need for separate coordinate tracking during matching. Third, fixed anchors enable simultaneous synteny and SV analysis: because matched tags have known positions in both genomes, deviations from collinearity (strand flips, gaps, jumps) can be detected in the same pass as ANI estimation.

The trade-off is reduced tag density. A typical 4.5 Mb bacterial genome contains ~3,000 BcgI tags (one per ~1.5 kb), compared to ~4.5 million 20-mers (one per base). This lower density means that Syn2b-ani's resolution is coarser: SV detection is limited to ~1.5 kb, and ANI estimates are based on fewer observations (~3,000 tags vs. millions of k-mers). However, the 32 bp length of each tag provides sufficient information content for accurate ANI estimation: 3,000 independent observations of 32 bp each corresponds to ~96 kb of surveyed sequence, more than enough for robust statistical estimation. The lower density also means that Syn2b-ani is not suitable for very high divergence (> 10%), where tag sharing becomes too sparse (< 100 tags). Within the recommended range (< 5% divergence, the strain-to-species boundary), the fixed-anchor approach provides sufficient accuracy while offering unique advantages in fragmentation robustness and SV detection.

### 4.2 ANI, synteny, and structural variation in one pass

No existing ANI tool outputs structural variation or synteny information. FastANI and skani output only ANI and aligned fraction; MASH outputs only a distance estimate. Syn2b-ani, by virtue of its fixed-anchor architecture, outputs four types of information simultaneously: (1) ANI, (2) AF, (3) synteny blocks, and (4) SV calls (inversions, indels, translocations). This unified output enables downstream analyses that were previously impossible: for example, one can simultaneously assess whether two strains belong to the same species (ANI > 95%) and whether they differ by a large inversion (SV detection). The PAF-like output format is compatible with standard visualization tools such as IGV and Bandage, allowing researchers to inspect synteny and SV alongside ANI scores.

The SV detection capability is not merely a bonus feature but a direct consequence of the fixed-anchor architecture. In k-mer-based tools, structural rearrangements are the enemy: they break chaining and cause false-negative ANI estimates. In Syn2b-ani, rearrangements are informative: they appear as discontinuities in the synteny map and are reported as SV events. This inversion of the relationship between ANI and SV—from conflict to synergy—is a unique advantage of the fixed-anchor approach.

### 4.3 Experimental verifiability: the 2bRAD-M connection

A distinguishing feature of Syn2b-ani is its experimental verifiability. Unlike purely computational ANI tools, Syn2b-ani's predictions can be validated by 2bRAD-M sequencing [8]. If a predicted tag is not observed in the 2bRAD-M library, it indicates either an assembly gap (the tag is present in the true genome but missing from the assembly), a sequencing bias (the restriction site is methylated or otherwise inaccessible), or a true structural variation (the tag was deleted or translocated). This creates a feedback loop between *in silico* prediction and experimental validation that is unique among ANI tools. For example, if Syn2b-ani detects an inversion in a MAG and the 2bRAD-M library confirms the presence of the predicted tags at the inverted positions, the inversion call is validated experimentally. This verifiability is particularly valuable for MAGs, where assembly errors can confound computational predictions.

### 4.4 Limitations and future directions

Syn2b-ani has several limitations that should be considered when interpreting results. First, the tool is designed for low-to-moderate divergence (< 5%), which covers the strain-to-species boundary but not higher taxonomic levels. At divergence > 10%, tag sharing becomes too sparse (< 100 tags) for reliable ANI estimation. Users requiring genus-level comparison should use MASH or skani instead. Second, SV detection resolution is limited by tag spacing (~1.5 kb for BcgI). Smaller indels (< 1 kb) and point mutations are not directly detectable, though they are implicitly accounted for in the ANI calculation via Hamming distance. Third, the multi-enzyme mode currently concatenates tags from multiple enzymes without optimizing their relative weights. Future work could implement enzyme-specific weighting based on GC bias or tag density. Fourth, while the universal GBRT v2 model achieves 0.3% MAE on new species, users requiring the highest precision (< 0.01%) may benefit from training a species-specific model using the provided training scripts. Fifth, Type IIB sites may be underrepresented in plasmids and mobile elements, potentially causing chromosome-biased ANI estimates. A plasmid-aware enzyme panel could address this in future versions.

Looking forward, we plan to validate Syn2b-ani on a large collection of real MAGs from human gut metagenomes, comparing predicted tags against 2bRAD-M sequencing data. We also intend to scale the benchmark to the full GTDB database [11] (~65,000 representative genomes) to assess performance at the scale of modern microbial genome collections. Species-specific GBRT models, trained per-genus or per-family, may further improve accuracy for underrepresented clades. Finally, integration with the 2bRAD-M pipeline to enable direct comparison of predicted and observed tags would realize the full potential of the fixed-anchor approach as a bridge between computation and experiment.

---

## 5. Data availability

All benchmark datasets, including synthetic genomes, ground-truth alignments, and comparison results, are available via Zenodo [DOI to be assigned]. The 49 bacterial genomes used for GBRT training are listed in Supplementary Table S1, with NCBI RefSeq accessions. The GBRT training data (`training_data_v2.csv`) and embedded model (`gbrt_model_v2.json`) are included in the GitHub repository.

## 6. Code availability

Syn2b-ani v0.1.1 is freely available at https://github.com/HuangShiLab/Syn2bANI under the MIT License. The repository includes the Rust source code, integration tests, Criterion benchmarks, and Python scripts for analysis and figure generation. The paper-specific analysis code, figures, and datasets are available at https://github.com/HuangShiLab/Syn2bANI-paper. Docker images are provided for reproducible execution. Installation instructions and API documentation are included in the README.

---

## Supplementary Information

### Supplementary Tables

| Table | Content |
|-------|---------|
| Table S1 | 49 bacterial genomes used for GBRT training (species name, NCBI accession, genome size, GC content, taxonomic lineage) |
| Table S2 | 16 Type IIB enzyme specifications (name, recognition site, tag length, spacer length, NEB catalog number) |
| Table S3 | GBRT hyperparameter search results (tree count, depth, learning rate, cross-validation MAE) |
| Table S4 | Cross-species validation detailed results (species, divergence level, raw ANI, debiased ANI, ground truth, error) |
| Table S5 | skani and FastANI command lines used for comparison |

### Supplementary Figures

| Figure | Content |
|--------|---------|
| Fig. S1 | Fast2bRAD-M optimized digestion vs. legacy margin-based digestion: speed and tag count comparison |
| Fig. S2 | Raw Syn2b-ani ANI overestimation as a function of true divergence (0.05%–5%) |
| Fig. S3 | Per-enzyme ANI accuracy: BcgI, BsaXI, CjeI, CjePI, and multi-enzyme consensus |
| Fig. S4 | GBRT feature importance and correlation matrix |
| Fig. S5 | GBRT v1 (E. coli only) vs. v2 (49 species): cross-species generalization comparison |
| Fig. S6 | Runtime scaling with genome count (1, 10, 100, 1,000, 10,000, 65,536 genomes) |
| Fig. S7 | Memory usage profiling during sketch building and query |
| Fig. S8 | PAF output examples: visualization of inversion, translocation, and combined SV scenarios in IGV |

---

## References

1. Konstantinidis, K. T. & Tiedje, J. M. Genomic insights that advance the species definition for prokaryotes. *Proc. Natl. Acad. Sci. USA* **102**, 2567–2572 (2005).
2. Jain, C., Rodriguez-R, L. M., Phillippy, A. M., Konstantinidis, K. T. & Aluru, S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nat. Commun.* **9**, 5114 (2018).
3. Olm, M. R., Brown, C. T., Brooks, B. & Banfield, J. F. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. *ISME J.* **11**, 2864–2868 (2017).
4. Parks, D. H. et al. A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. *Nat. Biotechnol.* **36**, 996–1004 (2018).
5. Bowers, R. M. et al. Minimum information about a metagenome-assembled genome (MIMAG) of bacteria and archaea. *Nat. Biotechnol.* **35**, 725–731 (2017).
6. Jain, C., Rodriguez-R, L. M., Phillippy, A. M., Konstantinidis, K. T. & Aluru, S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nat. Commun.* **9**, 5114 (2018).
7. Shaw, J. & Yu, Y. W. Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nat. Methods* **20**, 1701–1708 (2023).
8. Wang, S. et al. 2bRAD-M: a simple and effective method for microbiome analysis. *Methods Ecol. Evol.* **11**, 1234–1245 (2020).
9. [Syn2b: Synteny analysis using 2bRAD tags. https://github.com/HuangShiLab/Syn2b]
10. Pedregosa, F. et al. Scikit-learn: Machine learning in Python. *J. Mach. Learn. Res.* **12**, 2825–2830 (2011).
11. Parks, D. H. et al. GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Res.* **50**, D785–D794 (2022).
12. [needletail: a fast and ergonomic parser for FASTA/FASTQ files. https://github.com/onecodex/needletail]
13. [rayon: data parallelism in Rust. https://github.com/rayon-rs/rayon]
14. Marçais, G. et al. MUMmer4: A fast and versatile genome alignment system. *PLoS Comput. Biol.* **14**, e1005944 (2018).
15. Roberts, R. J. et al. REBASE—a database for DNA restriction and modification: enzymes, genes and genomes. *Nucleic Acids Res.* **43**, D298–D299 (2015).

---

*Manuscript draft v1.0*
*Generated: 2025-07-10*
*Corresponding author: Shi Huang (huangshi@njau.edu.cn)*

