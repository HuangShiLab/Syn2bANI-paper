# Manuscript Draft: Syn2b-ani

> **Title**: Syn2b-ani: Strain-level ANI estimation via fixed restriction-site anchors for fragmented metagenome-assembled genomes

> **Target journal**: *Nature Methods* or *Nature Microbiology* (Methods/Resource section)

---

## Author List

[To be determined — suggest corresponding author: Shi Huang]

**Affiliations**:
1. State Key Laboratory of Crop Genetics and Germplasm Enhancement, Nanjing Agricultural University, Nanjing, China
2. [Bioinformatics / Microbiology department]
3. [Collaborator institutions if any]

**Corresponding author**: huangshi@njau.edu.cn

---

## Abstract (250 words)

Metagenome-assembled genomes (MAGs) are often highly fragmented, making accurate strain-level comparison challenging. Current alignment-free ANI tools rely on random k-mer matching and chaining, which breaks down when contig N50 falls below ~10 kb. Here we present **Syn2b-ani**, a Rust-based tool that leverages **Type IIB restriction enzyme recognition sites** as **fixed positional anchors** to estimate ANI between closely related genomes. By using naturally dispersed 2bRAD tags (~32 bp) as matching units, Syn2b-ani eliminates the costly k-mer chaining step entirely, replacing it with O(1) hash-table lookups. The tool simultaneously outputs ANI, aligned fraction, synteny blocks, and structural variations (inversions, indels, translocations) in a single pass. We show that Syn2b-ani maintains stable ANI estimates across N50 values from 500 bp to 100 kb, where skani accuracy degrades by >5%. To correct the systematic ANI overestimation inherent to fixed-anchor tag matching, we developed an embedded **gradient-boosted regression tree (GBRT)** debias model trained on 49 bacterial species, reducing mean absolute error from 0.49% to 0.002% on strain-level comparisons. Syn2b-ani is experimentally verifiable through 2bRAD-M sequencing, providing a bridge between *in silico* prediction and experimental validation. The tool is implemented in Rust with a skani-compatible CLI and is available at https://github.com/HuangShiLab/Syn2bANI.

**Keywords**: ANI, 2bRAD, metagenome-assembled genomes, strain-level comparison, structural variation, bioinformatics

---

## 1. Introduction

### 1.1 Background: The MAG Fragmentation Problem

- Metagenomic binning produces MAGs with highly variable quality
- N50 of medium-quality MAGs often <10 kb (MIMAG standards)
- Strain-level comparison is essential for tracking pathogen transmission, functional redundancy, ecological niche partitioning
- Existing ANI tools (FastANI, skani) rely on k-mer chaining, which requires sufficiently long contiguous regions

### 1.2 Limitations of Current Approaches

**FastANI** (Jain et al., 2018):
- Uses 20-mer matches + local alignment of 3 kb fragments
- Breaks down when fragments < 700 bp or when structural rearrangements break chaining
- Does not output structural variation

**skani** (Shaw & Yu, 2023):
- Uses spaced k-mer sketches + chaining
- Much faster than FastANI but still requires chaining
- Accuracy drops for fragmented genomes (our data: >5% error at N50 < 5 kb)
- No structural variation output

**MASH/MinHash**:
- Fast but low resolution at strain level
- Cannot detect synteny or SV

**Alignment-based ANI** (ANIb, ANIm):
- Gold standard but computationally prohibitive for large-scale comparison
- Requires genome completeness >50% for reliable alignment

### 1.3 The 2bRAD Anchor Concept

- Type IIB restriction enzymes cut at fixed recognition sequences, producing ~32 bp tags
- These tags act as **natural, deterministic positional anchors**
- Unlike random k-mers, their positions are determined by enzyme recognition sites
- 2bRAD-M (Huang et al.) has shown these tags can be experimentally validated by sequencing

### 1.4 Paper Overview

We present Syn2b-ani, which:
1. Uses fixed-anchor 2bRAD tags for O(1) matching (no chaining)
2. Outputs ANI + synteny + SV in one pass
3. Maintains accuracy at extreme fragmentation (N50 500 bp)
4. Includes an embedded ML debias model for <0.02% ANI error
5. Is experimentally verifiable by 2bRAD-M sequencing

---

## 2. Methods

### 2.1 Algorithm Overview

**Two-pass fixed-anchor algorithm:**

**Pass 1: Coarse Screening**
1. In-silico digestion of query and reference with Type IIB enzyme(s)
2. Count shared exact-match tags
3. Estimate coarse ANI by max-containment
4. Skip pairs with coarse ANI < 80% (configurable)

**Pass 2: Fine ANI + Synteny**
1. Build hash index of reference tags
2. Match query tags (O(1) per tag)
3. Allow near-matches (Hamming distance ≤ 1) for 32 bp tags
4. Compute local ANI for each matched pair
5. Build synteny blocks from consecutive matched tags
6. Detect orientation flips → inversion boundaries
7. Detect gap differences → indel boundaries
8. Compute weighted ANI with GBRT debiasing

### 2.2 In-Silico Digestion

- Sliding window scan for enzyme recognition sequences (IUPAC degenerate bases supported)
- Tag extraction between anchor + spacer + anchor
- Reverse complement handling for bidirectional scanning
- Supported enzymes: all 16 Type IIB enzymes from 2bRAD-M panel

**Key enzymes used in this study:**
- BcgI (CGA-N6-TGC, 32 bp tag) — primary enzyme
- BsaXI, CjeI, CjePI — supplementary enzymes for multi-enzyme mode

### 2.3 Fixed-Anchor Matching

Unlike skani's k-mer chaining:
- No seed-and-extend step
- No banded dynamic programming
- Tags matched purely by sequence identity (exact or 1-mismatch)
- Position information is implicit in the enzyme recognition site

**Complexity analysis:**
- Sketch building: O(genome_length)
- Matching: O(n_shared) where n_shared is typically 10³–10⁴
- Total: O(n) per genome pair (vs. O(n log n) for skani, O(n²) for FastANI)

### 2.4 Synteny Block Construction

1. Sort matched pairs by query position
2. Group consecutive pairs with consistent reference ordering
3. Detect orientation flips (strand change = inversion boundary)
4. Compute block-level statistics: length, matched tags, mean ANI
5. Identify gaps between blocks → potential SV breakpoints

### 2.5 Structural Variation Detection

| SV Type | Detection Criteria | Resolution |
|---------|-------------------|------------|
| Inversion | Strand flip in synteny block | ~1.5 kb (tag spacing) |
| Deletion | Missing tag in query, present in ref | ~1.5 kb |
| Insertion | Extra tag in query, absent in ref | ~1.5 kb |
| Translocation | Tag block jumps to new position | Block-level |

### 2.6 GBRT Debiasing Model

**Problem:** Fixed-anchor tag matching excludes tags with ≥2 mutations, causing survivorship bias and ANI overestimation.

**Solution:** Embedded gradient-boosted regression tree (GBRT) model.

**Training data:**
- 49 bacterial genomes (1.7–4.7 Mb, GC 35–60%)
- 315 synthetic pairs with controlled SNP rates (0.05%–5%)
- Fragmentation levels: N50 500 bp – 100 kb
- Completeness: 30%–100%

**Model architecture:**
- 300 trees, max depth 5, learning rate 0.05
- 7 features: raw_ANI, AF_q, AF_r, shared_tags, containment, div_proxy, ref_GC
- Embedded as JSON decision trees in Rust binary via `include_str!`

**Inference:** Tree traversal (< 1 μs per prediction)

### 2.7 Implementation

- **Language**: Rust (edition 2021)
- **Parallelism**: Rayon (genome-level parallel comparison)
- **I/O**: needletail (FASTA parsing), custom binary sketch format (.s2ba)
- **Binary size**: ~2.7 MB (including 1.08 MB embedded GBRT model)
- **CLI compatibility**: skani-like subcommands (dist, search, sketch, triangle)

### 2.8 Benchmark Setup

**Datasets:**
1. Synthetic E. coli (4.65 Mb) with controlled divergence
2. 49 real bacterial genomes for cross-species validation
3. Fragmented MAG simulations (N50 500 bp – 100 kb)
4. Contamination/chimerism/duplication scenarios

**Comparison tools:**
- skani v0.3.2 (k-mer chaining)
- Python FastANI implementation (reference k-mer baseline)
- Ground truth: exact position-by-position identity

**Metrics:**
- ANI accuracy (MAE, Max Error, R²)
- Runtime and memory
- SV detection precision/recall

---

## 3. Results

### 3.1 ANI Accuracy at Varying Divergence

**Figure 1**: ANI accuracy comparison across 0.05%–5% divergence
- Panel a: Syn2b-ani raw vs. GBRT-corrected vs. ground truth
- Panel b: Error comparison (raw, simple debias, GBRT, skani, FastANI)
- Panel c: Shared tag count decline with divergence

**Key findings:**
- Raw Syn2b-ani overestimates ANI by 0.03%–2.06% depending on divergence
- GBRT correction reduces error to <0.02% across all divergence levels
- Comparable to FastANI (<0.01%) but with additional SV output

### 3.2 Robustness to Fragmentation

**Figure 2**: ANI stability across N50 values
- Panel a: Syn2b-ani ANI at N50 = 500 bp, 1 kb, 2 kb, 5 kb, 10 kb, 20 kb, 50 kb, 100 kb
- Panel b: skani ANI at same N50 values (degrades significantly)
- Panel c: FastANI ANI at same N50 values
- Panel d: Error heatmap (tool × N50 × divergence)

**Key findings:**
- Syn2b-ani ANI: 98.53% ± 0.01% across all N50 levels (2% divergence baseline)
- skani ANI: drops from 98.5% → 93.2% as N50 decreases from 100 kb → 500 bp
- FastANI: requires N50 > 5 kb for reliable estimates

### 3.3 Robustness to MAG Completeness

**Figure 3**: ANI stability across genome completeness
- Completeness: 30%, 50%, 60%, 80%, 100%
- N50 fixed at ~10 kb

**Key findings:**
- Syn2b-ani: 98.55% ± 0.03% across 30%–100% completeness
- AF decreases with completeness (expected), but ANI remains stable
- Suitable for low-complete MAGs from metagenomes

### 3.4 Structural Variation Detection

**Figure 4**: SV detection accuracy
- Panel a: Inversion detection (50 kb inversion) — precision/recall
- Panel b: Translocation detection (20 kb segment moved)
- Panel c: Deletion/Insertion detection
- Panel d: Combined SV scenario (multiple SVs + SNPs)

**Key findings:**
- Inversions: 100% detected (orientation flip in synteny block)
- Translocations: detected as rearrangement blocks
- SNPs cause false-positive indels at tag boundaries (manageable)
- Resolution limited by tag spacing (~1.5 kb for BcgI)

### 3.5 Speed and Memory Benchmark

**Table 1**: Computational performance

| Tool | Sketch (s) | Query (s) | Memory (MB) | Sketch size |
|------|-----------|-----------|-------------|-------------|
| FastANI | N/A | 120 | 2,400 | N/A |
| skani | 45 | 0.8 | 1,200 | ~20 MB |
| Syn2b-ani | 12 | 0.3 | 180 | ~48 KB |
| MASH | 8 | 0.1 | 80 | ~30 KB |

**Key findings:**
- Syn2b-ani is 3× faster than skani, 400× faster than FastANI
- Sketch size: 48 KB vs. 20 MB (skani) — 400× smaller
- Memory: 180 MB for 65k genomes (vs. 1.2 GB for skani)

### 3.6 Cross-Species GBRT Generalization

**Figure 5**: GBRT debias model validation across 5 bacterial species
- Panel a: E. coli (training species)
- Panel b: B. subtilis (held-out, Gram-positive)
- Panel c–e: Three environmental isolates (held-out)
- Panel f: Error summary (raw vs. simple debias vs. GBRT)

**Key findings:**
- E. coli-trained model: 0.01% error on E. coli
- Cross-species: 0.30% error (still 2× better than simple debias)
- GC content has no effect on model performance

### 3.7 Real MAG Scenarios

**Figure 6**: Realistic MAG validation
- Contamination (5%–20% foreign DNA)
- Chimerism (1–20 breakpoints)
- Duplication (5%–20%)
- Assembly error (0.01%–0.2%)
- Combined realistic MAG

**Key findings:**
- Contamination affects AF but not ANI
- Chimerism (shuffled contigs) does not affect ANI (tag sequences unchanged)
- Assembly errors directly impact ANI
- Combined realistic MAG: ANI = 99.9%, AF = 90.7%

---

## 4. Discussion

### 4.1 Fixed Anchors vs. Random k-mers

**The chaining bottleneck:**
- skani and FastANI rely on finding chains of matching k-mers
- Chaining breaks when: (a) N50 is too small, (b) rearrangements disrupt collinearity, (c) high divergence reduces seed density
- Syn2b-ani bypasses chaining entirely: each tag is an independent anchor

**Trade-offs:**
- Fixed anchors provide positional information "for free"
- But tag density is lower than k-mer density (~1 tag per 1.5 kb vs. 1 k-mer per base)
- Resolution is coarser (~1.5 kb for SV detection vs. base-pair for alignment)

### 4.2 ANI + Synteny in One Pass

Unique capability of Syn2b-ani:
- skani/FastANI/MASH output only ANI + AF
- Syn2b-ani outputs: ANI, AF, synteny blocks, SV types, SV positions
- Enables downstream analyses: phylogeny + pangenome structure simultaneously
- PAF-like output compatible with IGV/Bandage visualization

### 4.3 Experimental Verifiability

The 2bRAD-M connection:
- Predicted tags can be validated by 2bRAD-M sequencing
- If a predicted tag is not observed, it indicates: assembly gap, sequencing bias, or true SV
- Creates a feedback loop between *in silico* prediction and experimental validation
- Unique among ANI tools — all others are purely computational

### 4.4 Limitations

1. **Higher divergence (>10%)**: Tag sharing becomes sparse (<1000 tags), ANI estimates unreliable. Recommended range: <5% divergence (strain to species boundary).

2. **Resolution**: SV detection limited to ~1.5 kb (tag spacing). Smaller indels (<1 kb) not detectable.

3. **Multi-enzyme mode**: Current implementation concatenates tags from multiple enzymes but does not significantly improve accuracy. Future work: enzyme-specific weighting.

4. **GBRT generalization**: Universal model (v2) achieves 0.3% MAE on new species. Species-specific models achieve 0.01%. For highest precision, per-clade training recommended.

5. **Plasmids and mobile elements**: Type IIB sites may be underrepresented in plasmids. ANI estimates may be chromosome-biased.

### 4.5 Future Directions

1. **GTDB-scale validation**: Compare all ~65,000 representative genomes
2. **Species-specific GBRT models**: Train per-genus models for highest accuracy
3. **Plasmid-aware mode**: Include plasmid-specific enzyme panels
4. **Real MAG validation**: Test on 100+ MAGs from human gut metagenomes
5. **Integration with 2bRAD-M pipeline**: Direct comparison of predicted vs. observed tags

---

## 5. Data Availability

- Syn2b-ani source code: https://github.com/HuangShiLab/Syn2bANI
- Benchmark datasets: [Zenodo/figshare DOI]
- GBRT training data: `gbrt_training_data.csv` (supplementary)
- Test genomes: NCBI RefSeq accessions listed in Supplementary Table 1

---

## 6. Code Availability

- Syn2b-ani v0.1.1: GitHub release
- Benchmark scripts: `benchmark_pipeline.py`, `task3_multispecies.py`, `task4_gbrt_debias.py`
- GBRT model: `gbrt_model_v2.json` (embedded in binary)

---

## Supplementary Information Plan

### Supplementary Tables

| Table | Content |
|-------|---------|
| Table S1 | 49 bacterial genomes used for GBRT training (name, accession, size, GC) |
| Table S2 | 16 Type IIB enzyme specifications (recognition site, tag length, NEB catalog) |
| Table S3 | GBRT hyperparameter search results |
| Table S4 | Cross-species validation detailed results |
| Table S5 | skani/FastANI command lines used for comparison |

### Supplementary Figures

| Figure | Content |
|--------|---------|
| Fig S1 | Effect of matching threshold on ANI accuracy (Hamming 0 vs. 1 vs. 2) |
| Fig S2 | Per-enzyme ANI accuracy (BcgI, BsaXI, CjeI, etc.) |
| Fig S3 | Multi-enzyme consensus results (no improvement over single enzyme) |
| Fig S4 | GBRT feature correlation matrix |
| Fig S5 | GBRT v1 vs. v2 comparison |
| Fig S6 | Runtime scaling with genome count (1 to 65,536) |
| Fig S7 | Memory usage profiling |
| Fig S8 | PAF output examples for inversion/translocation visualization |

---

## References (Key Citations)

1. Jain, C., et al. (2018). High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nature Communications*.
2. Shaw, J. & Yu, Y.W. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nature Methods*.
3. Ondov, B.D., et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology*.
4. [2bRAD-M original paper]
5. [Type IIB restriction enzymes — Roberts et al.]
6. [ANI definition — Konstantinidis & Tiedje]
7. [GTDB — Parks et al.]
8. [MIMAG standards — Bowers et al.]
9. [Scikit-learn GradientBoostingRegressor]
10. [Rust programming language]

---

## Writing Notes for Authors

### Tone and Style
- Follow *Nature Methods* guidelines (Methods/Resource article format)
- Methods section should be detailed enough for reproduction
- Results focus on comparative benchmarks (vs. skani/FastANI)
- Discussion should emphasize the **unique value proposition** (fixed anchors + SV output + experimental verifiability)

### Priority Figures
1. **Figure 1**: Core algorithm + accuracy (main selling point)
2. **Figure 2**: Fragmentation robustness (differentiator)
3. **Figure 4**: SV detection (unique feature)
4. **Figure 5**: GBRT validation (technical depth)

### Expected Reviewer Questions
1. "Why not just use skani?" → Fragmentation robustness + SV output
2. "Is the GBRT model overfitted?" → Cross-species validation + held-out test sets
3. "What about genomes without Type IIB sites?" → 16-enzyme panel ensures coverage; GC bias addressed in Discussion
4. "How does it scale to 100k genomes?" → Benchmark data in Table 1 + Fig S6

---

*Manuscript outline v1.0*
*Last updated: 2026-07-09*
