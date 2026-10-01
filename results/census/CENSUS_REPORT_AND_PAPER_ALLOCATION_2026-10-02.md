# Three-database within-species census: status, results, and paper allocation

**Snapshot:** 2026-10-02 00:10 HKT.  
**Scope:** GTDB R207, HROM, and OAPGC.  
**Purpose:** consolidate the current census evidence, define quality-aware rules for cross-species comparisons, and assign results to the Syn2b and Syn2bANI papers without duplicating claims.

This report applies the Supervisor-Skills framework as follows:

- `benchmark-paper-template`: the census is treated as a resource/benchmark contribution, requiring a construction pipeline, evaluation framework, empirical findings, and quality controls.
- `tech-paper-template`: Syn2b is positioned as the scalable structural-inference method, whereas Syn2bANI is positioned as the joint ANI–structural interpretation layer.
- `vibe-research-workflow`: biological claims are restricted to what the current outputs support; incomplete runs are status, not conclusions.

---

## 1. Executive status

| Database | Biological scope | Multi-genome clusters | Genomes analysed | Unique unordered pairs planned | Structural census | ANI layer |
|---|---|---:|---:|---:|---|---|
| **GTDB R207** | global bacterial/archaeal diversity | 22,535 | 274,374 | **711,020,841** | Raw structural tasks complete: 25,582/25,582. Low-memory merge/QC running. | Streaming skani pass running; 21,965/22,535 cluster outputs done; two large clusters actively running. |
| **HROM** | human gut | 2,241 | 142,277 | **62,872,901** | Complete; 0 missing tasks or clusters. | Complete; 62,872,901/62,872,901 pairs matched to skani ANI. |
| **OAPGC** | oral and airway | 1,488 | 98,229 | **30,560,247** | Running locally: 521/2,120 tasks done; 24,380,902 pair rows already emitted. | Not started; run after structural merge. |
| **Combined planned** | — | — | — | **804,453,989** | — | — |

Current interpretation boundary:

- GTDB has a complete raw structural layer, but not yet a final merged/QC biological summary.
- HROM is the only fully completed dual-axis census and therefore currently carries the strongest biological conclusions.
- OAPGC has emitted roughly 80% of its planned pair rows by row count, but the run is not finished and has no merged/QC biological result.
- Cross-species conclusions must be quality-stratified; unadjusted species rankings are not acceptable.

---

## 2. GTDB R207 structural census update

### 2.1 Completed structural layer

The GTDB task plan contains:

| Item | Value |
|---|---:|
| Multi-genome species clusters | 22,535 |
| Genomes in those clusters | 274,374 |
| Unique unordered within-species pairs | 711,020,841 |
| Structural block-pair tasks | 25,582 |
| Completed structural outputs | 25,582 |
| Completed checkpoints | 25,582 |

Thus, the GTDB raw structural census is complete. The remaining work is computational consolidation, not genome collection or structural calculation.

The low-memory merge is now running after replacing the original implementation. The old merger retained global pair-state in memory, which was unsuitable for 711 million pairs. The replacement writes the final pair table streamingly, performs local duplicate suppression, keeps exact per-cluster row counts and exact `breakpoints >= 2` counts, and uses reservoir sampling only for per-cluster medians.

### 2.2 GTDB ANI pass

The first full ANI implementation accumulated a dense triangle matrix in Python and was OOM-killed on very large clusters. The revised implementation now:

1. runs skani triangle with one thread per worker;
2. parses skani output line-by-line;
3. writes only the lower-triangle non-self comparisons;
4. avoids retaining a dense Python matrix;
5. processes very large clusters at reduced concurrency;
6. uses a 96 GB SLURM allocation.

Current status:

```text
ANI cluster outputs: 21,965 / 22,535
Remaining:            570 clusters
Active tmp outputs:      2
Active examples: Escherichia coli and Staphylococcus aureus
```

The lack of new completed ANI outputs is expected because the revised scheduler intentionally starts with the largest clusters. This should not be interpreted as a stall while the job is running and temporary outputs exist.

### 2.3 Natural genome-quality variation inside GTDB species

A complete GTDB quality table was built for 317,542 census genomes, joining bare accessions to GTDB’s prefixed metadata accessions and retaining CheckM completeness, contamination, strain heterogeneity, contig count, N50, longest contig, and GC.

For the 22,535 multi-genome GTDB species clusters:

| Quality group | Genomes | Fraction |
|---|---:|---:|
| High: completeness ≥90%, contamination ≤5% | 243,727 | 88.85% |
| Medium: completeness 70–90%, contamination ≤5% | 20,553 | 7.49% |
| High contamination: completeness ≥70%, contamination >5% | 3,446 | 1.26% |
| Low completeness: completeness <70% | 6,648 | 2.42% |

Quality heterogeneity is not a marginal issue:

| Observation | Value |
|---|---:|
| Species containing more than one quality group | 8,887 / 22,535 |
| Fraction of species with mixed quality groups | 39.44% |
| Species containing at least one non-high genome | 11,668 / 22,535 |
| Species with high-quality fraction ≤50% | 8,360 / 22,535 |
| Species in which every genome is high quality | 10,867 / 22,535 |

**GTDB quality conclusion:** species labels do not define homogeneous measurement conditions. Almost 40% of GTDB species mix different genome-quality classes, so a species-level structural rate reflects both biology and assembly quality. Cross-species comparisons must therefore separate biological divergence from quality-dependent detection power.

---

## 3. HROM census: completed dual-axis result

HROM is complete and is currently the strongest biological evidence base.

### 3.1 Census size and QC

| Item | Value |
|---|---:|
| Total HROM genomes | 145,149 |
| Multi-genome clusters | 2,241 |
| Genomes included | 142,277 |
| Unique unordered pairs | 62,872,901 |
| Structural tasks | 3,525 |
| Missing tasks or clusters | 0 |
| Pairs with matched skani ANI | 62,872,901 / 62,872,901 |

### 3.2 Structural prevalence

```text
Pairs with breakpoints >= 2: 17,535,013
Fraction:                    27.8896%
```

**HROM Finding 1:** structural block-level variation is common within gut bacterial species; 27.9% of all within-species pairs carry at least two visible breakpoints.

### 3.3 ANI-stratified structural signal

| skani ANI bin | Pairs | % with breakpoints ≥2 |
|---|---:|---:|
| <95% | 2,376,169 | 45.128% |
| 95–97% | 28,516,234 | 31.085% |
| 97–98% | 22,576,854 | 20.837% |
| 98–99% | 8,946,967 | 31.610% |
| 99–99.5% | 363,255 | 14.104% |
| ≥99.5% | 93,418 | 15.751% |

**HROM Finding 2:** ANI rank is not a sufficient proxy for structural divergence. The non-monotonic rise in the 98–99% bin relative to 97–98% is a key motivation for the joint ANI–structural view.

### 3.4 Genome-quality sensitivity

The HROM quality analysis shows that SV detection is strongly assembly-dependent:

| Minimum completeness | Pairs | % with breakpoints ≥2 |
|---|---:|---:|
| <70% | 16,337,934 | 16.171% |
| 70–90% | 25,466,226 | 24.712% |
| ≥90% | 18,210,578 | 45.692% |

N50 shows the same direction:

| Minimum N50 | Pairs | % with breakpoints ≥2 |
|---|---:|---:|
| <10 kb | 36,999,517 | 19.260% |
| ≥10 kb | 23,015,221 | 44.014% |

**HROM Finding 3:** structural detection power depends strongly on completeness and contiguity. A quality-unadjusted cross-species comparison would confound genome biology with assembly fragmentation.

### 3.5 High-quality sensitivity result

Restricting pairs to completeness ≥90% and contamination ≤5% gives:

| Item | All qualified | High-quality subset |
|---|---:|---:|
| Pairs | 60,014,738 | 18,210,578 |
| Species tested | 521 | 354 |
| Top 5% per axis | 3,000,736 | 910,528 |
| SV/structural-enriched species | 96 | 60 |
| ANI-enriched species | 160/161 | 93 |
| Both-axis enriched species | 36 | 19 |

Mode stability from all pairs to the high-quality subset:

| Mode | All | High quality | Intersection | Jaccard |
|---|---:|---:|---:|---:|
| Structural-enriched | 96 | 60 | 53 | 0.515 |
| ANI-enriched | 160 | 93 | 79 | 0.451 |
| Both-enriched | 36 | 19 | 9 | 0.196 |

**HROM Finding 4:** the dual-axis pattern survives quality filtering, but exact species membership is less stable. Papers should report both all-pairs and high-quality-only results rather than selecting whichever list is more appealing.

### 3.6 Hypermode interpretation

Current HROM examples include:

- SV/structural-enriched: *Rothia mucilaginosa*, *Lautropia mirabilis*, *Actinomyces graevenitzii*.
- ANI-enriched / sequence-divergence-enriched: *Prevotella pallens*, *Porphyromonas gingivalis* in the oral catalogue context, *Streptococcus mutans*, *Tannerella forsythia*.
- Both-axis enriched: species where high ANI-distance and high SV burden co-occur, suggesting coupled divergence rather than a pure SNP-only or SV-only mode.

Claim boundary:

- `ANI-enriched` should not be described as `hypermutator` unless validated by direct SNP calling, mutation-spectrum analysis, or popANI/inStrain-type comparison.
- `SV-enriched` should not be described as `hyper-recombinator` unless validated by long-read assembly, dnadiff/minimap2, or phased structural haplotypes.

---

## 4. OAPGC census status

OAPGC is the oral and airway resource intended to answer whether same-species genomes from different body sites differ more by sequence divergence or structural variation.

| Item | Value |
|---|---:|
| HQMAG rows in annotation | 149,921 |
| Strain-level non-redundant representatives | 99,215 |
| Multi-genome SGB clusters | 1,488 |
| Genomes in structural census | 98,229 |
| Planned unique pairs | 30,560,247 |
| Structural tasks | 2,120 |
| Current completed tasks | 521 |
| Current emitted pair rows | 24,380,902 |

The local run is healthy, but the current partial table should not be used for final OAPGC conclusions. After structural completion, the required next steps are:

1. merge and QC the OAPGC structural table;
2. run skani ANI on all 30.56 million pairs;
3. join oral/airway sample metadata;
4. classify pairs as oral–oral, airway–airway, and oral–airway;
5. compare ANI distance and structural burden with species-level normalization;
6. repeat all key tests on the high-quality subset.

---

## 5. Quality-aware framework for cross-species comparison

Because GTDB and HROM both show strong quality dependence, cross-species analyses should use the following standard.

### 5.1 Pair-level quality covariates

For every pair:

```text
min_completeness = min(completeness_A, completeness_B)
max_contamination = max(contamination_A, contamination_B)
min_N50 = min(N50_A, N50_B)
min_length = min(total_length_A, total_length_B)
```

Primary high-quality pair definition:

```text
min_completeness >= 90
max_contamination <= 5
```

Report contiguity separately:

```text
min_N50 >= 10 kb
```

### 5.2 Species-level model

Do not rank species only by raw `fraction_of_pairs_with_breakpoints >= 2`. Instead:

1. fit a quality-adjusted binomial or mixed model;
2. include species as a random effect when pooling pairs;
3. include completeness, contamination, and N50 as fixed covariates;
4. report marginal means or shrunk species estimates;
5. show unadjusted, covariate-adjusted, and high-quality-only estimates.

### 5.3 Within-species quality contrast

For species with enough high- and lower-quality genomes, compare high-quality versus lower-quality genomes within the same species. This directly measures measurement bias while holding lineage and ecological context approximately constant.

### 5.4 Cross-species enrichment

For hypermode/enrichment analyses:

1. define top 5% structural and top 5% ANI pairs within each species;
2. require a minimum number of pairs and genomes;
3. compute enrichment within species;
4. repeat using only high-quality pairs;
5. require the same direction of effect in both all-pairs and high-quality analyses before naming a species;
6. treat unstable small clusters as exploratory.

### 5.5 Minimum reporting standard

Every cross-species figure or table should state:

- all-pairs estimate;
- high-quality-only estimate;
- quality-group composition of each species;
- number of genomes and pairs;
- whether species with mixed quality groups are included;
- whether the claim survives high-quality-only sensitivity analysis.

---

## 6. Paper allocation strategy

The two papers should answer different questions. The recommended division is:

> **Syn2b paper:** scalable structural census engine and global structural resource.  
> **Syn2bANI paper:** joint ANI–structural interpretation and biological/ecological applications.

### 6.1 Syn2b paper

**Core question:** Can Syn2b produce a scalable, reproducible, high-confidence structural census across hundreds of millions of within-species genome pairs?

**Primary evidence**

1. Syn2b method and structural channel.
2. Benchmark against skani/FastANI, dnadiff/minimap2, inStrain/SynTracker-style comparators where available.
3. Runtime, memory, sensitivity, precision, and scalability.
4. GTDB R207 structural census as the flagship global resource.
5. HROM structural census as an independent biological application.
6. OAPGC as a third-scope demonstration of portability to MAG-heavy oral/airway data.

**Primary results**

- GTDB: 711,020,841 within-species pairs across 22,535 species clusters.
- HROM: 62,872,901 pairs; 27.89% with `breakpoints >= 2`.
- OAPGC: 30,560,247 planned pairs; local feasibility demonstrated on Apple Silicon.
- Quality sensitivity: SV detection depends on completeness and N50.
- GTDB quality heterogeneity: 39.44% of species mix quality groups.

**Main figures / tables**

| Item | Suggested location |
|---|---|
| Syn2b structural inference concept and census pipeline | Figure 1 |
| Benchmark sensitivity, runtime, memory, and concordance | Figure 2 |
| GTDB structural census design and scale | Figure 3 |
| HROM structural prevalence and quality sensitivity | Figure 4 |
| OAPGC third-cohort structural summary | Main text table or supplementary figure |
| Full quality-stratified sensitivity | Supplementary tables |

**What Syn2b should not claim**

- It should not make the central hypermutator/hyper-recombinator biological claim.
- It should not become an ANI-methods paper.
- It should not spend main-text space on OAPGC oral–airway evolutionary interpretation.

### 6.2 Syn2bANI paper

**Core question:** Does joint analysis of ANI and structural variation reveal distinct modes of within-species evolution, and can these modes be mapped across human-associated habitats?

**Primary evidence**

1. HROM completed dual-axis census.
2. HROM SV-enriched, ANI-enriched, and both-axis enriched species.
3. High-quality-only sensitivity and mode-stability analysis.
4. OAPGC oral–airway contrast once complete.
5. GTDB joint ANI–structural examples or global summary once ANI completes.
6. Syn2bANI benchmark against ANI-only and structural-only rankings.

**Primary results**

- HROM 62.87 million pairs with complete structural + ANI matching.
- Non-monotonic ANI–SV relationship.
- 60–96 SV-enriched and 93–160 ANI-enriched species depending on quality set.
- 19–36 both-axis enriched species.
- Quality-aware high-quality subset preserves the dual-axis pattern.
- OAPGC will test whether cross-site strains are preferentially ANI-divergent, SV-divergent, or both.

**Main figures / tables**

| Item | Suggested location |
|---|---|
| Joint ANI–structural decision framework | Figure 1 |
| HROM dual-axis species enrichment | Figure 2 |
| High-quality sensitivity and mode stability | Figure 3 |
| OAPGC oral–airway SNP-vs-SV contrast | Figure 4 |
| GTDB global joint species map or examples | Figure 5 or supplementary figure |
| Full enriched species tables | Supplementary tables |

**What Syn2bANI should not duplicate**

- It should not re-derive the Syn2b algorithm or repeat the full runtime benchmark.
- It should not present GTDB structural prevalence as its primary new result.
- It should not repeat all raw HROM quality tables; it should cite the Syn2b paper and show only the sensitivity relevant to species calls.

### 6.3 Shared resources

Both papers should cite the same Zenodo-hosted raw pair tables:

```text
GTDB: gtdb_r207_within_species_sv.tsv.gz
HROM: hrom_within_species_sv.tsv.gz
HROM joint: hrom_within_species_sv_ani.tsv.gz
OAPGC: oapgc_within_species_sv.tsv.gz, after completion
```

Both papers may use the same summaries, but main-text claims should not overlap:

- Syn2b owns structural scalability and prevalence.
- Syn2bANI owns dual-axis interpretation and habitat-level evolutionary modes.
- Quality sensitivity is primary in Syn2b, with a condensed version repeated in Syn2bANI for robustness.

### 6.4 Submission strategy

Recommended order:

1. **Finish GTDB merge/QC**, then lock Syn2b census numbers.
2. **Finish OAPGC structural run**, merge and QC locally.
3. **Run OAPGC ANI and cross-site contrast.**
4. **Finish GTDB streaming ANI**, then run a quality-adjusted GTDB joint summary.
5. Submit **Syn2b** once GTDB structural merge and benchmark tables are final.
6. Submit **Syn2bANI** after OAPGC and GTDB joint analyses are complete.

If schedule pressure requires earlier submission, Syn2b can be submitted first without waiting for OAPGC because GTDB and HROM are sufficient for the structural-resource claim. Syn2bANI should wait for OAPGC unless it is reframed as an HROM-only dual-axis study, which would be weaker for Nature Methods.

---

## 7. Nature Methods positioning

### Syn2b as a Nature Methods paper

Strong framing:

> A scalable structural-variation census engine that measures within-species genome architecture across hundreds of millions of pairs and reveals quality-dependent structural diversity invisible to ANI-only workflows.

Key Nature Methods requirements:

- methodological advance;
- rigorous benchmark;
- broad applicability;
- large-scale demonstration;
- reproducible resource;
- clear quality controls.

The GTDB census supplies scale. HROM supplies biological realism. OAPGC supplies portability to oral/airway MAGs. This is a coherent method/resource paper.

### Syn2bANI as a Nature Methods paper

Strong framing:

> A joint ANI–structural framework that resolves distinct within-species divergence modes and identifies sequence-divergence-enriched, SV-enriched, and both-axis-enriched lineages across human-associated microbiomes.

This paper is stronger as a **resource/analysis workflow** than as a purely synthetic-method paper. It should emphasize that ANI-only workflows collapse two biologically distinct axes. The strongest version needs HROM plus OAPGC, with GTDB as a broad validation layer.

---

## 8. Immediate next actions

1. Wait for GTDB structural merge to finish; validate exact row count against 711,020,841 and inspect `census_qc.md`.
2. Monitor the streaming ANI large-cluster pass; do not cancel while tmp files and the SLURM job are active.
3. Finish local OAPGC structural tasks, then run low-memory merge/QC.
4. Run OAPGC skani ANI and attach site metadata.
5. Build a GTDB quality-stratified species table after merge.
6. Build the OAPGC oral–airway contrast after structural + ANI are complete.
7. Freeze Syn2b census numbers only after GTDB and OAPGC merge QC pass.
