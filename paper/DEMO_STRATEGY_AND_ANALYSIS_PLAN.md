# Demo strategy, user-need evaluation, and analysis plan

Date: 2026-09-25. Purpose: sharpen the Syn2bANI story for skani/FastANI users
and map each claim to a concrete demo, dataset, and analysis.

## 0. Positioning (one sentence)

Syn2bANI does not replace skani/FastANI on their axis — it adds the axis they
do not report (structure) at the same cost class: **ANI tells you two genomes
are the same; Syn2bANI also tells you when they are structurally different.**

---

## 1. User-need evaluation (microbiome + clinical + database curation)

| persona | current tool & workflow | pain point Syn2bANI addresses | strongest demo |
|---|---|---|---|
| Metagenome/MAG analyst | `skani dist/search` MAG vs GTDB; then dnadiff on top hits if structure matters | MAGs are fragmented: dnadiff breakpoint counts are dominated by contig boundaries; skani reports nothing structural | Fragmented MAG search on GTDB50k: Syn2b metrics unchanged by contig count (Fig 3d of Syn2b paper); rank changes when structure matters |
| Microbial taxonomist / database curator | FastANI 95% species delineation, GTDB clustering | Within-species structural divergence is invisible; mis-assemblies/orientation issues pass QC | GTDB census discordance catalog (running): list of "ANI-concordant, structurally divergent" pairs per species cluster |
| Clinical/public-health genomics (outbreaks) | core-genome SNP trees, cgMLST; ANI for species ID | Outbreak isolates are ANI-identical by construction; differences live in plasmids/prophages/large indels | PulseNet/FDA-ARGOS isolates: ANI ~100% but resistance plasmid gain/loss and structural changes flagged |
| AMR surveillance | CARD/ResFinder on assemblies + ANI clustering | Resistance genotype differences between "identical" isolates not surfaced during search | H. pylori 528 cohort: at ANI ≥ 99.9%, clarithromycin-resistance genotype differs in 5.3% of pairs (already computed) |
| Strain-collection QC | FastANI identity check vs reference | Orientation/origin artifacts and large inversions invisible | K-12 vs W3110 vignette (ANI 99.99, AF 100%; 782-kb inversion + 9 indels detected) |

Key insight from the H. pylori numbers: discordance persists at ANI ≥ 99.9%
(cagPAI state 10.5%, vacA 2.6%, clarithromycin resistance 5.3%), so the
problem is not limited to distant relatives — it is exactly where ANI tools
are *trusted most*.

---

## 2. Demonstration matrix (what we show, where it lives)

| demo | data | status | feeds |
|---|---|---|---|
| One-to-all search equivalence (time + extra columns) | GTDB50k reference set | pending (benchmark script exists) | main figure / Supp Table |
| GTDB census discordance catalog | full GTDB R207 (census running on HPC) | running | resource figure + online catalog |
| K-12 vs W3110 vignette | NC_000913.3 vs NC_007779.1 | **done** (case_studies/ecoli_k12_w3110, fig_k12_w3110_vignette) | vignette figure |
| H. pylori ANI-discordant biology | song 2026 cohort (528 genomes, metadata recovered) | **done** (pairwise_ani/) | table/figure S |
| B. longum abfA (ANI-identical isolates, functional island intact vs lost) | 185 JNU isolates + abfA cluster | existing results, reframe | case figure |
| MAG fragmentation robustness | CAMI2 marine/strain + GTDB50k | data exist; decide scope | Supp figure |
| PulseNet/FDA-ARGOS outbreak plasmid demo | FDA-ARGOS E. coli (129 genomes) | **done** (case_studies/fda_argos; see below) | case figure |

## 3. Candidate datasets for the missing demos

Clinical/outbreak (highest persuasion value):
1. **FDA-ARGOS** curated isolate database (all have phenotypic/clinical
   metadata; RefSeq-linked). Pick *Salmonella enterica*, *E. coli* ST131,
   *Klebsiella pneumoniae* ST258 clusters: same-sequence-type isolates with
   documented plasmid/resistance differences.
2. **Salmonella Typhi H58** or **E. coli ST131-H30Rx** isolate collections
   (large public sets with resistance phenotypes) — ANI ~99.99 within clone,
   AMR gains/losses are structural.
3. **B. longum abfA** (in-house, Cell Host Microbe 2023) — already have it;
   the only dataset where the structural difference is tied to a *measured
   human-treatment outcome*.

Microbiome:
4. **UHGG / proGenomes3** species-level collections for one-to-all search
   demos against MAGs (use GTDB50k where possible to avoid another dataset).
5. **CAMI2 marine/strain** — already on HPC; use only the reference genomes
   for MAG fragmentation robustness, not the read tarballs.

Rule: prefer datasets with (a) public genomes, (b) phenotype/metadata columns
that matter clinically (resistance, virulence, disease stage), (c) dense
within-clone sampling (so ANI ≥ 99.9 pairs exist).

## 4. Analysis plan (ordered)

1. **Now (local, done):** K-12/W3110 vignette figure + H. pylori ANI-discordant
   biology tables (`results/pairwise_ani/`).
2. **HPC, this week:** one-to-all search benchmark vs skani at GTDB50k scale
   (time, memory, output columns; script exists).
3. **HPC, census-dependent:** package census output as the discordance catalog
   (per-cluster discordant pair lists + summary stats); one figure + one
   downloadable table.
4. **Next (new data, small):** pick ONE clinical collection from section 3
   (FDA-ARGOS Salmonella or E. coli ST131), ~200–500 genomes; run
   syn2bani dist+struct; cross with AMR phenotype; target table: "ANI-identical
   pairs with differing resistance phenotype: n/X".
   **Status 2026-09-25: DONE (FDA-ARGOS E. coli, 129 genomes).**
   skani 0.3.2 triangle: 103 pairs >= 99.9 ANI (gate >= 50 passed).
   syn2bani struct on all 103 pairs, then minimap2 validation of every call:
   - 41/103 pairs are whole-chromosome flip/rotation convention artifacts
     (their Mb-scale INV calls sit at rotation junctions — the same artifact
     class as the H. pylori TRA:1459-1666206 lesson)
   - corrected headline: 101/103 pairs (98%) carry >= 1 genuine structural
     difference; 876 indels >= 1 kb verified in minimap2 cs strings
   - 16 pairs carry genuine large inversions
   - AMR genotypes recovered for 127/130 isolates via the NCBI Pathogen
     Detection isolates API (BioSample has none)
   - flagship pairs: FDAARGOS_1264/1265 (ANI 99.93, blaCTX-M-15 cassette +
     22 verified indels), 348/403 (O104:H4, blaCTX-M-15), 348/772 (ANI 100.0,
     tet cassettes + 48 kb insertion)
   Remaining: AMR x structure cross-table ("ANI-identical pairs with differing
   AMR genotype: n/X") for the paper.
5. **Write-up:** reframe B. longum and H. pylori cases as "what the search
   would have missed"; K-12 vignette as the every-microbiologist anchor.

## 5. Decision gates

- If census discordance prevalence confirms ≥10% at ≥97% ANI (preliminary:
  34% carry ≥2 inversions in the 50k sample), the catalog becomes a headline
  figure; otherwise it stays supplementary.
- Clinical collection (step 4) proceeds only if its within-clone ANI range
  contains ≥50 pairs ≥99.9 ANI — checked with skani before any Syn2b run.
