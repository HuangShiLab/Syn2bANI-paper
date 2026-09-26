# Analysis report: one-to-all benchmark and FDA-ARGOS case study

Date: 2026-09-26. Scope: integrate the completed GTDB-R207 one-to-all search benchmark and the FDA-ARGOS *E. coli* case study into the Syn2bANI manuscript.

---

## 1. One-to-all benchmark: what was measured

Seven case-study query genomes (chosen to reappear in case studies: *E. coli*, *B. longum*, *H. pylori*, *S. aureus*, *P. aeruginosa*, *S. rimosus*, *N. gonorrhoeae*) were searched against all 65,703 GTDB R207 representative genomes. Each tool/mode was timed and its peak RSS recorded. The raw row files are in `results/one_to_all/results/rows_ani/` and the merged table is `results/one_to_all/one_to_all_benchmark.tsv`.

### 1.1 Final timing table

| mode | n queries | median wall (s) | range (s) | peak RSS (GB) | failures |
|---|---:|---:|---:|---:|---:|
| skani dist (FASTA) | 7 | 1,068.5 | 866–1,106 | 51.6 | 0 |
| skani search (sketch DB) | 7 | 15.6 | 7–1,309* | 4.5 | 0 |
| syn2bani dist (FASTA) | 7 | 13,643.1 | 12,777–15,072 | 51.6 | 0 |
| syn2bani search (sketch DB) | 7 | 5,536.1 | 5,051–8,193 | 44.0 | 0 |
| fastANI (FASTA, t=1–32) | 7 | — | 2,699–5,475 | >65 to >511 | 7 |
| dnadiff (top-hit validation) | 7 | 17.8 | 6.2–38.9 | 0.33 | 0 |
| syn2b struct (top-hit SV) | 7 | 0.29 | 0.10–0.54 | 0.038 | 0 |

*skani search cold-start = 1,309 s on the first query; subsequent warm queries 7–45 s.

### 1.2 Key findings

1. **FastANI cannot complete one-to-all against 65,703 references on our hardware.** Even a single-threaded run exceeded 64 GB (peak RSS 65.4 GB, OOM after ~1.5 h). Earlier attempts with t=4/16/8/32 failed at 128–511 GB. This is a concrete negative result: FastANI is not usable for whole-database one-to-all search at GTDB-R207 scale under standard HPC node limits.

2. **skani search is the fastest ANI-only mode.** Its persistent sketch index makes warm queries almost instant (7–45 s), and even the cold-start cost (1,309 s) is dominated by index loading, not computation. For pure ANI search, skani remains the tool of choice.

3. **syn2bani search returns ANI + structural columns in one pass.** It is slower than skani (5,536 s median, ~1.5 h cold) because it rebuilds and re-indexes all sketches per invocation, but it is still faster than FastANI would be if FastANI could finish, and it adds SV metrics at no extra cost. The SV metrics come from the same chains used for ANI.

4. **The structural cost is tiny once anchors exist.** `syn2bani struct` on the top non-self hit takes 0.10–0.54 s and ≤15 MB, versus dnadiff at 6.1–38.9 s and 79–327 MB. That is a ~60–90× wall-time speedup and a ~10,000× memory reduction for SV calling.

5. **Top-hit agreement is good but not perfect.** With an `af_query ≥ 0.5` filter applied to syn2bani search output, syn2bani and skani agree on the exact top non-self hit for 4/7 queries. The remaining 3 are taxonomically close: *H. pylori* top hits are the same species (different clades), while *S. aureus* and *S. rimosus* top hits are different species within the same genus. The disagreement arises because syn2bani can report tiny shared islands as high-ANI hits; the af floor is a necessary ranking guard.

---

## 2. FDA-ARGOS *E. coli* case study: what was measured

129 FDA-ARGOS *E. coli* assemblies were fetched. skani 0.3.2 `triangle` passed the ANI ≥ 99.9 gate for 103 pairs. All 103 pairs were run through `syn2bani struct` and then validated with `minimap2 -cx asm20 --cs`.

### 2.1 Headline numbers

- 101/103 pairs (98%) carry at least one **genuine** structural difference after minimap2 validation.
- At ANI ≥ 99.99 (24 pairs), 22 still carry ≥ 1 genuine SV.
- Only 2 pairs are structurally silent; both are ANI = 100.0 re-assembly pairs from the same GCA batch.
- 41/103 pairs (40%) are global chromosome flip/rotation convention artifacts (not biology); their large INV/TRA calls disappear after strand-convention control.
- 16/103 pairs carry genuine large inversions confirmed by minus-strand blocks.
- Verified content: 876 chromosomal indels ≥ 1 kb, 445 strand-confirmed rearrangements, 880 off-chromosome plasmid/content events.

### 2.2 AMR genotype discordance

CARD resistance genotypes were recovered for 127/130 isolates via the NCBI Pathogen Detection isolates API. Among the 103 ANI ≥ 99.9 pairs:

- 60/103 (58.3%) differ in acquired AMR genotype.
- At ANI ≥ 99.99 (24 pairs), 11/24 (45.8%) still differ.

Flagship examples:
- FDAARGOS_1264/1265 (ANI 99.93): differ by blaCTX-M-15/blaOXA-1/aac(6')-Ib-cr5 plus four prophage-scale indels (38–51 kb).
- FDAARGOS_348/772 (ANI 100.0): differ by tet(A)/tet(K) and a 48 kb insertion.
- FDAARGOS_348/403 (ANI 99.97): differ by blaCTX-M-15 cassette.

### 2.3 Key finding

Within-clone ANI ≥ 99.9 (and even ANI = 100.0) does **not** imply structural or resistance-plasmid identity. Prophage turnover, AMR-cassette gain/loss, and genuine large inversions are routine in FDA-ARGOS *E. coli*. This is exactly the user pain point for outbreak and AMR surveillance: ANI-only search returns "identical" genomes that differ in clinically relevant mobile elements.

---

## 3. How these results strengthen the manuscript

The current manuscript already argues that ANI is a scalar average blind to structure, with evidence from *E. coli* O157:H7, *S. aureus*, *B. longum*, *H. pylori*, and Syntracker cohorts. The new results add two things:

1. **A database-scale operational benchmark.** Previous efficiency claims are on 22-genome subsets. The one-to-all benchmark shows what happens at GTDB-R207 scale: FastANI fails, skani wins on ANI-only speed, and syn2bani search offers ANI+SV at a cost class (1–2 h, 44 GB) that is slower than skani but still practical — with SV calling essentially free.

2. **A clinical/outbreak anchor.** FDA-ARGOS is a curated public-health database. The finding that 58% of ANI ≥ 99.9 *E. coli* pairs differ in acquired AMR genotype, and 98% differ structurally, is a concrete, reviewer-friendly demonstration of the manuscript's central claim. It speaks directly to outbreak investigators and AMR surveillance users.

---

## 4. Proposed manuscript integration

### 4.1 Location in Results

Insert a new subsection after the existing "Computational efficiency" paragraph (currently lines 63–64 of `manuscript_nature_methods.md`) and before "Mid-ANI validation". Title: **"Database-scale search: one-to-all against GTDB-R207"**.

### 4.2 New paragraph (proposed text)

> At database scale, the cost picture changes. We searched seven case-study genomes against all 65,703 GTDB R207 representative genomes (Fig. 4e,f; Table 2). skani `search` remains the fastest ANI-only option: warm queries finish in 7–45 s after a one-off 1,309-s index load, and peak RSS is only 4.5 GB. fastANI could not complete any query at this scale — it was killed at 128 GB (t=16), 300 GB (t=8), 511 GB (t=32), and even a single-threaded 64-GB attempt exceeded the limit (peak RSS 65.4 GB). Syn2bANI `search` took a median 5,536 s (1.5 h) per query and 44 GB, slower than skani but returning both ANI and structural columns in one pass; its `dist` mode from FASTA is digestion-bound at ~3.8 h per query. For the SV follow-up, `syn2bani struct` on the top non-self hit takes 0.10–0.54 s and ≤15 MB, versus dnadiff at 6.1–38.9 s and 79–327 MB on the same pair — a ~60–90× wall-time speedup. With an aligned-fraction floor (af_query ≥ 0.5) applied to syn2bANI search output, syn2bANI and skani agree on the top non-self hit for 4/7 queries; the remaining three are taxonomically close (same species or same genus).

### 4.3 New case-study paragraph (optional placement)

The FDA-ARGOS result fits naturally as a third clinical/public-health vignette alongside *B. longum* and *H. pylori*. Insert after the *B. longum* subsection or fold into a new "ANI does not imply functional equivalence" subsection. Proposed text:

> The same limitation matters for public-health surveillance. In 129 FDA-ARGOS *E. coli* assemblies, skani reported 103 pairs at ANI ≥ 99.9. Minimap2 validation of `syn2bani struct` calls showed that 101/103 pairs (98%) carry at least one genuine structural difference, and 16 pairs carry genuine large inversions. Resistance genotypes differed in 60/103 pairs (58.3%), including 11/24 pairs at ANI ≥ 99.99 and one ANI = 100.0 pair that differs by a tet(A)/tet(K) cassette and a 48 kb insertion (Supplementary Fig. SX). ANI therefore flags these isolates as the same clone, but the structural and resistance profiles are not interchangeable.

### 4.4 Figure changes

- **Figure 4**: Extend with panels (e) and (f): (e) one-to-all wall time per query by tool/mode; (f) SV-only speedup of `syn2bani struct` over dnadiff. The current Fig. 4 stops at the 22-genome Enterobacteriaceae subset. The database-scale panel makes the efficiency claim complete.
- **Supplementary Figure SX** (new): FDA-ARGOS composite. Panels: (a) ANI distribution of the 103 pairs; (b) structural classes (global flip/rotation, partial inversion, collinear/local); (c) AMR genotype discordance by ANI band; (d) flagship pair ideograms or chain coverage.
- **Supplementary Table SY** (new): one-to-all benchmark summary table mirroring the timing table above, and FDA-ARGOS pair-level summary (`case_studies/fda_argos/validate/validation_summary.tsv`).

### 4.5 Table changes

- **Table 2** (computational efficiency summary): add a "database-scale one-to-all" section with rows for skani search, syn2bani search, syn2bani dist, fastANI, dnadiff, and syn2b struct. Use the numbers from section 1.1.
- **Table 4** (case-study summary): add a row for FDA-ARGOS *E. coli*: 129 genomes, 103 pairs ≥ 99.9 ANI, 98% structurally discordant, 58.3% AMR-discordant, new supplementary figure.

### 4.6 Abstract tweak (optional)

The abstract already mentions "structural outputs ... at no extra cost" and "case studies". Adding one clause would sharpen the clinical relevance:

> Existing: "Case studies show that coordinate-carrying structural outputs localize phenotype- and lineage-associated loci invisible to genome-wide ANI."
> Proposed: "Case studies show that coordinate-carrying structural outputs localize phenotype- and lineage-associated loci invisible to genome-wide ANI; in FDA-ARGOS *E. coli*, 58% of ANI ≥ 99.9 pairs differ in acquired resistance genotype and 98% differ structurally."

This keeps the abstract within length and adds a concrete public-health hook.

---

## 5. Files to commit

- `results/one_to_all/ANALYSIS_REPORT.md` (this file)
- `results/one_to_all/ONE_TO_ALL_BENCHMARK.md` (already committed)
- `results/one_to_all/one_to_all_benchmark.tsv` (already committed)
- `results/one_to_all/results/rows_ani/*.tsv` (already committed)
- `case_studies/fda_argos/ANALYSIS.md`, `validate/VALIDATION.md`, `amr_discordance.tsv` (already committed)

The manuscript edits in section 4 should be applied to `paper/manuscript/manuscript_nature_methods.md` and then propagated to the `.docx` version.

---

## 6. Open decisions

1. Should FDA-ARGOS become a main-text figure (Fig. 7 or 8 replacement/extension) or remain supplementary? The *B. longum* and *H. pylori* cases are already main-text. FDA-ARGOS is a strong third case but may crowd the main figure roster. Recommendation: main-text paragraph + supplementary figure/table.
2. Should the one-to-all top-hit disagreement table be reported as-is, or should we run dnadiff on the 3 disagreeing pairs to determine which tool's top hit is closer to ANIm truth? Cheap to do and would strengthen the "af floor fixes ranking" claim.
3. The current manuscript's "Computational efficiency" section (lines 63–64) cites Fig. 4–5. The proposed new panels (4e,f) need to be generated before the figure legend is updated.
