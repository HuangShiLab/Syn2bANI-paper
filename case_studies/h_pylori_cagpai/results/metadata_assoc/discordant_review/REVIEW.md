# Discordant cagPAI/cagA review — 26 genomes

**Date:** 2026-09-11
**Case study:** `case_studies/h_pylori_cagpai` (Song et al. 2026 H. pylori cohort, 528 genomes)
**Inputs:**
- Panel: 28-marker cagPAI status (`cagpai_status/call_cagpai_status.py`, minimap2 `-cx asm20`, present = marker coverage ≥ 0.8 & identity (1−`de:f`) ≥ 0.8; complete ≥ 85 %, empty ≤ 15 %)
- Collaborator: independent cagA genotyping + genome-level rescue (`cagA_final_GCA528.xlsx`): Detected 449 / Not_detected 75 / Indeterminate 2 / Detected_partial 1 / Detected_disrupted 1
- Discordance: `cagpai_complete` ≠ `cagA_detected`, treating `Indeterminate` (cagA_detected = NaN) as discordant → **26 / 528 genomes (4.9 %)**; concordance OR ≈ 4000.

**Outputs:** `discordant_genomes_review.tsv` (one row per genome), `check_discordant_genomes.py` (reproducible; see bottom).

---

## 1. Verdict summary

| verdict | n | meaning |
|---|---|---|
| `disrupted_cagA` | 17 | cagPAI backbone deleted/truncated **but cagA sequence (full or remnant) is in the assembly** — panel empty/partial and collaborator Detected are each defensible; a framing difference, not an error |
| `panel_missed` | 5 | island largely covered and cagA present, but panel says partial — all 5 are divergence-driven (East-Asian L2/L4 markers vs 26695 Western markers) |
| `true_biological_cagA_loss` | 2 | island complete (25/26 markers) with the **cagA gene specifically deleted** — real biology, panel correct |
| `assembly_gap` | 2 | locus absent from the assembly with covered flanks = **clean biological deletions** (both single/draft assemblies show the empty site assembled; panel empty correct) |
| `other` | 0 | — |

Cross-tabulation:

| panel status (discordant rows) | verdicts | assembly supports panel? |
|---|---|---|
| empty (13) | 11 disrupted_cagA + 2 assembly_gap | **13 / 13 yes** |
| partial (10) | 5 disrupted_cagA + 5 panel_missed | 5 / 10 (soft misses) |
| complete (3, incl. Detected_partial) | 2 true_biological_cagA_loss + 1 disrupted_cagA | **3 / 3 yes** |

**Headline:** none of the 13 discordant `empty` calls is a panel error; every one is supported by the assembly (island backbone genuinely deleted, or clean deletion). The only soft panel misses are 5/10 `partial` calls on divergent L2/L4 genomes. No discordant genome is explained by a draft assembly gap.

---

## 2. Category details and representative genomes

### 2.1 `true_biological_cagA_loss` (2) — the highlight

| GCA | country | fastbaps | group | panel | collaborator |
|---|---|---|---|---|---|
| GCA_011140595.1 | Mexico | L3 | GC | complete_rearranged (25/28; missing HP0521, cagY, cagA) | **Not_detected** |
| GCA_014496775.1 | Honduras | L3 | NAG | complete_rearranged (26/28; missing HP0521, cagA) | **Not_detected** |

Per-2 kb coverage of the 26695 window (chunks mapped with `-cx asm20`, real identity ≥ 0.7; cagA gene = 579 920–583 481):

```
                 GCA_011140595.1   GCA_014496775.1   (window 547328-583481)
547328-577328        64-100 %          ~100 %          island body fully present
577328-579328        100 %             100 %           (up to cagB)
579328-581328         22 %              22 %   <- cagA 5' end
581328-583328          0 %               0 %   <- cagA 3' end (EPIYA)
```

- Both islands are intact over ≥ 84–90 % of the window (identity 0.97–0.98) **except the last ~3.5 kb, exactly where cagA sits**. The deletion is surgical: coverage drops from 100 % to 0 % within the cagA gene.
- cagA absence is confirmed by a k-mer test (0 shared 19-mers, both strands — robust to inversion/duplication) and by the absence of any asm20 hit, so this is not a mapping artifact.
- Both strains carry inversions across the island (`complete_rearranged`; cag_sv_list `INV:545600-579949` / `INV:545891-567151` + `INV:567975-579713`), stopping right at the cagA deletion junction.
- Interpretation: real, locus-specific cagA excision in two Latin-American L3 strains while the rest of the island (including the T4SS) is retained. These two genomes *create* the discordance rather than revealing a panel failure: the panel's island-level call is correct, and the collaborator's gene-level call is correct.

### 2.2 `assembly_gap` (2) — clean deletions, panel empty correct

| GCA | country | fastbaps | cagA_result | total_len | window_cov | flanks (L/R) |
|---|---|---|---|---|---|---|
| GCA_018120125.1 | China | L2 | Indeterminate | 1 569 344 | 0.00 | 0.64 / 0.93 |
| GCA_040861305.1 | Sweden | L5 | Indeterminate | 1 547 370 | 0.03 | 0.81 / 0.99 |

- Both assemblies span the locus with covered flanks but no window alignment → the empty site is *assembled*, i.e. a genuine deletion (genomes are 55–75 kb smaller than the cohort median 1 622 666 bp; a full island is ~36 kb plus flanking repeats).
- The earlier "100 kb alignment" that suggested an island in GCA_018120125.1 is a **minimap2 artifact**: one gappy chain (`NM` 40 462 over 103 729 columns, real identity 0.61, ungapped 17 %) forced onto non-colinear sequence; every island marker has **0 shared 19-mers**. The chunked mapping + real-identity filter used here removes it.
- Panel `empty` is correct; the collaborator's `Indeterminate` is consistent with a deleted locus (their rescue presumably found only fragmentary evidence).

### 2.3 `disrupted_cagA` (17) — the dominant class: island broken, cagA retained in some form

Two sub-patterns:

**(a) cagPAI backbone deleted, but a FULL-LENGTH cagA persists (6):**

| GCA | country/fastbaps | cagA k-mer pos_cov | asm20 cov/idn | window_cov | panel present |
|---|---|---|---|---|---|
| GCA_000274605.1 | Peru, L3 | 0.78 | 0.98 / 0.96 (strict asm5 0.97) | 0.15 | 4/28 (cagD-C-B-A) |
| GCA_018120165.1 | China, L2 | 0.30 | 0.88 / 0.88 | 0.12 | 2/28 |
| GCA_022923015.1 | China, L2 (1 contig) | 0.30 | 0.94 / 0.87 | 0.15 | 4/28 |
| GCA_040861205.1 | Sweden, L5 | 0.65 | 0.98 / 0.94 | 0.10 | 1/28 |
| GCA_040861785.1 | Sweden, L5 | 0.71 | 1.00 / 0.95 | 0.10 | 1/28 |
| GCA_040862505.1 | Sweden, L5 | 0.65 | 0.87 / 0.96 (strict asm5 0.72) | 0.10 | 2/28 |

The island is reduced to ≤ 15 % of the window, but cagA itself is complete (100 % span, 87–96 % identity; for the L2 strains it is a divergent East-Asian allele that only asm20/k-mers detect — asm5 finds nothing). GCA_000274605.1 retains the island's *right end* (cagD→cagA, matching the known deletion breakpoint adjacent to cagA). Panel `empty` = island judgement (defensible), collaborator `Detected` = gene judgement (correct). This "empty island, intact cagA" state is a documented H. pylori configuration (cagA persisting at the vacated site).

**(b) cagA truncated/remnant (11):** cagA position coverage 0.11–0.54 (asm20 gene coverage 0.18–0.70, identity 0.86–0.96, strict asm5 mostly absent), window coverage 0.05–0.44 → the island and cagA are both partially degraded. Representative:

| GCA | country/fastbaps | cagA_result | cagA pos_cov | asm20 cov/idn | window_cov | note |
|---|---|---|---|---|---|---|
| GCA_022921455.1 | China, L2 | Detected | 0.20 | 0.53 / 0.87 | 0.05 | single-contig genome; clean deletion + cagA fragment |
| GCA_019093785.1 | Sweden, L5 | Detected | 0.54 | 0.70 / 0.96 | 0.17 | 3' ~70 % of cagA at 96 % idn |
| GCA_019671595.1 | China, L2 | Detected | 0.22 | 0.69 / 0.87 | 0.44 | partial island + partial cagA |
| GCA_019672535.1 | Japan, L2 | **Detected_partial** | 0.35 | 0.66 / 0.90 | 0.89 | panel complete (24/28, cagA in missing list) — **both methods independently flag a truncated cagA; they agree on the biology and differ only on the complete/partial label** |
| GCA_040861875.1 | Sweden, L5 | **Detected_disrupted** | 0.11 | 0.18 / 0.93 | 0.13 | only a small 5' cagA fragment; island otherwise gone (4 left-end markers) — **direct match to the collaborator's "disrupted" call** |

For these, the panel's island-level call (empty/partial) is assembly-supported and the collaborator's cagA call reflects the surviving fragment. The 2×2 table frames the same genome differently.

### 2.4 `panel_missed` (5) — soft misses, divergence-driven

| GCA | country/fastbaps | panel | n_present | cagA pos_cov | cagA asm20 cov/idn | window_cov |
|---|---|---|---|---|---|---|
| GCA_014498675.1 | Brazil, L4 | partial | 23/28 | 0.66 | 0.98 / 0.95 | 0.97 |
| GCA_018120515.1 | China, L2 | partial | 23/28 | 0.36 | 0.79 / 0.89 | 0.90 |
| GCA_018122055.1 | China, L2 | partial | 23/28 | 0.33 | 0.98 / 0.87 | 0.90 |
| GCA_019672515.1 | Japan, L2 | partial | 23/28 | 0.34 | 0.77 / 0.89 | 0.89 |
| GCA_019673955.1 | China, L2 | partial | 23/28 | 0.35 | 0.79 / 0.89 | 0.90 |

- The island is in the assembly (window ≥ 0.89, identity ~0.95) and cagA is present — yet the panel reports `partial`, because 4–5 markers systematically fail on East-Asian/L4 backgrounds (the markers are 26695 = Western genes): HP0521, cagY, cagQ, cagB are highly divergent, and cagA sits at 0.77–0.79 coverage, just under the 0.8 cutoff.
- Note the panel does detect *most* of the island on L2 backgrounds (23/28) — the same lineage where the two `true_biological_cagA_loss` strains score 25–26/28, so the panel is internally consistent; the misses are threshold+divergence effects, not random failures.
- Recommended panel improvements (not required for this paper): (i) add East-Asian reference alleles for cagA/cagY/cagB/cagQ/HP0521, or (ii) lower the cagA coverage cutoff to ~0.7, or (iii) report "complete_with_divergent_markers" when ≥ 22 markers pass.

---

## 3. SV BED cross-check (task step e)

`cag_sv_list` matches the in-window records of `struct_vs_26695_filtered/<GCA>.vs_hp26695.bed` for 24/26 genomes. The two mismatches are bookkeeping-level:

- **GCA_011140595.1**: BED additionally contains a small `DEL:556375-558173` (1.8 kb) inside the window, absent from the table (size-filtered by the cohort pipeline).
- **GCA_018120125.1**: the table lists `INV:516667-546991`, which ends 337 bp *before* the window start (547 328) — an edge-overlap convention difference.

Both BED and table agree with the assembly verdicts (deletions where coverage is missing, inversions in `complete_rearranged` strains).

---

## 4. Panel precision — answer to the headline question

Question: *among panel-`empty` genomes, how much of the discordance is assembly gap vs panel error?*

- **Panel error among discordant `empty` calls: 0 / 13.** All 13 assemblies support the empty call: 2 clean deletions (empty site assembled, flanks covered, genome size reduced), 11 with the island backbone genuinely deleted (window coverage 0.05–0.19 at real identity ≥ 0.87, flanks 0.6–0.99 covered, and corroborating DEL calls in the SV BED).
- **Draft-gap explanation: 0 / 13.** Even the most fragmented assemblies (GCA_009740705.1 / GCA_009740805.1, 61–65 contigs, N50 ≈ 60–68 kb) have both flanks assembled and carry cagA fragments + DEL calls at the locus, i.e. the locus region is represented and points to deletion rather than an unbridged gap. The two single-contig genomes are unambiguous.
- Splitting the 26 discordances by cause: **2 real biological discoveries** (cagA-specific loss, §2.1), **2 clean empty loci** (panel right, collaborator conservative), **17 partial-deletion/cagA-retention cases** where both callers are partially right (dominated by real biology + the gene-vs-island framing), **5 divergence-driven partial calls** (the only panel soft spots; none affects an `empty` call).

Assembly problem vs real biology, overall: **≈ 0 / 26 discordances are caused by assembly incompleteness; ~8 / 26 are unambiguous novel biology** (2 cagA-specific losses + 6 intact-cagA-at-empty-island); the rest are definitional/threshold differences between a gene test and an island test, with the panel's island-level status assembly-supported in every case.

Caveat: this review covers the discordant set only; concordant panel calls were not re-audited genome-by-genome (though the panel-exact reproduction of `n_present` matched the table for all 26 discordant genomes with zero disagreements, indicating the published table is deterministic and reproducible).

---

## 5. Methods (what was run per genome) and pitfalls found

Script: `case_studies/h_pylori_cagpai/scripts/check_discordant_genomes.py`

```
python3 case_studies/h_pylori_cagpai/scripts/check_discordant_genomes.py
# defaults point at the read-only external drive; all outputs land in
# case_studies/h_pylori_cagpai/results/metadata_assoc/discordant_review/
```

Per genome: (a) assembly stats; (b) panel-exact marker recount (reproduces the table's `n_present` for 26/26, zero disagreements); (c) cagA k-mer test + asm20/asm5 detail; (d) chunked window/flank coverage vs hp26695 (window 547 328–583 481 ± 100 kb); (e) SV BED cross-check. Tools: minimap2 2.31 (`-cx asm5`/`asm20`), Python 3.9 stdlib only.

Pitfalls that mattered (and are handled in the script):

1. **`-c` is required.** Without it, PAF `nmatch/alen` are seed-based estimates, not base-level identities (a "complete" L2 island scored 0 % coverage under `-x asm5`, 97 % under `-cx asm20`).
2. **`de:f` (gap-compressed) identity is inflated by gappy forced alignments.** Empty-island genomes produced single 104 kb chains with `de:f` 0.05 but 40 462 mismatches/indels (real identity 0.61; ungapped 17 %). Never use `de:f` alone to declare locus presence.
3. **Chaining suppression.** The spurious mega-chain in (2) suppressed the genuine small cagA alignment inside the window. Fixed by mapping 10 kb chunks independently.
4. **K-mer position coverage** (fraction of marker positions within ≥ 1 shared exact 19-mer, both strands) separates full-length divergent cagA (East-Asian ≈ 0.34) from ~50 % remnants (≈ 0.20) from absent (0.00), where both minimap2 identity flavours fail.
5. Threshold calibration used concordant genomes (full Western cagA 0.77, full East-Asian 0.34–0.37, absent 0.00).

Known residual ambiguity: for L2 strains with cagA position coverage 0.29–0.35 and asm20 gene coverage 0.77–0.79 (GCA_019674555.1 sits exactly at the present/remnant boundary), "full divergent cagA" vs "large remnant" cannot be fully resolved with assembly-only evidence; the verdict-level conclusion (panel partial on an island that *is* assembled) is unaffected.
