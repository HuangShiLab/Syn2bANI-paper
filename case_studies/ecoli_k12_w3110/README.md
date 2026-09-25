# Case study: E. coli K-12 MG1655 × W3110 (the ANI-blind flagship vignette)

**Question.** What do ANI-only tools report for the textbook E. coli K-12 pair,
and what does Syn2b/Syn2bANI add?

**Pair.**
- MG1655: `ecoli_k12_MG1655.fasta` (NC_000913.3, 4,641,652 bp)
- W3110: `genomes/GCF_000010245.2.fna` (NC_007779.1, 4,646,332 bp), fetched with
  `case_studies/fetch_ncbi_fna.py`

**Headline numbers.**

| tool | readout |
|---|---|
| skani 0.3.2 | ANI 99.99, AF 100/100 — no structural readout |
| minimap2 asm20 (ANIm ≈ 99.6) | one inverted block, 3,423,516–4,216,755 (793 kb, 17.1% of the chromosome) |
| Syn2b `synteny` | breakpoints = 2 (exactly the two inversion boundaries), inverted fraction = 0.160, junctions at 3,423,157 / 4,207,508 |
| Syn2bANI `struct --bed --circular NC_000913.3` | 11 SVs: INV 3,428,880–4,211,231 (782 kb) + 9 short indels |

The 782-kb inversion is the classic W3110/MG1655 inversion between rrnD and rrnE
(Hill & Harnish 1981; Hill et al. 1990), covering ~700 genes. Syn2b's two
junctions sit within 5.4 kb of the alignment-derived boundaries; the inverted
fraction (0.160) matches the alignment-based inverted length (0.171).

**Story.** An ANI-only search reports these two strains as near-identical
(99.99%, 100% aligned). The structural axis reports a 782-kb inversion — a
documented source of phenotypic divergence between the two most-cited E. coli
strains. At ~19 ms per pair, the second axis is free.

**Reproduce.**

```bash
python3 ../fetch_ncbi_fna.py --accessions accessions.txt --outdir genomes/
skani dist genomes/ecoli_k12_MG1655.fasta genomes/GCF_000010245.2.fna
minimap2 -cx asm20 genomes/ecoli_k12_MG1655.fasta genomes/GCF_000010245.2.fna > out/minimap2.paf
syn2b digest --enzymes BcgI,AlfI,AloI,FalI -i genomes/ecoli_k12_MG1655.fasta -o out/mg1655.tgt
syn2b digest --enzymes BcgI,AlfI,AloI,FalI -i genomes/GCF_000010245.2.fna -o out/w3110.tgt
mkdir -p tgts && cp out/*.tgt tgts/
syn2b synteny -i tgts -o out/syn2b_synteny.tsv
syn2bani struct genomes/ecoli_k12_MG1655.fasta genomes/GCF_000010245.2.fna --bed --circular NC_000913.3 > out/struct_bed.tsv
python3 make_figure.py
```

**Figure.** `figures/fig_k12_w3110_vignette.{png,pdf}` (panel a: dotplot with
the inversion and Syn2b junctions; panel b: same pair, three readouts).
