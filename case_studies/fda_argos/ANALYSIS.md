# Structural discordance within ANI-identical FDA-ARGOS *E. coli* pairs

Analysis date: 2026-09-25. Pair set: 103 pairs with skani **0.3.2** ANI ≥ 99.9
(`pairs99.9_ecoli.tsv`, regenerated from `skani_triangle_ecoli_v032.tsv`).
Tool: local `syn2bani struct` v0.1.0 (`/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani`),
per pair: `syn2bani struct <query>.fna <ref>.fna --bed --circular <all ref contigs> -o struct_pairs99.9/<ref>__<query>.bed`
(reference = alphabetically first accession; query = second). All 103 runs exit 0.

## Headline (post-validation; see "Validation against minimap2" below)

- **101 of 103 ANI ≥ 99.9 pairs (98%) carry ≥ 1 structural-variant call — and ≥ 1
  *genuine* structural difference after minimap2 validation** (876 indels ≥ 1 kb
  cs-verified, 445 strand-confirmed rearrangements, 880 plasmid-content events).
- Total raw struct events: 2,718 — INV 909, INS 707, DEL 661, TRA 441 — of which
  **358 (13%) are orientation/origin artifacts** (mostly the largest INV/TRA calls);
  see validation section before quoting size statistics.
- Per pair: median 25 raw calls (range 0–61).
- Even at ANI ≥ 99.99 (24 pairs), 22 carry ≥ 1 SV.
- Only 2 pairs are structurally silent — both are ANI 100.0 pairs from the
  GCA_01904xxxxx re-assembly batch (FDAARGOS_1380/1381, FDAARGOS_1378/1379).

## Five most extreme raw pairs (by summed affected bp) — VALIDATED

| pair (ref vs query) | ANI | n SV | raw bp | validation verdict |
|---|---|---|---|---|
| GCF_016903335.1 vs GCF_016904195.1 | 99.93 | 29 | 8.1 Mb | **flip/rotation artifact** (mf=0.981); 10 artifact INV/TRA, but 11 cs-verified indels remain |
| GCF_016903835.1 vs GCF_016904075.1 | 99.90 | 28 | 6.5 Mb | **flip/rotation artifact** (mf=0.983); 12 artifact calls, 15 cs-verified indels remain |
| GCF_016904195.1 vs GCF_016904535.1 | 99.92 | 31 | 6.2 Mb | **flip/rotation artifact** (mf=0.974); 12 artifact calls, 8 cs-verified indels remain |
| GCF_002206405.2 vs GCF_017081115.1 | 99.93 | 58 | 5.8 Mb | **genuine**: partial inversion, 23 strand-confirmed rearrangements + 21 cs-verified indels |
| GCF_002208865.2 vs GCF_017081335.1 | 99.91 | 61 | 5.5 Mb | **genuine**: partial inversion, 23 strand-confirmed rearrangements + 27 cs-verified indels |

## Metadata concordance (from `metadata_ecoli.tsv`)

- The USA:WA cluster extreme pairs (top 3 above) are flip/rotation artifacts at the Mb
  scale, but still differ in AMR cassette content (`aadA2, dfrA12, mph(A), mrx(A),
  ompC_Q76Ter, sul1`) consistent with their 8–15 cs-verified indels.
- The two genuine extreme pairs (FDAARGOS_292/293 vs the FDAARGOS_1238/1355 cluster)
  differ in `catA1`/`emrD` carriage and geography (USA:DC vs USA:MD).
- Strongest verified flagships (all both-complete, collinear): see `validate/VALIDATION.md` —
  e.g. FDAARGOS_1264/1265 (ANI 99.93) with a verified **blaCTX-M-15/blaOXA-1/aac(6')-Ib-cr5**
  cassette difference plus four 38–51 kb prophage indels; FDAARGOS_348/772 at ANI **100.0**
  differing by **tet(A)/tet(K)** and a 48 kb insertion.
- Take-home for the paper: within-clone ANI ≥ 99.9 (even ≥ 99.99) does **not** imply
  structural identity — prophage turnover, plasmid and AMR-cassette gains/losses, and
  (in 16/103 pairs) genuine large inversions are routine. Large inversion *calls*
  additionally require strand-convention control, as this validation shows.

## Validation against minimap2 (2026-09-25)

Every pair was re-aligned with `minimap2 -cx asm20` (+ `--cs` for indel verification);
full details in `validate/VALIDATION.md`, per-pair numbers in `validate/validation_summary.tsv`.

- Pair classes (strand split on the reference chromosome): **41 global_flip_or_rotation**
  (minus_fraction > 0.8 — whole-chromosome reverse-complement/rotation conventions),
  **16 partial_inversion** (0.2–0.8 — genuine large inversions), **46 collinear_local**
  (≤ 0.2). Manual spot-check reproduced: FDAARGOS_401 vs 772 aligns in exactly two
  minus-strand blocks abutting at the struct INV breakpoint.
- Event accounting: of 2,718 raw calls, **358 (13%) are artifacts** (355 INV/TRA inside
  flipped background or at rotation junctions + 3 junction indels). Artifacts hold most
  of the *span*: raw "affected bp" sums up to 8.1 Mb are convention artifacts.
- Genuine content: 876 chromosomal indels ≥ 1 kb verified directly in PAF cs strings,
  445 strand-confirmed/uncontradicted rearrangements (all 16 partial-inversion pairs
  confirmed), 880 off-chromosome plasmid/content events, 159 small unverified indels.
- **Corrected headline: 101/103 pairs (98%) carry ≥ 1 genuine structural difference**;
  41/41 flipped pairs retain genuine indel/plasmid differences; only the 2 ANI-100
  re-assembly pairs are fully silent.

## Caveats

- `total_affected_bp` in `struct_summary.tsv` sums raw per-event reference spans
  (BED end−start) **including artifact-zone calls** — use `validate/validation_summary.tsv`
  (`n_rearr_artifact`, `n_indel_cs_verified`, `n_sv_local_genuine_estimate`) for
  validated numbers. Overlapping events also double-count shared bases.
- The `--circular` origin filter (default threshold 0.5 of contig length) does **not**
  catch rotation artifacts presenting below 50% of the contig; 41/103 pairs (40%) are
  global flip/rotation conventions (minimap2 minus_fraction > 0.8). struct
  rearrangement calls on such chromosomes need strand-convention control before
  biological interpretation; indel calls are largely unaffected (392/454 cs-verified
  in flipped pairs).
- cs-verification threshold: indels ≥ 1 kb on the chromosome; smaller indels (159
  events) are "likely genuine, unverified". 17 pairs have >30 alignment blocks
  (contig-level assemblies) where TRA calls remain partly unresolvable.
- BED score column = anchor support (left+right), not event size; see
  `struct_summary.tsv` and per-pair BEDs in `struct_pairs99.9/`.
