# H. pylori cagPAI case study

This case study demonstrates how Syn2bANI combines fast ANI search with structural profiling using a clinically important pathogenicity island.

## Data

- 528 H. pylori genomes from Song et al. (2026), located in the paper data area.
- Reference strain 26695 (NC_000915.1) with annotated cagPAI locus tags HP0520–HP0547.
- Isolate metadata include disease stage (NAG/AG/IM/GC), country, FastBAPS lineage, and phylogenetic population.

## Scripts

- `call_cagpai_status.py` — classify each genome as `empty` / `partial` / `complete` based on 28 cagPAI marker genes.
- `run_struct_vs_26695.py` — run `syn2bani struct --bed` for every cohort genome vs 26695 in parallel.
- `classify_rearrangement.py` — extend the marker-based state with structural calls inside the cagPAI coordinates.
- `associate_cagpai_metadata.py` — correlate extended cagPAI states with metadata and generate stacked-bar figures.

## Results

See `results/struct_extended/`:

- `cagpai_states_extended.tsv` — per-genome extended state (`empty`, `partial`, `complete_collinear`, `complete_rearranged`).
- `cagpai_association.tsv` — contingency tables and χ² tests by group, country, FastBAPS, and population.
- `cagpai_summary_extended.md` — interpretation of the findings.
- `cagpai_state_by_*.png` — publication-style stacked-bar visualizations.

## Key findings

- 79.6% of strains that retain cagPAI show a detectable rearrangement (inversion or translocation) relative to 26695.
- Only 2.3% are fully collinear.
- cagPAI state is significantly associated with disease stage (χ² p = 9.9 × 10⁻⁴) and FastBAPS lineage (χ² p = 3.9 × 10⁻⁸).
- Engineered controls (wt, hp26695, mut1, inv, transloc, mut1_inv, mut1_transloc, del, mut1_del) are classified consistently with their known architectures.
