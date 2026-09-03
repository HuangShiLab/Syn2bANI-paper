#!/usr/bin/env python3
"""Generate ANI vs breakpoint scatter plots for E. coli O157:H7 and S. aureus case studies."""
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path('/Users/macstudio/Downloads/Syn2bANI-paper/case_studies')
OUT = Path('/Users/macstudio/Downloads/Syn2bANI-paper/paper/supplementary/figures')
OUT.mkdir(parents=True, exist_ok=True)

# E. coli O157:H7
eco = BASE / 'ecoli_o157_fitzgerald_2021'
tri_eco = pd.read_csv(eco / 'results' / 'triangle.tsv', sep='\t')
meta_eco = pd.read_csv(eco / 'results' / 'metadata_with_lineage.tsv', sep='\t')
meta_eco = meta_eco.rename(columns={'nucleotide_acc': 'query', 'assigned_lineage': 'lineage'})

def norm_acc(x):
    if pd.isna(x):
        return x
    x = str(x).split('.')[0]
    if x.startswith('NZ_'):
        x = x[3:]
    return x

tri_eco['query'] = tri_eco['query'].apply(norm_acc)
meta_eco['query'] = meta_eco['query'].apply(norm_acc)
eco_merged = tri_eco.merge(meta_eco[['query', 'lineage', 'host_category']], on='query', how='left')

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, col, title in zip(axes, ['lineage', 'host_category'], ['Lineage', 'Host category']):
    for val, grp in eco_merged.groupby(col):
        ax.scatter(grp['ani'], grp['breakpoint_count'], label=val, alpha=0.6, s=20)
    ax.set_xlabel('ANI (%)')
    ax.set_ylabel('Breakpoint count')
    ax.set_title(f'E. coli O157:H7 by {title}')
    ax.legend(title=title, fontsize=7)
plt.tight_layout()
fig.savefig(OUT / 'Fig_S15_ecoli_o157_ani_vs_breakpoints.png', dpi=300)
plt.close(fig)

# S. aureus
sau = BASE / 'fda_argos_s_aureus'
tri_sau = pd.read_csv(sau / 'results' / 'triangle.tsv', sep='\t')

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].scatter(tri_sau['ani'], tri_sau['breakpoint_count'], alpha=0.4, s=10, c='steelblue')
axes[0].set_xlabel('ANI (%)')
axes[0].set_ylabel('Breakpoint count')
axes[0].set_title('S. aureus: ANI vs breakpoint count')

axes[1].hist(tri_sau['breakpoint_count'], bins=60, color='steelblue', edgecolor='white')
axes[1].set_xlabel('Breakpoint count')
axes[1].set_ylabel('Number of pairs')
axes[1].set_title('S. aureus: breakpoint count distribution')
plt.tight_layout()
fig.savefig(OUT / 'Fig_S16_saureus_ani_vs_breakpoints.png', dpi=300)
plt.close(fig)

print('Wrote', OUT / 'Fig_S15_ecoli_o157_ani_vs_breakpoints.png')
print('Wrote', OUT / 'Fig_S16_saureus_ani_vs_breakpoints.png')
