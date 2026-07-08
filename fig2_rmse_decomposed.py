#!/usr/bin/env python3
# ============================================================
# fig2_rmse_decomposed.py
#
# Generates Figure 2: RMSE decomposed by true-association status
# (true-nonzero signal pairs vs true-zero null pairs)
#
# UPDATE SCENARIO 2 VALUES when new simulation run completes.
# Scenarios 1 and 3 are final (n=100 replicates).
#
# Run: python fig2_rmse_decomposed.py
# Output: fig2_rmse_decomposed.png (300 DPI)
# ============================================================

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Global style ─────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':       'DejaVu Sans',
    'font.size':         11,
    'axes.titlesize':    12,
    'axes.labelsize':    11,
    'xtick.labelsize':   10,
    'ytick.labelsize':   10,
    'savefig.dpi':       300,
    'savefig.bbox':      'tight',
    'axes.spines.top':   False,
    'axes.spines.right': False,
})

# ── Colours (colour-blind friendly) ─────────────────────────────────────
COL = {
    'full':         '#0072B2',   # blue
    'no_network':   '#E69F00',   # orange
    'no_horseshoe': '#CC79A7',   # pink
}
LAB = {
    'full':         'Full model',
    'no_network':   'No-network',
    'no_horseshoe': 'No-horseshoe',
}
MODELS = ['full', 'no_network', 'no_horseshoe']

# ── Data ─────────────────────────────────────────────────────────────────
# Each entry: {'signal': RMSE on true nonzero pairs,
#              'null':   RMSE on true zero pairs}
#
# UPDATE sc2 values below when Scenario 2 rerun completes.
# Replace the four numbers marked with  <-- UPDATE

data = {
    # ── Scenario 1: Sparse Independent (FINAL, n=100) ────────────────
    'sc1': {
        'full':         {'signal': 0.562, 'null': 0.066},
        'no_network':   {'signal': 0.263, 'null': 0.102},
        'no_horseshoe': {'signal': 0.415, 'null': 0.166},
    },
    # ── Scenario 2: Network Structured (UPDATE WHEN RERUN COMPLETES) ─
    'sc2': {
        'full':         {'signal': 0.586, 'null': 0.067},   # <-- UPDATE
        'no_network':   {'signal': 0.263, 'null': 0.099},   # <-- UPDATE
        'no_horseshoe': {'signal': 0.418, 'null': 0.166},   # <-- UPDATE
    },
    # ── Scenario 3: Dense Weak Noise (FINAL, n=100) ───────────────────
    'sc3': {
        'full':         {'signal': 0.562, 'null': 0.073},
        'no_network':   {'signal': 0.259, 'null': 0.104},
        'no_horseshoe': {'signal': 0.411, 'null': 0.167},
    },
}

SC_LABELS = {
    'sc1': 'Sc.1: Sparse\nIndependent',
    'sc2': 'Sc.2: Network\nStructured',
    'sc3': 'Sc.3: Dense\nWeak Noise',
}

# ── Plot ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(13, 5))
fig.suptitle(
    'Figure 2  |  RMSE Decomposed by True-Association Status (n = 100 replicates)',
    fontweight='bold', y=1.01
)

x     = np.arange(len(MODELS))
width = 0.35

for col, sc in enumerate(['sc1', 'sc2', 'sc3']):
    ax = axes[col]

    sig_vals  = [data[sc][m]['signal'] for m in MODELS]
    null_vals = [data[sc][m]['null']   for m in MODELS]
    colors    = [COL[m] for m in MODELS]

    # Solid bars = true nonzero (signal) pairs
    bars_sig = ax.bar(x - width / 2, sig_vals, width,
                      color=colors, alpha=0.87)

    # Hatched bars = true zero (null) pairs
    bars_null = ax.bar(x + width / 2, null_vals, width,
                       color=colors, alpha=0.40, hatch='//')

    # Annotate signal bars (bold)
    for bar, v in zip(bars_sig, sig_vals):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.005,
                f'{v:.3f}',
                ha='center', va='bottom',
                fontsize=8.5, fontweight='bold')

    # Annotate null bars (regular weight)
    for bar, v in zip(bars_null, null_vals):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.005,
                f'{v:.3f}',
                ha='center', va='bottom',
                fontsize=8.5)

    ax.set_title(SC_LABELS[sc], fontweight='bold', pad=8)
    ax.set_xticks(x)
    ax.set_xticklabels(['Full\nmodel', 'No-\nnetwork', 'No-\nhorseshoe'],
                       fontsize=9.5)
    ax.set_ylabel('RMSE' if col == 0 else '')
    ax.set_ylim(0, 0.70)
    ax.grid(axis='y', alpha=0.25, zorder=0)

# ── Legends ──────────────────────────────────────────────────────────────
# Legend 1: solid vs hatched (subset type)
from matplotlib.patches import Patch
subset_legend = [
    Patch(facecolor='grey', alpha=0.87,
          label='True nonzero pairs (3.3% of all pairs)'),
    Patch(facecolor='grey', alpha=0.40, hatch='//',
          label='True zero pairs (96.7% of all pairs)'),
]
leg1 = fig.legend(handles=subset_legend,
                  loc='lower center', ncol=2, frameon=False,
                  bbox_to_anchor=(0.5, -0.10), fontsize=10)

# Legend 2: model colours
model_handles = [mpatches.Patch(color=COL[m], label=LAB[m], alpha=0.87)
                 for m in MODELS]
fig.legend(handles=model_handles,
           loc='lower center', ncol=3, frameon=False,
           bbox_to_anchor=(0.5, -0.22), fontsize=10)

fig.add_artist(leg1)   # keep both legends visible

# ── Save ─────────────────────────────────────────────────────────────────
plt.tight_layout(rect=[0, 0.12, 1, 1])
fig.savefig('fig2_rmse_decomposed.png', dpi=300, bbox_inches='tight')
plt.close()

print("Saved: fig2_rmse_decomposed.png")
print()
print("To update Scenario 2:")
print("  Edit the three lines marked  # <-- UPDATE  in the data dict above,")
print("  then re-run:  python fig2_rmse_decomposed.py")
