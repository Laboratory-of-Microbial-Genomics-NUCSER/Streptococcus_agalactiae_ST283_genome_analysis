import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, fisher_exact
from statsmodels.stats.multitest import multipletests
from itertools import combinations
import warnings

# =======================================================
# STEP 1 — LOAD DATA
# =======================================================
df = pd.read_csv(
    r"C:\Users\Yajnesh\Documents\Correlation\AMR_VGs_ST283_binary.tsv",
    sep="\t",
    index_col=0
)

# Ensure numeric
df = df.apply(pd.to_numeric, errors="coerce")

# =======================================================
# PREVALENCE FILTER
# Remove genes present in fewer than 3 isolates
# =======================================================
min_presence = 3

df = df.loc[:, df.sum(axis=0) >= min_presence]

print("After prevalence filtering:", df.shape)
print("Data shape:", df.shape)

# =======================================================
# STEP 2 — PHI CORRELATION MATRIX
# =======================================================
corr = df.corr(method="pearson")

corr.to_csv(
    r"C:\Users\Yajnesh\Documents\Correlation\correlation_matrix.tsv",
    sep="\t"
)

print("✅ Correlation matrix saved")

# =======================================================
# STEP 3 — PAIRWISE FISHER EXACT TEST
# =======================================================
genes = corr.columns.tolist()

# Empty matrix for Fisher p-values
pvals_fisher = pd.DataFrame(
    np.nan,
    index=genes,
    columns=genes
)

# Store pairwise p-values
pval_list = []
pairs = []

for g1, g2 in combinations(genes, 2):

    # Skip invariant genes
    if df[g1].nunique() < 2 or df[g2].nunique() < 2:
        continue

    # 2x2 contingency table
    table = pd.crosstab(df[g1], df[g2])

    # Force full 2x2 structure
    table = table.reindex(
        index=[0, 1],
        columns=[0, 1],
        fill_value=0
    )

    # Fisher exact test
    _, p = fisher_exact(table)

    pval_list.append(p)
    pairs.append((g1, g2))

# =======================================================
# STEP 4 — FDR CORRECTION
# =======================================================
p_adj = multipletests(
    pval_list,
    method="fdr_bh"
)[1]

# Map adjusted p-values back
for (g1, g2), p_corr in zip(pairs, p_adj):

    pvals_fisher.loc[g1, g2] = p_corr
    pvals_fisher.loc[g2, g1] = p_corr

# Diagonal
np.fill_diagonal(pvals_fisher.values, 0)

# Save
pvals_fisher.to_csv(
    r"C:\Users\Yajnesh\Documents\Correlation\fisher_FDR_pvalues.tsv",
    sep="\t"
)

print("✅ Fisher exact test + FDR completed")

# =======================================================
# STEP 5 — FILTERED LOWER TRIANGLE DATA
# =======================================================
x, y, vals = [], [], []

cutoff = 0.25

n = len(corr)

for i in range(n):
    for j in range(n):

        if i > j:   # lower triangle only

            r = corr.iloc[i, j]
            p = pvals_fisher.iloc[i, j]

            # Skip missing p-values
            if pd.isna(p):
                continue

            # FILTER:
            # |Phi| >= 0.25
            # AND FDR-adjusted Fisher p < 0.05
            if abs(r) >= cutoff and p < 0.05:

                x.append(j)
                y.append(i)
                vals.append(r)

# =======================================================
# STEP 6 — PLOT
# =======================================================
plt.figure(figsize=(16, 14), dpi=300)

ax = plt.gca()
ax.set_axisbelow(True)

# -------------------------------------------------------
# GRIDLINES (LOWER TRIANGLE ONLY)
# -------------------------------------------------------
for i in range(n):
    for j in range(n):

        if i >= j:

            ax.plot(
                [j, j],
                [i, i+1],
                linestyle="--",
                linewidth=0.5,
                alpha=0.6,
                zorder=0
            )

            ax.plot(
                [j, j+1],
                [i+1, i+1],
                linestyle="--",
                linewidth=0.5,
                alpha=0.6,
                zorder=0
            )

# -------------------------------------------------------
# BUBBLE PLOT
# -------------------------------------------------------
sc = ax.scatter(
    x,
    y,
    s=[abs(v) * 900 for v in vals],
    c=vals,
    cmap="coolwarm",
    vmin=-1,
    vmax=1,
    edgecolor="black",
    linewidth=0.3,
    alpha=0.9,
    zorder=3
)

# -------------------------------------------------------
# TEXT LABELS
# -------------------------------------------------------
for xi, yi, v in zip(x, y, vals):

    ax.text(
        xi,
        yi,
        f"{v:.2f}",
        ha="center",
        va="center",
        fontsize=8,
        fontweight="bold"
    )

# -------------------------------------------------------
# AXES
# -------------------------------------------------------
plt.xticks(
    range(n),
    corr.columns,
    rotation=60,
    fontsize=9,
    fontweight="bold"
)

plt.yticks(
    range(n),
    corr.index,
    fontsize=9,
    fontweight="bold"
)

ax.invert_yaxis()
ax.set_aspect("equal")

plt.title(
    "Pairwise Correlation of Gene Presence–Absence",
    fontsize=14,
    fontweight="bold"
)

plt.tight_layout()

# =======================================================
# SAVE FIGURES
# =======================================================
plt.savefig(
    r"C:\Users\Yajnesh\Documents\Correlation\correlation_plot_filtered.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    r"C:\Users\Yajnesh\Documents\Correlation\correlation_plot_filtered.svg"
)

plt.show()