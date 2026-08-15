#!/usr/bin/env python3
# make_FigS5.py — Figure S5: robustness of the conserved co-expression core.
#   (a) Conservation permutation test: the observed all-pairwise core size against the
#       degree-preserving rewired null.
#   (b) Cluster coherence in independent data (STRING): observed within-cluster co-expression
#       edges among the core Arabidopsis genes against the label-permuted null.
#
# This step only plots. The null distributions and observed values are produced by the two
# analysis scripts in scripts/validation/, which require inputs that are NOT in the figshare
# deposit and must be obtained separately:
#   panel (a): scripts/validation/conservation_null_test.py  -> conservation_null_null.csv,
#              conservation_null_obs.txt   (needs the per-species co-expression edge lists)
#   panel (b): scripts/validation/string_external_validation.py -> string_null.csv,
#              string_obs.txt   (needs the STRING v12 A. thaliana co-expression download)
# Run those first, then run this from the same directory.
#
# Output: Figure_S5.pdf, Figure_S5.png

import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

nullA = pd.read_csv("conservation_null_null.csv")["null_core_size"].values
obsA  = int(open("conservation_null_obs.txt").read().strip())
nullB = pd.read_csv("string_null.csv")["null"].values
obsB  = int(open("string_obs.txt").read().strip())
pA = (np.sum(nullA >= obsA) + 1) / (len(nullA) + 1)
pB = (np.sum(nullB >= obsB) + 1) / (len(nullB) + 1)
foldB = obsB / max(nullB.mean(), 1e-9)

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
ax[0].hist(nullA, bins=np.arange(-0.5, max(obsA, nullA.max()) + 2), color="#add8e6", edgecolor="none")
ax[0].axvline(obsA, color="#d62728", lw=2.5)
ax[0].annotate(f"observed = {obsA}", xy=(obsA, 5), xytext=(obsA * 0.55, len(nullA) * 0.62),
               color="#d62728", fontsize=13, va="center",
               arrowprops=dict(color="#d62728", arrowstyle="->", lw=1.5))
ax[0].set_title(f"(a) Conservation permutation test (null max = {nullA.max()}; P = {pA:.4g})", fontsize=12)
ax[0].set_xlabel("Orthogroups conserved across all 15\ncomparisons (permuted null)", fontsize=11)
ax[0].set_ylabel("Permutations", fontsize=11); ax[0].set_xlim(-3, obsA + 6)

ax[1].hist(nullB, bins=30, color="#90ee90", edgecolor="white")
ax[1].axvline(obsB, color="#d62728", lw=2.5)
h = np.histogram(nullB, bins=30)[0].max()
ax[1].annotate(f"observed = {obsB}", xy=(obsB, h * 0.02), xytext=(obsB * 0.60, h * 0.62),
               color="#d62728", fontsize=13, va="center",
               arrowprops=dict(color="#d62728", arrowstyle="->", lw=1.5))
ax[1].set_title(f"(b) Cluster coherence in independent data ({foldB:.1f}-fold; P = {pB:.4g})", fontsize=12)
ax[1].set_xlabel("Within-cluster co-expression edges in\nSTRING (clusters randomised)", fontsize=11)
ax[1].set_ylabel("Permutations", fontsize=11); ax[1].set_xlim(nullB.min() - 20, obsB + 30)

plt.tight_layout()
plt.savefig("Figure_S5.pdf"); plt.savefig("Figure_S5.png", dpi=150)
print(f"wrote Figure_S5.pdf/.png ; panel a P={pA:.4g}, panel b P={pB:.4g}, fold={foldB:.1f}")
