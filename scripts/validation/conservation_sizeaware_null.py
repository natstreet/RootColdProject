#!/usr/bin/env python3
"""Family-size-aware null for the conserved co-expression core (companion to the label-permutation null).

The published null (conservation_label_permutation.R) replaces each comparison's co-expressolog
orthogroups by a uniform random draw of the same size, so every orthogroup is equally likely to be
drawn. Large orthogroups contain more gene pairs and are therefore more often co-expressologs in any
one comparison. This null keeps that dependence: within each of 15 comparisons, an orthogroup is
included with the observed inclusion rate of orthogroups of its size class in that comparison, and
comparisons are independent. The statistic is the number of orthogroups included in all 15.

Inputs (deposit v4 or later, unpacked into ./data): data/open/co_expressologs.tsv, data/*/Orthogroups_20240613.tsv
Outputs: results/conservation_sizeaware_null.csv (null_core_size), results/conservation_sizeaware_summary.tsv
Usage:   python3 scripts/validation/conservation_sizeaware_null.py
"""
import glob, sys
import numpy as np, pandas as pd
P, N_SIM, SEED = 0.1, 5000, 1
EDGES = [2, 6, 10, 15, 20, 30, 50, 75, 100, 150, 250, 500, 10**6]
SP = ["Arabidopsis_thaliana", "Populus_tremula", "Betula_pendula", "Picea_abies", "Pinus_sylvestris"]
co = pd.read_csv("data/open/co_expressologs.tsv", sep="\t"); co = co[co.MaxpVal < P]
co["cmp"] = [" ".join(sorted(x)) for x in zip(co.Species1, co.Species2)]
per = co.groupby("cmp").OrthologGroup.apply(set); assert len(per) == 15
og = pd.read_csv(glob.glob("data/*/Orthogroups_20240613.tsv")[0], sep="\t", dtype=str)
n = og[SP].fillna("").apply(lambda c: c.str.count(",") + (c.str.len() > 0).astype(int))
keep = (n > 0).all(axis=1); U = og.loc[keep, "Ortholog_Group"].values; size = n.loc[keep].sum(axis=1).values
M = np.array([[u in per[c] for u in U] for c in per.index])
obs = int(M.all(axis=0).sum())
b = np.digitize(size, EDGES); R = np.zeros(M.shape)
for k in np.unique(b): R[:, b == k] = M[:, b == k].mean(axis=1, keepdims=True)
rng = np.random.default_rng(SEED)
null = np.array([int((rng.random(R.shape) < R).all(axis=0).sum()) for _ in range(N_SIM)])
pd.DataFrame({"null_core_size": null}).to_csv("results/conservation_sizeaware_null.csv", index=False)
out = pd.DataFrame({"quantity": ["observed", "universe", "expected", "null_mean", "null_max", "fold", "empirical_P", "n_sim"],
                    "value": [obs, len(U), R.prod(axis=0).sum(), null.mean(), null.max(), obs / R.prod(axis=0).sum(),
                              ((null >= obs).sum() + 1) / (N_SIM + 1), N_SIM]})
out.to_csv("results/conservation_sizeaware_summary.tsv", sep="\t", index=False); print(out.to_string(index=False))
