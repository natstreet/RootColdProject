#!/usr/bin/env python3
# make_FigS7.py — Figure S7: worked example of resource use, conserved-core orthogroup
#   OG0000118 (NAC transcription factors, ATAF/NAM clade).
#   (a) Cold response of the orthogroup's member gene in each genotype across the time
#       course, as log2 fold change against the matched control; an open ring marks a timepoint
#       at which it is differentially expressed (adjusted P <= 0.05 and |log2FC| >= 1).
#   (b) The six members and their co-expressolog relationships: all 15 pairwise species
#       comparisons are linked, so the orthogroup forms a complete conserved clique.
#   (c) Number of genes co-expressed with each member in that species' whole-genome
#       co-expression network, with the transcription-factor subset overlaid.
#
# Inputs, all from the figshare deposit (version 4 or later; the open-format tables and
# the per-species networks were introduced in version 4). Run from the directory that
# contains data/ :
#   data/open/conserved_core.tsv          the 116 conserved-core orthogroups
#   data/open/gene_index.tsv              per-gene fold changes, adjusted P, TF family
#   data/open/co_expressologs.tsv         co-expressolog pairs, used for panel (b)
#   data/open/transcription_factors.tsv   TF assignments, used for panel (c)
#   data/networks/network_<species>_hrr100.tsv   per-species co-expression networks
#   DE/<species>/DEGs_<species>_C_vs_<timepoint>.txt   per-timepoint DEG lists written by
#                                         scripts/DE/differential_expression.R (the rings in panel a)
#
# Output: Figure_S7.pdf, Figure_S7.png
#
# ORTHOGROUP IDENTIFIERS ARE RUN-SPECIFIC. OG0000118 refers to the orthology in
# data/orthologs/; it is not comparable with identifiers from any other OrthoFinder run.

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OG = "OG0000118"
SP = ["col0", "ost0", "aspen", "birch", "spruce", "pine"]
LAB = {"col0": "Col-0", "ost0": "Ost-0", "aspen": "aspen",
       "birch": "birch", "spruce": "spruce", "pine": "pine"}
COL = {"col0": "#0072B2", "ost0": "#56B4E9", "aspen": "#009E73",
       "birch": "#7FBF3F", "spruce": "#D55E00", "pine": "#CC79A7"}
TPS = ["6h", "24h", "3d", "10d"]
GREY = "#666666"

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "figure.dpi": 110,
})

core = pd.read_csv("data/open/conserved_core.tsv", sep="\t")
gi = pd.read_csv("data/open/gene_index.tsv", sep="\t")
coex = pd.read_csv("data/open/co_expressologs.tsv", sep="\t")
tf = pd.read_csv("data/open/transcription_factors.tsv", sep="\t")

row = core.loc[core.OrthoGroup == OG].iloc[0]
mem = {s: row[s] for s in SP if pd.notna(row[s])}
tfset = {(r.species, r.gene) for r in tf.itertuples()}

# Differentially expressed at a timepoint as defined in the Methods (Wald test against no change,
# adjusted P <= 0.05 and |log2FC| >= 1): membership of that timepoint's DEG list from
# differential_expression.R. (The deposit's padj columns test against a two-fold change instead.)
_deg_cache = {}
def is_deg(species, gene, tp):
    key = (species, tp)
    if key not in _deg_cache:
        hits = [f for f in os.listdir(f"DE/{species}") if f.startswith(f"DEGs_{species}_C_vs_{tp}")]
        if len(hits) != 1:
            raise SystemExit(f"expected one DEG list for {species} {tp} in DE/{species}, found {hits}")
        _deg_cache[key] = set(open(f"DE/{species}/{hits[0]}").read().split())
    return gene in _deg_cache[key]

info = {}
for s, g in mem.items():
    r = gi[(gi.species == s) & (gi.gene == g)].iloc[0]
    info[s] = {"gene": g,
               "lfc": [r.lfc_6h, r.lfc_24h, r.lfc_3d, r.lfc_10d],
               "deg": [is_deg(s, g, tp) for tp in ("6h", "24h", "3d", "10d")]}

nbr = {}
for s, g in mem.items():
    f = f"data/networks/network_{s}_hrr100.tsv"
    n = pd.read_csv(f, sep="\t")
    sub = n[(n.source == g) | (n.target == g)].copy()
    sub["partner"] = np.where(sub.source == g, sub.target, sub.source)
    nbr[s] = sub

# panel (b): every species pair linked between the six members
pairs = set()
sel = coex[coex.OrthologGroup == OG]
for r in sel.itertuples():
    if r.GeneSpecies1 in mem.values() and r.GeneSpecies2 in mem.values():
        pairs.add(frozenset([r.Species1, r.Species2]))

fig = plt.figure(figsize=(7.2, 5.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.95], width_ratios=[1.25, 1],
                      hspace=0.42, wspace=0.34)

ax = fig.add_subplot(gs[0, :])
x = np.arange(4)
for s in SP:
    y = [np.nan if pd.isna(v) else v for v in info[s]["lfc"]]
    ax.plot(x, y, "-o", color=COL[s], ms=4, lw=1.6, label=LAB[s], zorder=3)
    for k, (v, d) in enumerate(zip(y, info[s]["deg"])):
        if d and not pd.isna(v):
            ax.plot(x[k], v, "o", ms=9, mfc="none", mec=COL[s], mew=1.3, zorder=4)
ax.axhline(0, color=GREY, lw=0.8, zorder=1)
ax.set_xticks(x)
ax.set_xticklabels(["6 h", "24 h", "3 d", "10 d"])
ax.set_ylabel("log$_2$ fold change vs control")
ax.set_title("Repressed in both $\\it{Arabidopsis}$ ecotypes, "
             "induced in aspen, spruce and pine")
ax.margins(0.06)
ax.text(0.015, 0.97, "open ring: differentially expressed (adjusted $P$ ≤ 0.05, |log$_2$FC| ≥ 1)", transform=ax.transAxes,
        ha="left", va="top", fontsize=6, color=GREY)
ax.legend(frameon=False, ncol=6, loc="lower center", bbox_to_anchor=(0.5, -0.30),
          handlelength=1.2, columnspacing=1.1)
ax.text(-0.075, 1.06, "a", transform=ax.transAxes, fontweight="bold", fontsize=10)

axb = fig.add_subplot(gs[1, 0])
ang = np.linspace(90, 450, 7)[:-1] * np.pi / 180
pos = {s: (np.cos(a), np.sin(a)) for s, a in zip(SP, ang)}
for i, s1 in enumerate(SP):
    for s2 in SP[i + 1:]:
        if frozenset([s1, s2]) in pairs:
            axb.plot(*zip(pos[s1], pos[s2]), color="#BBBBBB", lw=0.9, zorder=1)
for s in SP:
    axb.plot(*pos[s], "o", ms=13, color=COL[s], zorder=3)
    axb.text(pos[s][0] * 1.30, pos[s][1] * 1.30, LAB[s],
             ha="center", va="center", fontsize=6)
axb.set_xlim(-1.7, 1.7)
axb.set_ylim(-1.55, 1.6)
axb.axis("off")
axb.set_title(f"All {len(pairs)} pairs are linked:\na complete conserved clique", fontsize=8)
axb.text(-0.02, 1.02, "b", transform=axb.transAxes, fontweight="bold", fontsize=10)

axc = fig.add_subplot(gs[1, 1])
nn = [len(nbr[s]) for s in SP]
ntf = [sum(1 for p in nbr[s].partner if (s, p) in tfset) for s in SP]
yp = np.arange(len(SP))[::-1]
axc.barh(yp, nn, color=[COL[s] for s in SP], height=0.62, zorder=2)
axc.barh(yp, ntf, color="#333333", height=0.26, zorder=3)
for y_, a in zip(yp, nn):
    axc.text(a + 2, y_, f"{a}", va="center", fontsize=6)
axc.set_yticks(yp)
axc.set_yticklabels([LAB[s] for s in SP])
axc.set_xlabel("co-expressed genes")
axc.margins(x=0.16)
axc.set_title("Its local network differs in size;\ndark = transcription factors", fontsize=8)
axc.text(-0.02, 1.02, "c", transform=axc.transAxes, fontweight="bold", fontsize=10)

fig.savefig("Figure_S7.pdf", bbox_inches="tight")
fig.savefig("Figure_S7.png", dpi=300, bbox_inches="tight")
print(f"wrote Figure_S7.pdf/.png for {OG}")
print("  species pairs linked:", len(pairs), "of 15")
print("  co-expressed genes:", dict(zip(SP, nn)))
print("  of which transcription factors:", dict(zip(SP, ntf)))
