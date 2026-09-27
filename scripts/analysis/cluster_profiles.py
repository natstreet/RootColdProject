#!/usr/bin/env python3
"""Cluster-level cold response of the conserved-core cliques (Results, clusters 1 and 2).

Produces the numbers quoted in the Results for manuscript clusters 1 and 2: per cluster,
genotype and timepoint, the median log2 fold change of the clique member genes, the
equivalent fold, and how many genes rise, fall and change at least two-fold. Also lists
cluster 2's members with their Arabidopsis annotation, which is the basis for naming the
member classes in the text.

Inputs (deposit, data/open/): gene_index.tsv, conserved_core.tsv, orthogroup_annotation.tsv.
Cluster numbers: conserved_core.tsv uses the analysis numbering; the manuscript numbering
is analysis {1:4, 2:1, 3:3, 4:2}, as in scripts/GO/make_TableS4.R.
Fold changes are the deposit's log2FoldChange (DESeq2 MLE, not shrunken).

    python3 cluster_profiles.py <deposit data/open dir> <output dir>
"""
import os, sys
import pandas as pd

OPEN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
AN2MS = {1: 4, 2: 1, 3: 3, 4: 2}
SP = ["col0", "ost0", "aspen", "birch", "spruce", "pine"]
TPS = ["6h", "24h", "3d", "10d"]

gi = pd.read_csv(f"{OPEN}/gene_index.tsv", sep="\t")
core = pd.read_csv(f"{OPEN}/conserved_core.tsv", sep="\t")
ann = pd.read_csv(f"{OPEN}/orthogroup_annotation.tsv", sep="\t")
core["ms_cluster"] = core.Cluster.map(AN2MS)
assert core.ms_cluster.value_counts().to_dict() == {1: 50, 2: 30, 3: 20, 4: 16}

# clique member genes: exactly the gene named per genotype in conserved_core.tsv
long = core.melt(id_vars=["OrthoGroup", "ms_cluster"], value_vars=SP,
                 var_name="species", value_name="gene").dropna()
m = long.merge(gi[["species", "gene"] + [f"lfc_{t}" for t in TPS]],
               on=["species", "gene"], how="left", validate="one_to_one")
assert m[[f"lfc_{t}" for t in TPS]].notna().all(axis=1).mean() > 0.99

rows = []
for (c, sp), d in m.groupby(["ms_cluster", "species"]):
    for t in TPS:
        x = d[f"lfc_{t}"].dropna()
        rows.append(dict(ms_cluster=c, species=sp, timepoint=t, n=len(x),
                         median_log2fc=round(x.median(), 3),
                         fold=round(2 ** x.median(), 2) if x.median() >= 0 else round(-2 ** -x.median(), 2),
                         n_up=int((x > 0).sum()), n_down=int((x < 0).sum()),
                         n_abs_log2fc_ge1=int((x.abs() >= 1).sum())))
P = pd.DataFrame(rows)
P["species"] = pd.Categorical(P.species, SP, ordered=True)
P["timepoint"] = pd.Categorical(P.timepoint, TPS, ordered=True)
P = P.sort_values(["ms_cluster", "species", "timepoint"])
P.to_csv(f"{OUT}/cluster_profiles.tsv", sep="\t", index=False)

mem = (core[core.ms_cluster == 2][["OrthoGroup"] + SP]
       .merge(ann[["OrthoGroup", "Arabidopsis", "Symbol", "Name"]], on="OrthoGroup", how="left"))
mem.to_csv(f"{OUT}/cluster2_members.tsv", sep="\t", index=False)
print(f"cluster_profiles.tsv: {len(P)} rows | cluster2_members.tsv: {len(mem)} orthogroups")
