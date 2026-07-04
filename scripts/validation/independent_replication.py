#!/usr/bin/env python3
"""
independent_replication.py

Does the conserved core reproduce in an INDEPENDENT expression dataset?

RATIONALE
---------
The strongest fast validation is to show that the conserved co-expressologs
re-emerge in data you did not generate. NB: the spruce data here is the same
as the published Vergara spruce data, so spruce CANNOT serve as the independent
set. Genuinely independent options:
  - Wang et al. 2023 (Sci Data) 11-species angiosperm cold phylotranscriptome
    (ENA/their repository) - overlaps Arabidopsis + adds independent angiosperms.
  - ATTED-II / public Arabidopsis co-expression compendia for the Arabidopsis arm.
  - any independent Populus / Betula expression compendium for the broadleaf arm.

WHAT IT DOES
------------
1. Builds (or loads) a co-expression network from the INDEPENDENT dataset for a
   species/arm, using the same density/mutual-rank rule as the main paper.
2. For each conserved-core orthogroup, tests whether its co-expression
   neighbourhood in the independent network overlaps the neighbourhood seen in
   the original network more than expected (hypergeometric), i.e. whether the
   co-expressolog relationship replicates.
3. Reports the replication rate of the conserved core and its significance
   against a degree-matched random set of orthogroups (so replication isn't just
   a hub effect).

INPUTS (edit CONFIG)
--------------------
independent_expr.tsv : genes x samples expression matrix (VST/log) for the
                       independent dataset, gene IDs mappable to orthogroups.
original_network.tsv : gene_a <tab> gene_b for the matching species in the paper.
orthogroups.tsv      : orthogroup <tab> species <tab> gene  (long format).
core_orthogroups.txt : one conserved-core orthogroup id per line.
SPECIES              : the species/arm being replicated (must match labels).

DEPS: numpy, pandas, scipy, networkx
Run:  python independent_replication.py
"""
from __future__ import annotations
import numpy as np, pandas as pd, networkx as nx
from collections import defaultdict
from scipy.stats import hypergeom

# ------------- CONFIG -------------
INDEP_EXPR   = "independent_expr.tsv"
ORIG_NET     = "original_network.tsv"
ORTHO_FILE   = "orthogroups.tsv"
CORE_FILE    = "core_orthogroups.txt"
SPECIES      = "arabidopsis"
EDGE_DENSITY = 0.03      # top fraction of edges retained (match the paper)
MIN_CORR     = None      # or set an absolute Pearson cutoff instead of density
N_RANDOM     = 1000      # degree-matched random orthogroup draws for the null
SEED         = 1
OUT          = "replication_result.csv"
# ----------------------------------
rng = np.random.default_rng(SEED)

def build_network(expr_path, density):
    X = pd.read_csv(expr_path, sep="\t", index_col=0)
    genes = X.index.to_numpy()
    C = np.corrcoef(X.to_numpy())
    np.fill_diagonal(C, 0.0)
    iu = np.triu_indices_from(C, k=1)
    w = C[iu]
    if MIN_CORR is not None:
        keep = w >= MIN_CORR
    else:
        thr = np.quantile(np.abs(w), 1 - density)
        keep = np.abs(w) >= thr
    g = nx.Graph()
    g.add_nodes_from(genes)
    for a, b in zip(*(idx[keep] for idx in iu)):
        g.add_edge(genes[a], genes[b])
    return g

def load_net(path):
    g = nx.Graph()
    for line in open(path):
        if line.strip() and not line.startswith("#"):
            a, b = line.split("\t")[:2]
            g.add_edge(a.strip(), b.strip())
    return g

def main():
    print(f"Replicating conserved core in independent data for: {SPECIES}")
    df = pd.read_csv(ORTHO_FILE, sep="\t", names=["og","species","gene"], header=0)
    df = df[df.species == SPECIES]
    gene2og = dict(zip(df.gene, df.og))
    og2genes = defaultdict(set)
    for r in df.itertuples(): og2genes[r.og].add(r.gene)
    core = [l.strip() for l in open(CORE_FILE) if l.strip()]

    g_orig  = load_net(ORIG_NET)
    g_indep = build_network(INDEP_EXPR, EDGE_DENSITY)
    print(f"  original net: {g_orig.number_of_edges()} edges; "
          f"independent net: {g_indep.number_of_edges()} edges")

    def nbr_ogs(g, node):
        return {gene2og.get(n) for n in g.neighbors(node)} - {None} if node in g else set()

    universe = {gene2og.get(n) for n in g_indep.nodes()} - {None}
    M = len(universe)

    def replicates(og):
        genes = [x for x in og2genes[og] if x in g_orig and x in g_indep]
        for gene in genes:
            A = nbr_ogs(g_orig, gene) & universe
            B = nbr_ogs(g_indep, gene) & universe
            k = len(A & B)
            if k and len(A) and len(B):
                p = hypergeom.sf(k - 1, M, len(A), len(B))
                if p < 0.05:
                    return True
        return False

    core_rep = [og for og in core if replicates(og)]
    rate = len(core_rep) / max(1, len(core))
    print(f"  conserved-core replication: {len(core_rep)}/{len(core)} = {rate:.1%}")

    # degree-matched null: random orthogroups with similar node degree
    deg = dict(g_orig.degree())
    core_deg = [np.mean([deg.get(x,0) for x in og2genes[og] if x in deg]) for og in core]
    all_ogs = [og for og in universe if any(x in g_orig for x in og2genes[og])]
    null_rates = []
    for _ in range(N_RANDOM):
        samp = rng.choice(all_ogs, size=len(core), replace=False)
        null_rates.append(np.mean([replicates(og) for og in samp]))
    null_rates = np.array(null_rates)
    p_emp = (np.sum(null_rates >= rate) + 1) / (N_RANDOM + 1)
    print(f"  null replication rate: {null_rates.mean():.1%} +/- {null_rates.std():.1%}")
    print(f"  empirical p (core > random): {p_emp:.4g}")

    pd.DataFrame({"orthogroup": core, "replicates": [og in set(core_rep) for og in core]}
                 ).to_csv(OUT, index=False)
    print(f"  wrote {OUT}")

if __name__ == "__main__":
    main()
