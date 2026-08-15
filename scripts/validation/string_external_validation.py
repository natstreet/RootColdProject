#!/usr/bin/env python3
"""
string_external_validation.py

Independent test that the conserved-core CLUSTERS are coherent co-regulated modules.

Tests whether the Arabidopsis members of the conserved core (the clique clusters of Fig. 5)
are more co-expressed WITHIN clusters than BETWEEN clusters in STRING v12 (an external
A. thaliana co-expression compendium). This validates the cluster structure as genuine
co-regulation, complementing the permutation test which establishes the cross-species
conservation itself.

Inputs (download once from STRING; https://string-db.org):
  3702.protein.links.detailed.v12.0.txt.gz   (co-expression channel = column 6)
  3702.protein.aliases.v12.0.txt.gz          (STRING protein -> TAIR locus)
  cliques.xlsx                               (conserved-core cliques: Col0 gene + Cluster)

Method: build the co-expression network at combined score >= THRESHOLD among genes mappable
to TAIR loci; among the conserved-core Arabidopsis genes (with cluster labels), count
co-expression edges that fall WITHIN a cluster; compare to N permutations that randomly
re-assign the same genes to clusters (preserving cluster sizes).

Result (this study, THRESHOLD=400, N=5000): 692 within-cluster edges observed vs ~203
expected; 3.4-fold; empirical P < 0.001.
"""
import numpy as np, re, gzip, pandas as pd

THRESHOLD = 400
N_PERM = 5000
LINKS = "3702.protein.links.detailed.v12.0.txt.gz"
ALIASES = "3702.protein.aliases.v12.0.txt.gz"
CLIQUES = "cliques.xlsx"      # columns: Col0 (AT id), Cluster

def main():
    at2s = {}
    with gzip.open(ALIASES, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and re.fullmatch(r"AT[1-5MC]G\d{5}", p[1].upper()):
                at2s.setdefault(p[1].upper(), p[0])
    mappable = set(at2s.values())

    idx, us, vs = {}, [], []
    with gzip.open(LINKS, "rt") as fh:
        next(fh)
        for line in fh:
            a, b, *c = line.split()
            if int(c[3]) < THRESHOLD or a not in mappable or b not in mappable:
                continue
            us.append(idx.setdefault(a, len(idx))); vs.append(idx.setdefault(b, len(idx)))
    us, vs = np.array(us), np.array(vs)

    cl = pd.read_excel(CLIQUES)[["Col0", "Cluster"]].dropna()
    g2c = {}
    for _, r in cl.iterrows():
        prot = at2s.get(str(r["Col0"]).upper())
        if prot in idx:
            g2c[idx[prot]] = int(r["Cluster"])
    nodeset = set(g2c)
    keep = np.array([(u in nodeset and v in nodeset) for u, v in zip(us, vs)])
    eu, ev = us[keep], vs[keep]
    nodes = np.array(sorted(g2c)); lab = np.array([g2c[n] for n in nodes])
    nix = {n: i for i, n in enumerate(nodes)}
    ei = np.array([nix[u] for u in eu]); ej = np.array([nix[v] for v in ev])
    obs = int(np.sum(np.array([g2c[u] for u in eu]) == np.array([g2c[v] for v in ev])))

    rng = np.random.default_rng(1); null = np.empty(N_PERM, int)
    for i in range(N_PERM):
        p = rng.permutation(lab); null[i] = int(np.sum(p[ei] == p[ej]))
    pval = (np.sum(null >= obs) + 1) / (N_PERM + 1)
    np.savetxt("string_null.csv", null, fmt="%d", header="null", comments="")   # for make_FigS5.py
    open("string_obs.txt", "w").write(str(obs))
    print(f"core genes with cluster + STRING: {len(nodes)}; core-core edges: {len(eu)}")
    print(f"observed within-cluster edges: {obs}")
    print(f"null mean={null.mean():.1f} sd={null.std():.1f} max={null.max()}")
    print(f"fold={obs/max(null.mean(),1e-9):.1f}x ; empirical P={pval:.4g}")

if __name__ == "__main__":
    main()
