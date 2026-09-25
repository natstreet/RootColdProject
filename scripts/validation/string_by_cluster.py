#!/usr/bin/env python3
"""
string_by_cluster.py

Breaks the Fig. S5b STRING test down by clique cluster and repeats it without manuscript cluster 2 (the
translation / ribosome-biogenesis cluster), to ask whether cluster coherence in STRING holds beyond that cluster.

Same data and settings as string_external_validation.py: STRING v12 A. thaliana co-expression channel, subscore
>= 400, each gene pair counted once, core genes that have at least one such link, 5,000 permutations of cluster
labels (cluster sizes kept), seed 1.

  1. per cluster, from the full permutation: observed and null-mean within-cluster pairs, empirical P
  2. manuscript cluster 2 excluded: the same test on the remaining genes
  3. the within-cluster pairs of every cluster, with Arabidopsis gene symbols

Run from a directory holding data/ (the deposit plus this repository's data/), core_clusters.tsv (from
make_core_clusters.R) and the two STRING files; manuscript cluster numbers are read from `MAP` in
scripts/GO/make_TableS3.R.
Outputs: string_by_cluster.tsv (tests 1 and 2), string_within_cluster_pairs.tsv (3)
"""
import csv, gzip, re
import numpy as np, pandas as pd

THRESHOLD, N_PERM, SEED = 400, 5000, 1
LINKS = "3702.protein.links.detailed.v12.0.txt.gz"
ALIASES = "3702.protein.aliases.v12.0.txt.gz"


def manuscript_map():
    """analysis cluster -> manuscript cluster, from MAP <- c(`ms` = analysis, ...) in make_TableS3.R."""
    m = re.search(r"MAP\s*<-\s*c\(([^)]*)\)", open("scripts/GO/make_TableS3.R").read())
    ms2an = {int(k): int(v) for k, v in re.findall(r"`?(\d+)`?\s*=\s*(\d+)", m.group(1))}
    return {v: k for k, v in ms2an.items()}


def perm_test(g2c, pairs, rng):
    """Observed within-cluster pairs (total and per cluster) and their label-permutation nulls."""
    nodes = sorted(g2c); lab = np.array([g2c[n] for n in nodes]); ix = {n: i for i, n in enumerate(nodes)}
    ei = np.array([ix[a] for a, b in pairs], int); ej = np.array([ix[b] for a, b in pairs], int)
    ks = sorted(set(lab))
    def counts(l):
        same = l[ei] == l[ej]
        return int(same.sum()), {k: int(np.sum(same & (l[ei] == k))) for k in ks}
    obs, obs_k = counts(lab)
    null = np.empty(N_PERM, int); null_k = {k: np.empty(N_PERM, int) for k in ks}
    for i in range(N_PERM):
        t, tk = counts(rng.permutation(lab)); null[i] = t
        for k in ks: null_k[k][i] = tk[k]
    size = {k: int(np.sum(lab == k)) for k in ks}
    return obs, null, obs_k, null_k, size


def main():
    at2s = {}
    with gzip.open(ALIASES, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and re.fullmatch(r"AT[1-5MC]G\d{5}", p[1].upper()):
                at2s.setdefault(p[1].upper(), p[0])
    s2at = {v: k for k, v in at2s.items()}
    mappable = set(at2s.values())
    linked, edges = set(), set()
    with gzip.open(LINKS, "rt") as fh:
        next(fh)
        for line in fh:
            a, b, *c = line.split()
            if a >= b or int(c[3]) < THRESHOLD or a not in mappable or b not in mappable:
                continue
            linked.update((a, b)); edges.add((a, b))

    an2ms = manuscript_map()
    cl = pd.read_csv("core_clusters.tsv", sep="\t")
    g2c = {}
    for g, c in zip(cl["Col0"], cl["Cluster"]):
        s = at2s.get(str(g).upper())
        if s in linked:
            g2c[s] = an2ms[int(c)]                       # manuscript cluster number
    core_pairs = [(a, b) for a, b in edges if a in g2c and b in g2c]

    rows = []
    obs, null, obs_k, null_k, size = perm_test(g2c, core_pairs, np.random.default_rng(SEED))
    p = lambda o, n: (np.sum(n >= o) + 1) / (N_PERM + 1)
    rows.append(dict(test="all clusters", cluster="all", genes=len(g2c), core_pairs=len(core_pairs), observed=obs,
                     null_mean=round(null.mean(), 2), null_sd=round(null.std(), 2), fold=round(obs / null.mean(), 2),
                     empirical_P=p(obs, null)))
    for k in sorted(obs_k):
        nm = null_k[k].mean()
        rows.append(dict(test="all clusters, per cluster", cluster=k, genes=size[k], core_pairs="", observed=obs_k[k],
                         null_mean=round(nm, 2), null_sd=round(null_k[k].std(), 2),
                         fold=round(obs_k[k] / nm, 2) if nm > 0 else "", empirical_P=p(obs_k[k], null_k[k])))

    g2c_x = {g: c for g, c in g2c.items() if c != 2}
    pairs_x = [(a, b) for a, b in core_pairs if a in g2c_x and b in g2c_x]
    obs, null, obs_k, null_k, size = perm_test(g2c_x, pairs_x, np.random.default_rng(SEED))
    rows.append(dict(test="cluster 2 excluded", cluster="1,3,4", genes=len(g2c_x), core_pairs=len(pairs_x), observed=obs,
                     null_mean=round(null.mean(), 2), null_sd=round(null.std(), 2),
                     fold=round(obs / null.mean(), 2) if null.mean() > 0 else "", empirical_P=p(obs, null)))
    for k in sorted(obs_k):
        rows.append(dict(test="cluster 2 excluded, genes per cluster", cluster=k, genes=size[k], core_pairs="",
                         observed=obs_k[k], null_mean=round(null_k[k].mean(), 2), null_sd=round(null_k[k].std(), 2),
                         fold="", empirical_P=p(obs_k[k], null_k[k])))
    pd.DataFrame(rows).to_csv("string_by_cluster.tsv", sep="\t", index=False)

    sym = {}
    for r in csv.DictReader(open("data/annotation/gene_aliases_20140331.txt"), delimiter="\t"):
        sym.setdefault(r["locus_name"].upper(), r["symbol"])
    with open("string_within_cluster_pairs.tsv", "w") as fh:
        fh.write("cluster\tgene1\tsymbol1\tgene2\tsymbol2\n")
        for a, b in sorted(core_pairs, key=lambda e: (g2c[e[0]], s2at[e[0]], s2at[e[1]])):
            if g2c[a] == g2c[b]:
                fh.write(f"{g2c[a]}\t{s2at[a]}\t{sym.get(s2at[a], '')}\t{s2at[b]}\t{sym.get(s2at[b], '')}\n")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
