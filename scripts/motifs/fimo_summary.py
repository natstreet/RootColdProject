#!/usr/bin/env python3
"""fimo_summary.py — gene x matrix hit counts from FIMO --text output (promoter motif track, D062 item 4, D063).
Hits are FIMO's own (p below the --thresh given to FIMO, both strands). Writes gene, matrix_id, n_hits (gzipped TSV).
Usage: fimo_summary.py <fimo.tsv(.gz)> <out.tsv.gz>
"""
import collections, gzip, sys
src, out = sys.argv[1:3]
op = gzip.open if src.endswith(".gz") else open
n = collections.Counter()
with op(src, "rt") as fh:
    head = fh.readline().rstrip("\n").split("\t")
    im, iseq = head.index("motif_id"), head.index("sequence_name")
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        n[(f[iseq], f[im])] += 1
with gzip.open(out, "wt") as fo:
    fo.write("gene\tmatrix_id\tn_hits\n")
    for (g, m), k in sorted(n.items()):
        fo.write(f"{g}\t{m}\t{k}\n")
print(f"{len(n)} gene-matrix pairs with hits -> {out}")
