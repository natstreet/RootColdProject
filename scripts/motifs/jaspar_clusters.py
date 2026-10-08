#!/usr/bin/env python3
"""jaspar_clusters.py — motif-to-cluster table for the promoter motif track (DECISIONS.md D062 item 4, D063; A035).

A035 scores motifs at the level of JASPAR 2024 plant matrix clusters, so a family's near-identical matrices count once.
JASPAR publishes the plant CORE clustering only as an interactive tree page; this script reads cluster membership from
that page (one `d3.select("#cluster_NN").append("svg")` block per cluster, followed by its `var jsontree`, whose leaves
carry `"label":"JASPAR_plants_CORE_<matrix>_n<k>"` and `"name":"<TF>"`) and checks it against the motif file.

Inputs (downloaded 2026-10-08, see external/jaspar2024/SOURCE.tsv):
  JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt
  JASPAR_2024_matrix_clustering_plants_archive.zip (plants_CORE_D3_dynamic_tree.html)
Output: TSV with matrix_id, tf_name, cluster. Stops if any matrix of the motif file has no cluster, or vice versa.
Usage: python3 jaspar_clusters.py <meme file> <cluster archive zip> <out tsv>
"""
import re, sys, zipfile

meme, arch, out = sys.argv[1:4]
with zipfile.ZipFile(arch) as z:
    html = z.read("plants_CORE_D3_dynamic_tree.html").decode()
starts = [(m.start(), m.group(1)) for m in re.finditer(r'd3\.select\("#(cluster_\d+)"\)\.append\("svg"\)', html)]
rows = []
for i, (pos, cl) in enumerate(starts):
    end = starts[i + 1][0] if i + 1 < len(starts) else len(html)
    block = html[pos:end]
    for lab, name in re.findall(r'"label"\s*:\s*"JASPAR_plants_CORE_(MA\d+\.\d+)_n\d+".*?"name"\s*:\s*"([^"]*)"', block, re.S):
        rows.append((lab, name, cl))
mat = re.findall(r"^MOTIF\s+(\S+)\s+(\S+)", open(meme).read(), re.M)
ids_meme = {m for m, _ in mat}
ids_cl = [r[0] for r in rows]
if len(ids_cl) != len(set(ids_cl)):
    sys.exit("a matrix appears in more than one cluster")
miss = ids_meme - set(ids_cl); extra = set(ids_cl) - ids_meme
if miss or extra:
    sys.exit(f"cluster table and motif file disagree: {len(miss)} matrices without cluster, {len(extra)} clustered matrices not in file")
with open(out, "w") as fh:
    fh.write("matrix_id\ttf_name\tcluster\n")
    for r in sorted(rows):
        fh.write("\t".join(r) + "\n")
print(f"{len(rows)} matrices in {len(set(r[2] for r in rows))} clusters -> {out}")
