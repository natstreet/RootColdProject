#!/usr/bin/env Rscript
# fig5_caption_candidates.R — proposed GO terms to caption each manuscript cluster of Figure 5, from the Table S4
# output: per cluster, the terms significant (padj < 0.05) in the most species, ties broken by the smallest padj,
# the top N (default 2). The output is a proposal to be checked and edited by the authors before it is copied to
# scripts/figures/Fig5_captions.tsv (cluster, GO.ID), which make_Fig5.R reads.
# Usage (repository root): Rscript scripts/figures/fig5_caption_candidates.R [N] [out.tsv]
a <- commandArgs(trailingOnly = TRUE)
N <- if (length(a) >= 1) as.integer(a[1]) else 2L
out <- if (length(a) >= 2) a[2] else "Fig5_caption_candidates.tsv"
s4 <- read.delim("tables/TableS4_cliqueClusterGO_classic.tsv", check.names = FALSE)
res <- do.call(rbind, lapply(sort(unique(s4[["Cluster (manuscript)"]])), function(cl) {
  x <- s4[s4[["Cluster (manuscript)"]] == cl, ]
  agg <- do.call(rbind, lapply(split(x, x$GO.ID), function(t)
    data.frame(cluster = cl, GO.ID = t$GO.ID[1], Term = t$Term[1], species = length(unique(t$Species)), best_padj = min(t$padj))))
  head(agg[order(-agg$species, agg$best_padj), ], N)
}))
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
