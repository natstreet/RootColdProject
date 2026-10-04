#!/usr/bin/env Rscript
# core_neighbourhoods.R — the co-expression neighbourhoods of the conserved-core clique members, as ComPlEx defines
# them, for one species pair. ComPlEx.R keeps only the size of each neighbourhood overlap; this script recovers
# which genes they are.
#
# For the pair, both species' networks are rebuilt exactly as in ComPlEx.R: genes with SD > 0, at least three
# non-zero values and mean expression > 1, that have an orthologue in the other species; |Pearson r|; ranks within
# each row; mutual rank sqrt(R * t(R)); a gene's neighbours are the genes at or above the network's 3% density
# threshold (the gene itself excluded, as in ComPlEx.R). For each core clique, with members g1 (species 1) and
# g2 (species 2), a neighbour of g1 is "overlapping" when it is an orthologue of a neighbour of g2 (and vice versa):
# the genes counted by Species1.neigh.overlap / Species2.neigh.overlap in the comparison table. The script checks
# that its overlap counts equal the table's for every core clique and stops otherwise.
#
# Inputs (repository root): data/expression/<sp>_root_cold_stress_expression.txt, data/ComPlEx/orthologs-<sp1>-<sp2>.RData,
#   data/ComPlEx/PairwiseComparisons/comparison-<sp1>-<sp2>-abspearsonMR0.03-table.RData,
#   data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS
# Output: data/ComPlEx/core_neighbourhoods/<sp1>-<sp2>.tsv.gz
#   (OrthoGroup, species, member, neighbour, overlap: whether the neighbour is in the overlap with the other member)
# Usage:  Rscript scripts/ComPlEx/core_neighbourhoods.R <sp1> <sp2>   (the pair as ordered in ComPlEx.R)
suppressPackageStartupMessages({ library(matrixStats); library(data.table) })
a <- commandArgs(trailingOnly = TRUE); s1 <- a[1]; s2 <- a[2]
density_thr <- 0.03
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
ortho <- ld(sprintf("data/ComPlEx/orthologs-%s-%s.RData", s1, s2))
core <- as.data.frame(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))

expressed <- function(sp) {                                   # ComPlEx.R's two expression filters
  x <- read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", sp), check.names = TRUE)
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
  m <- m[rowSds(m) > 0 & rowCounts(m != 0) >= 3, , drop = FALSE]
  m[rowMeans(m) > 1.0, , drop = FALSE]
}
e1 <- expressed(s1); e2 <- expressed(s2)
ortho <- ortho[ortho$Species1 %in% rownames(e1) & ortho$Species2 %in% rownames(e2), ]
e1 <- e1[rownames(e1) %in% ortho$Species1, , drop = FALSE]
e2 <- e2[rownames(e2) %in% ortho$Species2, , drop = FALSE]

# neighbourhoods of the given genes in one species' mutual-rank network
neighbourhoods <- function(e, genes) {
  net <- abs(cor(t(e), method = "pearson")); dimnames(net) <- list(rownames(e), rownames(e))
  R <- t(apply(net, 1, rank)); rm(net); gc()
  net <- sqrt(R * t(R)); rm(R); gc()
  diag(net) <- 0                                              # as ComPlEx.R: a gene is not its own neighbour
  v <- sort(net[upper.tri(net, diag = FALSE)], decreasing = TRUE)
  thr <- v[round(density_thr * length(v))]; rm(v); gc()
  genes <- intersect(genes, rownames(net))
  out <- lapply(genes, function(g) { r <- net[g, ]; names(r)[r >= thr] }); names(out) <- genes
  rm(net); gc()
  list(neigh = out, thr = thr, n = nrow(e))
}
cl <- core[!is.na(core[[s1]]) & !is.na(core[[s2]]), c("OrthoGroup", s1, s2)]
n1 <- neighbourhoods(e1, cl[[s1]]); n2 <- neighbourhoods(e2, cl[[s2]])
cat(sprintf("%s-%s: %d and %d genes; thresholds %.1f and %.1f\n", s1, s2, n1$n, n2$n, n1$thr, n2$thr))

tab <- as.data.frame(ld(sprintf("data/ComPlEx/PairwiseComparisons/comparison-%s-%s-abspearsonMR0.03-table.RData", s1, s2)))
rows <- list(); bad <- 0
for (i in seq_len(nrow(cl))) {
  g1 <- cl[[s1]][i]; g2 <- cl[[s2]][i]
  if (is.null(n1$neigh[[g1]]) || is.null(n2$neigh[[g2]])) stop(cl$OrthoGroup[i], ": member not in the pair's network")
  nb1 <- n1$neigh[[g1]]; nb2 <- n2$neigh[[g2]]
  ov1 <- unique(intersect(nb1, ortho$Species1[ortho$Species2 %in% nb2]))
  ov2 <- unique(intersect(nb2, ortho$Species2[ortho$Species1 %in% nb1]))
  t <- tab[tab$Species1 == g1 & tab$Species2 == g2, ]
  if (nrow(t) != 1 || t$Species1.neigh.overlap != length(ov1) || t$Species2.neigh.overlap != length(ov2)) {
    bad <- bad + 1; cat("MISMATCH", cl$OrthoGroup[i], g1, g2, length(ov1), length(ov2), "table:",
                        t$Species1.neigh.overlap, t$Species2.neigh.overlap, "\n") }
  rows[[length(rows) + 1]] <- data.table(OrthoGroup = cl$OrthoGroup[i], species = s1, member = g1, neighbour = nb1, overlap = nb1 %in% ov1)
  rows[[length(rows) + 1]] <- data.table(OrthoGroup = cl$OrthoGroup[i], species = s2, member = g2, neighbour = nb2, overlap = nb2 %in% ov2)
}
if (bad) stop(bad, " of ", nrow(cl), " core cliques do not reproduce the comparison table's overlap counts")
dir.create("data/ComPlEx/core_neighbourhoods", showWarnings = FALSE, recursive = TRUE)
f <- sprintf("data/ComPlEx/core_neighbourhoods/%s-%s.tsv.gz", s1, s2)
fwrite(rbindlist(rows), f, sep = "\t")
cat(sprintf("%d core cliques; overlap counts equal the comparison table for all; wrote %s\n", nrow(cl), f))
