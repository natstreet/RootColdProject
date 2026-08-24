#!/usr/bin/env Rscript
# make_Fig5.R — Figure 5: heatmap of the 116 conserved-core cliques across species.
# Reproduces the published heatmap: for each clique, the per-gene z-scored MEDIAN expression of
# its member in each species over the cold time course (Control, 6h, 24h, 3d, 10d), with rows
# clustered (Euclidean distance, Ward.D2) and cut into four clusters. This follows the original
# analysis (E. van Zalen). The per-cluster GO-term captions, "Cluster n" labels and species names
# on the published figure were added in post-processing and are not part of the plotted object.
#
# Inputs:
#   data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS  116 cliques, one gene per species
#   data/expression/med_expression/<sp>_med.txt.gz             per-timepoint median expression
# Output: Figure5.pdf, Figure5.png (and the four-cluster assignment, printed)
# Usage:  Rscript scripts/figures/make_Fig5.R

suppressPackageStartupMessages({ library(dplyr); library(readr); library(pheatmap); library(grDevices) })

sp    <- c("col0","ost0","aspen","birch","spruce","pine")
labs  <- c(col0="Col-0", ost0="Ost-0", aspen="Aspen", birch="Birch", spruce="Spruce", pine="Pine")
conds <- c("Control","6h","24h","3d","10d")

cl <- readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS")
template <- data.frame(cliqueID = cl$cliqueID, stringsAsFactors = FALSE)

for (s in sp) {
  e  <- as.data.frame(read_delim(sprintf("data/expression/med_expression/%s_med.txt.gz", s), show_col_types = FALSE))
  cc <- vapply(conds, function(k) grep(paste0("_", k, "$"), colnames(e), value = TRUE)[1], character(1))
  ex <- e[, c("Genes", cc)]; colnames(ex) <- c("Genes", conds)
  g  <- data.frame(cliqueID = cl$cliqueID, Genes = cl[[s]], stringsAsFactors = FALSE)
  m  <- left_join(g, ex, by = "Genes")
  v  <- as.matrix(sapply(m[, conds], as.numeric))
  z  <- t(apply(v, 1, function(r){ sdv <- sd(r); if (is.na(sdv) || sdv == 0) rep(0, length(r)) else (r - mean(r)) / sdv }))
  colnames(z) <- paste(labs[[s]], conds, sep = ".")
  template <- cbind(template, z)
}
M <- as.matrix(template[, -1]); rownames(M) <- template$cliqueID
M <- M[stats::complete.cases(M), ]
cat("Cliques in heatmap:", nrow(M), "\n")

pal    <- colorRampPalette(c("#313695", "#f8f8f8", "#A50026"))(100)
breaks <- seq(-3, 3, length.out = 101)                    # white anchored at 0 (matches original)
ph <- pheatmap(M, cluster_rows = TRUE, cluster_cols = FALSE,
               clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
               cutree_rows = 4, gaps_col = c(5, 10, 15, 20, 25),
               color = pal, breaks = breaks, border_color = NA,
               angle_col = 90, fontsize_col = 7, show_rownames = FALSE, show_colnames = TRUE,
               legend_breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
               filename = "Figure5.pdf", width = 9, height = 10)
pheatmap(M, cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_distance_rows = "euclidean", clustering_method = "ward.D2",
         cutree_rows = 4, gaps_col = c(5, 10, 15, 20, 25),
         color = pal, breaks = breaks, border_color = NA,
         angle_col = 90, fontsize_col = 7, show_rownames = FALSE, show_colnames = TRUE,
         legend_breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
         filename = "Figure5.png", width = 9, height = 10)

# report the four-cluster assignment (top-to-bottom heatmap order) for verification
clu <- cutree(ph$tree_row, k = 4)
cat("Cluster sizes (k=4):", paste(table(clu), collapse = ", "), "\n")
saveRDS(data.frame(cliqueID = names(clu), Cluster = unname(clu)), "Fig5_clusters.RDS")
