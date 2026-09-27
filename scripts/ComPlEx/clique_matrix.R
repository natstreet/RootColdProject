# clique_matrix.R — the matrix behind Figure 5 and the clique clusters: for each conserved-core clique, the per-gene
# z-scored MEDIAN expression of its member in each genotype over the cold time course (Control, 6h, 24h, 3d, 10d).
# Sourced by scripts/ComPlEx/clique_clusters.R and scripts/figures/make_Fig5.R, so both use the same matrix.
#
# Inputs: data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS, data/expression/med_expression/<sp>_med.txt.gz
clique_sp    <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
clique_labs  <- c(col0 = "Col-0", ost0 = "Ost-0", aspen = "Aspen", birch = "Birch", spruce = "Spruce", pine = "Pine")
clique_conds <- c("Control", "6h", "24h", "3d", "10d")

clique_matrix <- function(core = "data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS") {
  cl <- readRDS(core)
  template <- data.frame(cliqueID = cl$cliqueID, stringsAsFactors = FALSE)
  for (s in clique_sp) {
    e  <- as.data.frame(readr::read_delim(sprintf("data/expression/med_expression/%s_med.txt.gz", s), show_col_types = FALSE))
    cc <- vapply(clique_conds, function(k) grep(paste0("_", k, "$"), colnames(e), value = TRUE)[1], character(1))
    ex <- e[, c("Genes", cc)]; colnames(ex) <- c("Genes", clique_conds)
    g  <- data.frame(cliqueID = cl$cliqueID, Genes = cl[[s]], stringsAsFactors = FALSE)
    m  <- dplyr::left_join(g, ex, by = "Genes")
    v  <- as.matrix(sapply(m[, clique_conds], as.numeric))
    z  <- t(apply(v, 1, function(r){ sdv <- sd(r); if (is.na(sdv) || sdv == 0) rep(0, length(r)) else (r - mean(r)) / sdv }))
    colnames(z) <- paste(clique_labs[[s]], clique_conds, sep = ".")
    template <- cbind(template, z)
  }
  M <- as.matrix(template[, -1]); rownames(M) <- template$cliqueID
  M[stats::complete.cases(M), ]
}

# rows clustered on Euclidean distance with Ward.D2, as in the original analysis (E. van Zalen)
clique_tree <- function(M) hclust(dist(M, method = "euclidean"), method = "ward.D2")
