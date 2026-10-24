# clique_matrix.R — the matrix behind Figure 4 and the clique clusters: for each conserved-core clique, the per-gene
# z-scored MEDIAN expression of its member in each genotype over the cold time course (Control, 6h, 24h, 3d, 10d).
# Sourced by scripts/ComPlEx/clique_clusters.R and scripts/figures/make_Fig4.R, so both use the same matrix.
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

# D082: one row per core orthogroup built from ALL genes that occur in any of its complete cliques. Per genotype, each
# such gene's profile is computed as above (z-scored median expression over the five conditions); the cell is the
# median of those profiles at each condition. Genes without a median-expression row are left out of the median.
# Rows are named by the orthogroup's representative clique (core_cliques.R), so the cluster table keeps its cliqueID
# key for the scripts that read it. Also returns the gene sets (attribute "genes": OrthoGroup, genotype, gene).
clique_genes_by_genotype <- function(core = "data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS",
                                     wide = "data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS") {
  cl <- readRDS(core); w <- readRDS(wide); w <- w[w$OrthoGroup %in% cl$OrthoGroup, ]
  do.call(rbind, lapply(clique_sp, function(s) unique(data.frame(OrthoGroup = w$OrthoGroup, genotype = s, gene = w[[s]],
                                                                    stringsAsFactors = FALSE))))
}
gene_profiles <- function(s, genes) {          # z-scored median profile per gene, rows named by gene
  e  <- as.data.frame(readr::read_delim(sprintf("data/expression/med_expression/%s_med.txt.gz", s), show_col_types = FALSE))
  cc <- vapply(clique_conds, function(k) grep(paste0("_", k, "$"), colnames(e), value = TRUE)[1], character(1))
  e  <- e[e$Genes %in% genes, ]
  v  <- as.matrix(sapply(e[, cc], as.numeric)); rownames(v) <- e$Genes; colnames(v) <- clique_conds
  z  <- t(apply(v, 1, function(r){ sdv <- sd(r); if (is.na(sdv) || sdv == 0) rep(0, length(r)) else (r - mean(r)) / sdv }))
  colnames(z) <- clique_conds; z
}
family_matrix <- function(core = "data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS",
                          wide = "data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS") {
  cl <- readRDS(core); G <- clique_genes_by_genotype(core, wide)
  M <- matrix(NA_real_, nrow(cl), length(clique_sp) * length(clique_conds), dimnames = list(cl$cliqueID,
              as.vector(t(outer(clique_labs[clique_sp], clique_conds, paste, sep = ".")))))
  for (s in clique_sp) {
    z <- gene_profiles(s, G$gene[G$genotype == s])
    for (i in seq_len(nrow(cl))) {
      g <- intersect(G$gene[G$genotype == s & G$OrthoGroup == cl$OrthoGroup[i]], rownames(z))
      if (length(g)) M[i, paste(clique_labs[[s]], clique_conds, sep = ".")] <- apply(z[g, , drop = FALSE], 2, median)
    }
  }
  ok <- stats::complete.cases(M)
  if (!all(ok)) message("family_matrix: rows left out (a genotype without a profiled gene): ", paste(rownames(M)[!ok], collapse = " "))
  M <- M[ok, , drop = FALSE]; attr(M, "genes") <- G
  M
}

# rows clustered on Euclidean distance with Ward.D2, as in the original analysis (E. van Zalen)
clique_tree <- function(M) hclust(dist(M, method = "euclidean"), method = "ward.D2")
