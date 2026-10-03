#!/usr/bin/env Rscript
# clique_clusters.R — the four clusters of the conserved-core cliques (analysis numbering): the Figure 5 matrix
# (clique_matrix.R) clustered with Ward.D2 on Euclidean distance and cut into four, rows in dendrogram order.
# The manuscript numbers these clusters differently: data/ComPlEx/cliques/manuscript_cluster_map.tsv (manuscript_cluster_map.R).
#
# Inputs: data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS (from core_cliques.R),
#         data/expression/med_expression/<sp>_med.txt.gz
# Output: data/ComPlEx/cliques/clique_HMcluster_list.rda (object clique_cluster_order: ID = cliqueID, Cluster), or
#         the path given as the first argument.
# Read by: make_Fig5.R, make_TableS4.R
# Usage:  Rscript scripts/ComPlEx/clique_clusters.R [output.rda]
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args)) args[1] else "data/ComPlEx/cliques/clique_HMcluster_list.rda"
source("scripts/ComPlEx/clique_matrix.R")

M   <- clique_matrix()
tr  <- clique_tree(M)
cut <- cutree(tr, k = 4)
clique_cluster_order <- data.frame(ID = rownames(M)[tr$order], Cluster = unname(cut[tr$order]))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
save(clique_cluster_order, file = out)
cat("clique clusters:", nrow(clique_cluster_order), "cliques; sizes", paste(table(clique_cluster_order$Cluster), collapse = "/"),
    "->", out, "\n")
