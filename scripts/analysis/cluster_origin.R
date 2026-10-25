#!/usr/bin/env Rscript
# cluster_origin.R — per core orthogroup: its manuscript cluster (D095 numbering, manuscript_cluster_map.tsv), its analysis cluster
# (clique_clusters.R) and its cluster in the published version (scripts/ComPlEx/reference/published_manuscript_clusters.tsv; NA when
# the orthogroup was not in the published core). Replaces the hand-made R059 file of the same name (no producer, F032).
# Output: results/cluster_origin.tsv; prints the manuscript x published cross-table.
# Usage (working tree): Rscript scripts/analysis/cluster_origin.R
suppressPackageStartupMessages(library(data.table))
core <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))[, .(OrthoGroup, cliqueID)]
e <- new.env(); load("data/ComPlEx/cliques/clique_HMcluster_list.rda", envir = e); hm <- as.data.table(get(ls(e)[1], e))
mp <- fread("data/ComPlEx/cliques/manuscript_cluster_map.tsv"); pub <- fread("scripts/ComPlEx/reference/published_manuscript_clusters.tsv")
core[, analysis_cluster := hm$Cluster[match(cliqueID, hm$ID)]]
core[, manuscript_cluster := mp$manuscript_cluster[match(analysis_cluster, mp$analysis_cluster)]]
core[, published_cluster := pub$manuscript_cluster[match(OrthoGroup, pub$orthogroup)]]
stopifnot(!anyNA(core$analysis_cluster), !anyNA(core$manuscript_cluster), !anyDuplicated(core$OrthoGroup))
out <- core[order(manuscript_cluster, OrthoGroup), .(OrthoGroup, manuscript_cluster, analysis_cluster, published_cluster)]
fwrite(out, "results/cluster_origin.tsv", sep = "\t", na = "NA")
print(table(manuscript = out$manuscript_cluster, published = out$published_cluster, useNA = "ifany"))
