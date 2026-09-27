#!/usr/bin/env Rscript
# make_core_clusters.R — the Col-0 member and heatmap cluster of each conserved-core clique, the input of
# string_external_validation.py (Fig. S6b).
#
# Inputs (run from a directory containing the unpacked figshare deposit `data/`):
#   data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS  116 cliques, one gene per species
#   data/ComPlEx/cliques/clique_HMcluster_list.rda              clique -> cluster (analysis numbering)
# Output: core_clusters.tsv (columns Col0, Cluster)
# Usage:  Rscript scripts/validation/make_core_clusters.R

ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
cg <- as.data.frame(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
hm <- ld(Sys.glob("data/*/cliques/clique_HMcluster_list.rda")[1])
d  <- data.frame(Col0 = cg$col0, Cluster = hm$Cluster[match(cg$cliqueID, hm$ID)])
stopifnot(!anyNA(d$Cluster))
write.table(d, "core_clusters.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("core_clusters.tsv: %d cliques; cluster sizes %s\n", nrow(d), paste(table(d$Cluster), collapse = "/")))
