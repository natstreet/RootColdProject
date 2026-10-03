#!/usr/bin/env Rscript
# manuscript_cluster_map.R — the manuscript number of each analysis cluster of the conserved core: each manuscript
# cluster (as published, scripts/ComPlEx/reference/published_manuscript_clusters.tsv) is matched to the analysis
# cluster (clique_HMcluster_list.rda) sharing most of its orthogroups; the match must be one-to-one.
# Read by make_TableS4.R, make_Fig5.R, manuscript_clusters.R, cluster_profiles.R, string_by_cluster.R.
# Usage (repository root): Rscript scripts/ComPlEx/manuscript_cluster_map.R
# Output: data/ComPlEx/cliques/manuscript_cluster_map.tsv (manuscript_cluster, analysis_cluster, shared, analysis_size)
ref <- read.delim("scripts/ComPlEx/reference/published_manuscript_clusters.tsv")
e <- new.env(); load("data/ComPlEx/cliques/clique_HMcluster_list.rda", envir = e); x <- get("clique_cluster_order", e)
cur <- data.frame(orthogroup = sub("_[0-9]+$", "", x$ID), analysis_cluster = x$Cluster)
m <- merge(ref, cur, by = "orthogroup")
tab <- table(m$manuscript_cluster, m$analysis_cluster)
best <- apply(tab, 1, function(r) as.integer(colnames(tab)[which.max(r)]))
if (anyDuplicated(best) || length(best) != length(unique(cur$analysis_cluster)))
  stop("manuscript and analysis clusters do not match one-to-one:\n", paste(capture.output(print(tab)), collapse = "\n"))
out <- data.frame(manuscript_cluster = as.integer(rownames(tab)), analysis_cluster = best,
                  shared = tab[cbind(rownames(tab), as.character(best))],
                  analysis_size = as.integer(table(cur$analysis_cluster)[as.character(best)]))
write.table(out, "data/ComPlEx/cliques/manuscript_cluster_map.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
