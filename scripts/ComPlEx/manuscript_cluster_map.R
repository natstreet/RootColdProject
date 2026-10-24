#!/usr/bin/env Rscript
# manuscript_cluster_map.R — the manuscript number of each analysis cluster of the conserved core: each manuscript
# cluster (as published, scripts/ComPlEx/reference/published_manuscript_clusters.tsv) is matched to the analysis
# cluster (clique_HMcluster_list.rda) sharing most of its orthogroups; the match must be one-to-one.
# Read by make_TableS3.R, make_Fig4.R, manuscript_clusters.R, cluster_profiles.R, string_by_cluster.R.
# Usage (repository root): Rscript scripts/ComPlEx/manuscript_cluster_map.R [--renumber [a1,a2,a3,a4]]   (D081: numbering only)
# Output: data/ComPlEx/cliques/manuscript_cluster_map.tsv (manuscript_cluster, analysis_cluster, shared, analysis_size)
ref <- read.delim("scripts/ComPlEx/reference/published_manuscript_clusters.tsv")
e <- new.env(); load("data/ComPlEx/cliques/clique_HMcluster_list.rda", envir = e); x <- get("clique_cluster_order", e)
cur <- data.frame(orthogroup = sub("_[0-9]+$", "", x$ID), analysis_cluster = x$Cluster)
m <- merge(ref, cur, by = "orthogroup")
tab <- table(m$manuscript_cluster, m$analysis_cluster)
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) && cargs[1] == "--renumber") {
  # D081: the clusters were rebuilt and are not matched to the published four. Numbering only (membership unchanged):
  # an optional second argument gives the analysis cluster for manuscript clusters 1..4 (e.g. 4,2,3,1); without it each
  # cluster keeps its analysis number. shared = NA. The published-cluster overlap is printed for reference.
  print(tab)
  sz <- table(cur$analysis_cluster)
  ana <- if (length(cargs) >= 2) as.integer(strsplit(cargs[2], ",")[[1]]) else as.integer(names(sz))
  if (!setequal(ana, as.integer(names(sz))) || anyDuplicated(ana)) stop("--renumber: give each analysis cluster once")
  out <- data.frame(manuscript_cluster = seq_along(ana), analysis_cluster = ana,
                    shared = NA_integer_, analysis_size = as.integer(sz[as.character(ana)]))
  write.table(out, "data/ComPlEx/cliques/manuscript_cluster_map.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  print(out, row.names = FALSE); quit(save = "no")
}
best <- apply(tab, 1, function(r) as.integer(colnames(tab)[which.max(r)]))
if (anyDuplicated(best) || length(best) != length(unique(cur$analysis_cluster)))
  stop("manuscript and analysis clusters do not match one-to-one:\n", paste(capture.output(print(tab)), collapse = "\n"))
out <- data.frame(manuscript_cluster = as.integer(rownames(tab)), analysis_cluster = best,
                  shared = tab[cbind(rownames(tab), as.character(best))],
                  analysis_size = as.integer(table(cur$analysis_cluster)[as.character(best)]))
write.table(out, "data/ComPlEx/cliques/manuscript_cluster_map.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
