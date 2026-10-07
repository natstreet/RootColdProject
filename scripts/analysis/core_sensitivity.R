#!/usr/bin/env Rscript
# core_sensitivity.R — core size, manuscript cluster sizes and transcription-factor cliques under different size rules
# (the Methods sensitivity clause). Each argument is label=working directory, a tree in which core_cliques.R,
# clique_clusters.R and manuscript_cluster_map.R have been run (e.g. all orthogroups, the even-odds rule, the rule at the
# lower bootstrap bound via size_rule.R's threshold argument).
# Output (on standard output and in the file named by --out=): label, core_orthogroups, cluster sizes 1-4 in manuscript
# numbering, tf_cliques and their orthogroups.
# Usage (repository root or any tree): Rscript scripts/analysis/core_sensitivity.R --out=results/core_sensitivity.tsv label=dir ...
a <- commandArgs(trailingOnly = TRUE)
out <- sub("^--out=", "", grep("^--out=", a, value = TRUE))
trees <- grep("^--out=", a, value = TRUE, invert = TRUE)
here <- getwd(); tf_src <- normalizePath("scripts/analysis/tf_cliques.R")
rows <- lapply(trees, function(x) {
  lab <- sub("=.*", "", x); dir <- normalizePath(sub("^[^=]*=", "", x))
  setwd(dir); on.exit(setwd(here))
  core <- readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS")
  e <- new.env(); load("data/ComPlEx/cliques/clique_HMcluster_list.rda", envir = e)
  cl <- get("clique_cluster_order", e)
  map <- read.delim("data/ComPlEx/cliques/manuscript_cluster_map.tsv")
  ms <- map$manuscript_cluster[match(cl$Cluster, map$analysis_cluster)]
  sz <- table(factor(ms, levels = sort(map$manuscript_cluster)))
  source(tf_src, local = TRUE); tf <- tf_cliques()
  data.frame(label = lab, core_orthogroups = nrow(core),
             clusters = paste(sprintf("%s:%d", names(sz), as.integer(sz)), collapse = " "),
             tf_cliques = nrow(tf), tf_orthogroups = paste(tf$OrthoGroup, collapse = ","))
})
res <- do.call(rbind, rows)
if (length(out)) write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
