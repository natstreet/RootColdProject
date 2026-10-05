#!/usr/bin/env Rscript
# manuscript_clusters.R — add the manuscript numbering of the conserved-core clique clusters to the deposit tables
# that carry the analysis numbering.
#
# The analysis (clique_HMcluster_list.rda) numbers the four clusters differently from the manuscript. The mapping is
# MAP in scripts/GO/make_TableS4.R, which is also what make_Fig5.R reads; it is read from there, not restated.
#
# Edits in place, leaving every other column byte for byte as it was:
#   <open>/conserved_core.tsv   adds ManuscriptCluster after Cluster (analysis numbering)
#   <open>/gene_index.tsv       adds core_cluster_manuscript after core_cluster (analysis numbering; empty for
#                               genes outside the conserved core); skipped when the directory has no gene_index.tsv
# Stops if a table already has the column, or if a cluster number has no manuscript counterpart.
#
# Usage: Rscript scripts/deposit/manuscript_clusters.R <deposit data/open dir>

open <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(open), dir.exists(open))

mm <- read.delim("data/ComPlEx/cliques/manuscript_cluster_map.tsv")               # scripts/ComPlEx/manuscript_cluster_map.R
MAP <- setNames(mm$analysis_cluster, mm$manuscript_cluster)                                # manuscript -> analysis
stopifnot(length(MAP) == 4, setequal(MAP, 1:4))
to_ms <- function(a) {
  a <- trimws(a); out <- rep("", length(a)); k <- a != ""
  m <- match(as.integer(as.numeric(a[k])), MAP); stopifnot(!anyNA(m))
  out[k] <- names(MAP)[m]; out
}

add_after <- function(file, after, new) {
  p <- file.path(open, file)
  x <- read.delim(p, colClasses = "character", na.strings = character(0), check.names = FALSE, quote = "")
  if (new %in% names(x)) stop(file, " already has ", new)
  i <- match(after, names(x)); stopifnot(!is.na(i))
  x[[new]] <- to_ms(x[[after]])
  x <- x[, append(names(x)[-ncol(x)], new, after = i)]
  write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat(sprintf("%s: %s added (%s)\n", file, new,
              paste(sprintf("%s->%s", names(table(x[[after]][x[[after]] != ""])),
                            to_ms(names(table(x[[after]][x[[after]] != ""])))), collapse = ", ")))
}

add_after("conserved_core.tsv", "Cluster", "ManuscriptCluster")
if (file.exists(file.path(open, "gene_index.tsv"))) add_after("gene_index.tsv", "core_cluster", "core_cluster_manuscript") else
  cat("gene_index.tsv: not in", open, "- skipped\n")
