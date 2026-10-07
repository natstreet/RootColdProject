#!/usr/bin/env Rscript
# core_cliques.R — the conserved core: one complete clique per orthogroup, the one with the lowest MeanCliqueSum.
# Given the size-rule table (scripts/ComPlEx/size_rule.R; D026, D043), orthogroups with counted_as_conserved FALSE are
# left out of the core and written to a separate file of flagged orthogroups (reported, not hidden).
#
# Input:  data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS   every complete clique (one gene per genotype)
#         [results/orthogroup_size_rule.tsv]                         optional second argument
# Output: data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS (one row per orthogroup, the columns of the
#         input), or the path given as the first argument; with a size-rule table also <output>_flagged.RDS, the
#         same columns for the orthogroups the rule does not count.
# Read by: scripts/ComPlEx/clique_clusters.R, make_Fig5.R, make_Fig6.R, make_core_clusters.R, tf_clique_directions.R,
#          complete_cliques_per_orthogroup.R
# Usage:  Rscript scripts/ComPlEx/core_cliques.R [output.RDS] [size-rule table]
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args)) args[1] else "data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"

wider <- readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS")
core  <- wider %>% filter(conservedType == "complete") %>%
  group_by(OrthoGroup) %>% slice_min(MeanCliqueSum, n = 1, with_ties = FALSE) %>% ungroup()
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
if (length(args) >= 2) {
  rule <- read.delim(args[2])
  if (!setequal(rule$OrthoGroup[rule$has_complete_clique], core$OrthoGroup))
    stop("the size-rule table does not list the same orthogroups with a complete clique")
  keep <- rule$OrthoGroup[rule$counted_as_conserved]
  flagged <- core %>% filter(!OrthoGroup %in% keep)
  core <- core %>% filter(OrthoGroup %in% keep)
  saveRDS(flagged, sub("[.]RDS$", "_flagged.RDS", out))
  cat("not counted by the size rule:", nrow(flagged), "orthogroups:", paste(flagged$OrthoGroup, collapse = ", "), "\n")
}
saveRDS(core, out)
cat("conserved core:", nrow(core), "orthogroups ->", out, "\n")
