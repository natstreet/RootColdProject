#!/usr/bin/env Rscript
# core_cliques.R — the conserved core: one complete clique per orthogroup, the one with the lowest MeanCliqueSum.
#
# Input:  data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS   every complete clique (one gene per genotype)
# Output: data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS (one row per orthogroup, the columns of the
#         input), or the path given as the first argument.
# Read by: scripts/ComPlEx/clique_clusters.R, make_Fig5.R, make_Fig6.R, make_core_clusters.R, tf_clique_directions.R,
#          complete_cliques_per_orthogroup.R
# Usage:  Rscript scripts/ComPlEx/core_cliques.R [output.RDS]
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args)) args[1] else "data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"

wider <- readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS")
core  <- wider %>% filter(conservedType == "complete") %>%
  group_by(OrthoGroup) %>% slice_min(MeanCliqueSum, n = 1, with_ties = FALSE) %>% ungroup()
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(core, out)
cat("conserved core:", nrow(core), "orthogroups ->", out, "\n")
