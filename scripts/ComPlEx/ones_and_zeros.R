#!/usr/bin/env Rscript
# ones_and_zeros.R — per orthogroup and species pair, whether ComPlEx tested at least one gene pair (1) or none (0).
#
# The construction follows the ones_and_zeros_all step of Ellen Chapple's cliques workflow
# (CliquesWorkflow/basicFiles.R in the input data and code for "An evo-devo resource for wood: Comparative
# regulomics across dicot and conifer trees", Zenodo, doi:10.5281/zenodo.21025760), applied to this study's six
# genotypes and fifteen pairwise comparisons: every orthogroup of the OrthoFinder table, and for each comparison
# a 1 where the comparison table has any row for the orthogroup (whatever its P value).
#
# Inputs:
#   data/ComPlEx/Orthogroups_20240613.tsv                                    every orthogroup (first column)
#   data/ComPlEx/PairwiseComparisons/comparison-<a>-<b>-abspearsonMR0.03-table.RData   the 15 ComPlEx tables
# Output: data/ComPlEx/ones_and_zeros_all.RData (object ones_and_zeros_all: a data frame of integers, one row per
#         orthogroup named by its ID, one column per comparison "<a>-<b>"), or the path given as the first argument.
# Read by: scripts/ComPlEx/cliques_step1.R, scripts/ComPlEx/cliques_step2.R (both superseded; kept for the record)
# Usage:  Rscript scripts/ComPlEx/ones_and_zeros.R [output.RData]

args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args)) args[1] else "data/ComPlEx/ones_and_zeros_all.RData"

# the comparisons in the column order of the deposited object: angiosperm pairs, then angiosperm-spruce,
# angiosperm-pine and spruce-pine
pairs <- c("col0-ost0", "col0-aspen", "col0-birch", "ost0-aspen", "ost0-birch", "aspen-birch",
           "col0-spruce", "ost0-spruce", "aspen-spruce", "birch-spruce",
           "col0-pine", "ost0-pine", "aspen-pine", "birch-pine", "spruce-pine")

og  <- read.delim("data/ComPlEx/Orthogroups_20240613.tsv", check.names = FALSE, colClasses = "character")[[1]]
if (anyDuplicated(og)) stop("duplicated orthogroup IDs in Orthogroups_20240613.tsv")

ones_and_zeros_all <- data.frame(row.names = og)
for (p in pairs) {
  f <- sprintf("data/ComPlEx/PairwiseComparisons/comparison-%s-abspearsonMR0.03-table.RData", p)
  e <- new.env(); load(f, envir = e)
  tested <- unique(as.character(e$comparison_table$OrthoGroup))
  extra  <- setdiff(tested, og)
  if (length(extra)) stop(f, ": ", length(extra), " orthogroups not in Orthogroups_20240613.tsv")
  ones_and_zeros_all[[p]] <- as.integer(og %in% tested)
}
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
save(ones_and_zeros_all, file = out)
cat("ones_and_zeros_all:", nrow(ones_and_zeros_all), "orthogroups x", ncol(ones_and_zeros_all), "comparisons ->", out, "\n")
