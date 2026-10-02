#!/usr/bin/env Rscript
# sample_names.R — the column name each sample has in the co-expression input matrices
# (data/expression/<species>_root_cold_stress_expression.txt, e.g. Pt1.01) and its DESeq2 sample name (e.g. PtRWT01).
# The matrices name samples by replicate and time point; the names are recovered by matching each matrix column to
# the DESeq2 sample whose VST profile it equals up to a constant (the smallest standard deviation of the difference),
# and every match must be one-to-one.
# Run from the repository root:  Rscript scripts/expression/sample_names.R <dds dir: the deposit's data/dds>
# Output: data/expression/sample_names.tsv (species, dds_sample, matrix_name)
suppressPackageStartupMessages({ library(DESeq2); library(matrixStats) })
DDS <- commandArgs(trailingOnly = TRUE)[1]
dds <- c(col0 = "dds_col0_soil.rda", ost0 = "dds_ost0_soil.rda", aspen = "dds_Pt_new.rda", birch = "dds_Bp_new.rda",
         spruce = "dds_Pa.rda", pine = "dds_Ps.rda")
out <- list()
for (s in names(dds)) {
  f <- read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", s), check.names = FALSE)
  X <- as.matrix(f[, -1]); rownames(X) <- f[[1]]
  e <- new.env(); load(file.path(DDS, dds[[s]]), envir = e); d <- get(ls(e)[1], e)
  v <- suppressMessages(assay(varianceStabilizingTransformation(d, blind = TRUE)))[rownames(X), ]
  j <- vapply(colnames(X), function(cx) which.min(colSds(v - X[, cx])), integer(1))
  if (anyDuplicated(j)) stop(s, ": a DESeq2 sample matches two matrix columns")
  out[[s]] <- data.frame(species = s, dds_sample = colnames(v)[j], matrix_name = colnames(X))
}
res <- do.call(rbind, out)
write.table(res, "data/expression/sample_names.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote data/expression/sample_names.tsv:", nrow(res), "samples\n")
