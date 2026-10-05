#!/usr/bin/env Rscript
# core_size_columns.R — add the orthogroup size flag to the deposit's conserved-core table: expressed_genes,
# complete_cliques and p_all15_given_size (scripts/validation/orthogroup_size_flag.R), appended as the last columns.
# Every other column is left byte for byte as it was. Stops if the table already has the columns, or if a core
# orthogroup has no size flag.
#
# Usage (repository root): Rscript scripts/deposit/core_size_columns.R <deposit data/open dir>
open <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(open), dir.exists(open))
p <- file.path(open, "conserved_core.tsv")
x <- read.delim(p, colClasses = "character", na.strings = character(0), check.names = FALSE, quote = "")
f <- read.delim("results/orthogroup_size_flag.tsv", colClasses = "character")
cols <- c("expressed_genes", "complete_cliques", "p_all15_given_size")
if (any(cols %in% names(x))) stop("conserved_core.tsv already has the size columns")
m <- match(x$OrthoGroup, f$OrthoGroup)
if (anyNA(m)) stop("no size flag for ", paste(x$OrthoGroup[is.na(m)], collapse = ", "))
x[cols] <- f[m, cols]
write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("conserved_core.tsv: %s added for %d orthogroups\n", paste(cols, collapse = ", "), nrow(x)))
