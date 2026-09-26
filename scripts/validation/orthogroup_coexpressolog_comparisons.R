#!/usr/bin/env Rscript
# orthogroup_coexpressolog_comparisons.R — in how many of the 15 pairwise species comparisons an orthogroup carries a
# co-expressolog (MaxpVal below the threshold), for the orthogroups named on the command line.
#
# For each orthogroup the output also lists its Arabidopsis members and their symbols, from the deposit's
# orthogroup_annotation.tsv, so that the name used for it in the text can be checked against the data.
#
# Inputs (unpack data.tar.gz from figshare into ./data; the annotation table is in deposit v4 or later):
#   data/ComPlEx/co_expressologs.RData
#   data/open/orthogroup_annotation.tsv
# Output: results/orthogroup_coexpressolog_comparisons.tsv
#   (orthogroup, arabidopsis_loci, arabidopsis_symbols, comparisons_with_coexpressolog, of_comparisons, comparisons)
#
# Usage: Rscript scripts/validation/orthogroup_coexpressolog_comparisons.R [threshold] OG...
#        e.g. Rscript scripts/validation/orthogroup_coexpressolog_comparisons.R 0.1 OG0000037 OG0001369

a <- commandArgs(trailingOnly = TRUE)
thr <- as.numeric(a[1]); ogs <- a[-1]
stopifnot(!is.na(thr), length(ogs) > 0)
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

co <- ld("data/ComPlEx/co_expressologs.RData")
cmp <- apply(cbind(as.character(co$Species1), as.character(co$Species2)), 1, function(v) paste(sort(v), collapse = "-"))
all_cmp <- sort(unique(cmp)); stopifnot(length(all_cmp) == 15)
ann <- read.delim("data/open/orthogroup_annotation.tsv", stringsAsFactors = FALSE, quote = "")

out <- do.call(rbind, lapply(ogs, function(o) {
  k <- sort(unique(cmp[co$OrthologGroup == o & co$MaxpVal < thr]))
  r <- ann[ann$OrthoGroup == o, ]
  data.frame(orthogroup = o,
             arabidopsis_loci = if (nrow(r)) r$Arabidopsis[1] else "",
             arabidopsis_symbols = if (nrow(r)) r$Symbol[1] else "",
             comparisons_with_coexpressolog = length(k), of_comparisons = length(all_cmp),
             comparisons = paste(k, collapse = ";"))
}))
dir.create("results", showWarnings = FALSE)
write.table(out, "results/orthogroup_coexpressolog_comparisons.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out[, 1:5], row.names = FALSE)
