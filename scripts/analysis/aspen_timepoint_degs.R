#!/usr/bin/env Rscript
# aspen_timepoint_degs.R — aspen DEG counts per time point under the study's DEG call: Wald test against no change
# (lfcThreshold = 0), adjusted P <= 0.05 and |log2 fold change| >= 1 (scripts/DE/differential_expression.R, Methods).
# Counted from two sources that must agree: DE/de_results.tsv (column deg) and the per-timepoint lists
# data/DEGs/per_timepoint_DEG_lists/aspen/DEGs_pt_C_vs_<tp>_p0.05_lfc0.txt. Also reported, labelled as such: genes with
# adjusted P < 0.05 against a two-fold null (data/DEGs/padj_pt.RData), a different test that is not the study's call.
# (Ledger A052: an earlier version reported the two-fold-null count as the study's call; 3,141 at 3 d instead of 8,815.)
# Inputs (working tree): DE/de_results.tsv, data/DEGs/per_timepoint_DEG_lists/aspen/, data/DEGs/padj_pt.RData.
# Output: results/aspen_timepoint_degs.tsv
# Usage (working tree): Rscript scripts/analysis/aspen_timepoint_degs.R
suppressPackageStartupMessages(library(data.table))
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
de <- fread("DE/de_results.tsv")[species == "aspen"]; pa <- ld("data/DEGs/padj_pt.RData")
res <- rbindlist(lapply(c("6h", "24h", "3d", "10d"), function(t) {
  from_results <- de[sub("_5C$", "", timepoint) == t & deg %in% c(TRUE, "TRUE", "True"), uniqueN(gene)]
  l <- readLines(sprintf("data/DEGs/per_timepoint_DEG_lists/aspen/DEGs_pt_C_vs_%s_p0.05_lfc0.txt", t)); from_list <- sum(nzchar(l))
  if (from_results != from_list) stop("aspen ", t, ": DE results give ", from_results, " DEGs, the lfc0 list ", from_list)
  p <- pa[[paste0("padj_", t)]]
  data.table(timepoint = t, study_degs = from_results, twofold_null_padj_below_0.05 = sum(!is.na(p) & p < 0.05)) }))
fwrite(res, "results/aspen_timepoint_degs.tsv", sep = "\t"); print(res)
