#!/usr/bin/env Rscript
# aspen_timepoint_degs.R — aspen DEG counts per time point (ledger R045 item 3: the aspen 3 d legend clause). Two definitions,
# both reported: padj < 0.05 alone (as SCI counted), and the study's DEG call (padj <= 0.05 and |log2FC| >= 1, de_calls.R /
# differential_expression.R).
# Inputs (working tree): data/DEGs/padj_pt.RData, data/DEGs/log2FoldChange_pt.RData. Output: results/aspen_timepoint_degs.tsv
# Usage (working tree): Rscript scripts/analysis/aspen_timepoint_degs.R
suppressPackageStartupMessages(library(data.table))
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
pa <- ld("data/DEGs/padj_pt.RData"); lf <- ld("data/DEGs/log2FoldChange_pt.RData"); stopifnot(identical(rownames(pa), rownames(lf)))
res <- rbindlist(lapply(c("6h", "24h", "3d", "10d"), function(t) { p <- pa[[paste0("padj_", t)]]; l <- lf[[paste0("log2FoldChange_", t)]]
  data.table(timepoint = t, padj_below_0.05 = sum(!is.na(p) & p < 0.05), study_degs = sum(!is.na(p) & p <= 0.05 & !is.na(l) & abs(l) >= 1)) }))
fwrite(res, "results/aspen_timepoint_degs.tsv", sep = "\t"); print(res)
