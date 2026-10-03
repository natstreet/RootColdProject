#!/usr/bin/env Rscript
# deposit_de_objects.R — the deposited DE objects the figure and table producers read, written from
# differential_expression.R's per-gene results (DE/de_results.tsv):
#   DEGs/per_timepoint_DEG_lists/<species>/DEGs_<species>_C_vs_<tp>_p0.05_lfc0.txt   the study's DEG call per time point
#   DEGs/DEGs_<species>.RData            DEGs_<species>: the union over time points
#   DEGs/log2FoldChange_<code>.RData     log2FoldChange_<code>: one column per time point (log2FoldChange_<tp>) and
#                                        TopL2FC, the value with the largest absolute fold change over time points
#                                        (a list column of named values, the time point as the name, as deposited)
# Usage (repository root): Rscript scripts/DE/deposit_de_objects.R <DE/de_results.tsv> <out data dir> <species ...>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
res <- fread(a[1]); OUT <- a[2]; SPECIES <- a[-(1:2)]
CODE <- c(col0 = "col0", ost0 = "ost0", aspen = "pt", birch = "bp", spruce = "pa", pine = "ps")
TP <- c("6h", "24h", "3d", "10d")
for (sp in SPECIES) {
  r <- res[species == sp]
  tps <- unique(r$timepoint); tp_order <- tps[order(match(sub("_5C$", "", tps), TP))]
  dir.create(file.path(OUT, "DEGs/per_timepoint_DEG_lists", sp), recursive = TRUE, showWarnings = FALSE)
  for (tp in tp_order) {
    g <- r[timepoint == tp & deg == TRUE, sort(gene)]
    writeLines(g, file.path(OUT, "DEGs/per_timepoint_DEG_lists", sp, sprintf("DEGs_%s_C_vs_%s_p0.05_lfc0.txt", CODE[[sp]], sub("_5C$", "", tp))))
  }
  degs <- sort(unique(r[deg == TRUE, gene]))
  assign(paste0("DEGs_", sp), degs); save(list = paste0("DEGs_", sp), file = file.path(OUT, sprintf("DEGs/DEGs_%s.RData", sp)))
  w <- dcast(r, gene ~ timepoint, value.var = "log2FoldChange")
  m <- as.matrix(w[, ..tp_order]); rownames(m) <- w$gene
  lf <- as.data.frame(m); names(lf) <- paste0("log2FoldChange_", sub("_5C$", "", tp_order))
  j <- max.col(abs(replace(m, is.na(m), 0)), ties.method = "first")      # stored as the deposit stores it: a list
  lf$TopL2FC <- setNames(lapply(seq_len(nrow(m)), function(i) setNames(m[i, j[i]], names(lf)[j[i]])), rownames(m))
  nm <- paste0("log2FoldChange_", CODE[[sp]]); assign(nm, lf); save(list = nm, file = file.path(OUT, sprintf("DEGs/%s.RData", nm)))
  cat(sprintf("%-6s %d DEGs over %s -> %s\n", sp, length(degs), paste(tp_order, collapse = ","), OUT))
}
