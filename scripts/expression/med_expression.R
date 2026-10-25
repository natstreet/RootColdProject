#!/usr/bin/env Rscript
# med_expression.R — per-timepoint median expression per gene (ledger R068 F02): for each genotype, the median over that time point's
# samples of the deposited VST matrix (data/expression/<genotype>_root_cold_stress_expression.txt), samples assigned to time points by
# data/open/samples.tsv. Read by clique_matrix.R (Fig. 4, the D082 clusters, paralogue agreement), make_Fig1.R and make_FigS1.R.
# Replaces carried-over tables (dated 14 Oct, no producer) that did not equal these medians for Col-0, Ost-0, birch and Scots pine.
# Output: data/expression/med_expression/<genotype>_med.txt.gz (Genes, <Prefix>_Control, _6h, _24h, _3d, _10d)
# Usage (working tree): Rscript scripts/expression/med_expression.R
suppressPackageStartupMessages(library(data.table))
smp <- fread("data/open/samples.tsv")
tp_lab <- c(control = "Control", "6h" = "6h", "24h" = "24h", "3d" = "3d", "10d" = "10d")
pre <- c(col0 = "Col", ost0 = "Ost", aspen = "Asp", birch = "Bir", spruce = "Spr", pine = "Pin")
dir.create("data/expression/med_expression", showWarnings = FALSE, recursive = TRUE)
for (s in names(pre)) {
  e <- fread(sprintf("data/expression/%s_root_cold_stress_expression.txt", s)); m <- as.matrix(e[, -1])
  tp <- smp[species == s][match(colnames(m), sample_id), timepoint]
  if (anyNA(tp) || !all(names(tp_lab) %in% tp)) stop(s, ": samples without a time point, or a time point without samples")
  out <- data.table(Genes = e[[1]])
  for (t in names(tp_lab)) out[[paste0(pre[[s]], "_", tp_lab[[t]])]] <- apply(m[, tp == t, drop = FALSE], 1, median)
  f <- gzfile(sprintf("data/expression/med_expression/%s_med.txt.gz", s), "w"); write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE); close(f)
  cat(sprintf("%s: %d genes; samples per time point %s\n", s, nrow(out), paste(table(factor(tp, names(tp_lab))), collapse = "/")))
}
