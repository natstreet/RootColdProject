#!/usr/bin/env Rscript
# export_network_inputs.R — the inputs of build_coexpression_network.R for the six root cold-stress experiments: each
# DESeq2 object's raw counts (<sp>_counts.tsv: Genes + one column per sample) and sample groups (<sp>_groups.tsv:
# sample, group = the last variable of the design, here the treatment).
# Usage (repository root): Rscript scripts/networks/export_network_inputs.R <dds dir> <out dir>
#   e.g. Rscript scripts/networks/export_network_inputs.R data/dds data/networks/inputs
# then, per species:
#   Rscript scripts/networks/build_coexpression_network.R data/networks/inputs/<sp>_counts.tsv \
#       data/networks/inputs/<sp>_groups.tsv data/networks/network_<sp>_hrr100
suppressMessages({ library(data.table); library(DESeq2) })
a <- commandArgs(trailingOnly = TRUE); dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
DDS <- c(col0 = "dds_col0_soil", ost0 = "dds_ost0_soil", aspen = "dds_Pt_new", birch = "dds_Bp_new",
         spruce = "dds_Pa", pine = "dds_Ps")
for (sp in names(DDS)) {
  e <- new.env(); load(file.path(a[1], paste0(DDS[[sp]], ".rda")), envir = e); d <- get(ls(e)[1], e)
  v <- all.vars(design(d)); v <- v[length(v)]
  fwrite(data.table(Genes = rownames(d), as.data.table(counts(d))), file.path(a[2], paste0(sp, "_counts.tsv")), sep = "\t")
  fwrite(data.table(sample = colnames(d), group = as.character(colData(d)[[v]])), file.path(a[2], paste0(sp, "_groups.tsv")), sep = "\t")
  cat(sprintf("%-6s %d genes x %d samples; groups by %s\n", sp, nrow(d), ncol(d), v))
}
