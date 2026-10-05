#!/usr/bin/env Rscript
# deposit_de_table.R — the deposit's long-format differential-expression table (differential_expression.tsv:
# species, gene, timepoint, log2FoldChange, padj) from the DESeq2 objects: DESeq() on each object, each cold timepoint
# against control, log2FoldChange the maximum-likelihood estimate and padj the Benjamini-Hochberg adjusted P of the
# Wald test against a two-fold change (lfcThreshold = 1), as in the deposit. The Methods DEG call (lfcThreshold = 0)
# is a separate table (comparative_loader/de_methods_call.R).
# Usage: Rscript deposit_de_table.R <dds dir> <out tsv> [species ...]   (default: all six)
suppressMessages({ library(DESeq2); library(data.table) })
a <- commandArgs(trailingOnly = TRUE)
SP <- c(col0 = "dds_col0_soil", ost0 = "dds_ost0_soil", aspen = "dds_Pt_new", birch = "dds_Bp_new", spruce = "dds_Pa", pine = "dds_Ps")
sp <- if (length(a) > 2) a[-(1:2)] else names(SP)
TPS <- c("6h", "24h", "3d", "10d")
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
out <- rbindlist(lapply(sp, function(s) {
  dds <- DESeq(ld(file.path(a[1], paste0(SP[[s]], ".rda"))), quiet = TRUE)
  dv <- all.vars(design(dds))[1]; levs <- setdiff(levels(colData(dds)[[dv]]), "control")
  rbindlist(lapply(levs, function(lv) {
    r <- results(dds, contrast = c(dv, lv, "control"), lfcThreshold = 1)
    data.table(species = s, gene = rownames(r), timepoint = sub("_5C$", "", lv), log2FoldChange = r$log2FoldChange, padj = r$padj)
  }))
}))
out[, timepoint := factor(timepoint, TPS)]; setorder(out, species, gene, timepoint)
fwrite(out, a[2], sep = "\t", na = "NA")
cat(sprintf("%d rows for %s -> %s\n", nrow(out), paste(sp, collapse = ","), a[2]))
