#!/usr/bin/env Rscript
# check_padj_objects.R — do the deposited data/DEGs/padj_<code>.RData objects follow from the deposited DESeq2 objects?
# For each genotype: DESeq() on data/dds (default settings), then for each cold time point the Wald test against the
# control with lfcThreshold = 1 (the "two-fold null" the data dictionary describes for padj). Compared with the deposited
# padj_<tp> columns (NA compared as a pattern: genes DESeq2 filters out). Method check: a genotype whose counts did not change must
# reproduce its object exactly before a difference elsewhere means anything.
# Usage: Rscript check_padj_objects.R <unpacked deposit root> <out.tsv>
suppressMessages({ library(DESeq2); library(data.table) })
a <- commandArgs(trailingOnly = TRUE); R <- a[1]; OUT <- a[2]
G <- list(col0 = c("dds_col0_soil.rda", "col0"), ost0 = c("dds_ost0_soil.rda", "ost0"), aspen = c("dds_Pt_new.rda", "pt"),
          birch = c("dds_Bp_new.rda", "bp"), spruce = c("dds_Pa.rda", "pa"), pine = c("dds_Ps.rda", "ps"))
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
rows <- list()
for (g in names(G)) {
  dds <- ld(file.path(R, "data/dds", G[[g]][1])); old <- ld(file.path(R, "data/DEGs", sprintf("padj_%s.RData", G[[g]][2])))
  v <- all.vars(design(dds))[1]; dds <- suppressMessages(DESeq(dds, quiet = TRUE))
  for (col in grep("^padj_", names(old), value = TRUE)) {
    tp <- sub("^padj_", "", col)
    lv <- levels(colData(dds)[[v]]); hit <- lv[lv == tp | sub("_5C$", "", lv) == tp]
    if (length(hit) != 1) { rows[[length(rows) + 1]] <- data.table(genotype = g, column = col, status = "no matching level"); next }
    r <- results(dds, contrast = c(v, hit, "control"), lfcThreshold = 1, alpha = 0.05)
    p <- r$padj; names(p) <- rownames(r)
    o <- old[[col]]; names(o) <- rownames(old); gg <- intersect(names(p), names(o))
    both <- gg[!is.na(p[gg]) & !is.na(o[gg])]; d <- abs(p[both] - o[both])
    rows[[length(rows) + 1]] <- data.table(genotype = g, column = col, level = hit, genes = length(gg),
      na_old = sum(is.na(o[gg])), na_new = sum(is.na(p[gg])), na_pattern_differs = sum(is.na(o[gg]) != is.na(p[gg])),
      max_abs_diff = if (length(both)) max(d) else NA_real_, cells_differing_1e8 = sum(d > 1e-8),
      old_below_0.05 = sum(o[gg] < 0.05, na.rm = TRUE), new_below_0.05 = sum(p[gg] < 0.05, na.rm = TRUE))
  }
}
out <- rbindlist(rows, fill = TRUE); fwrite(out, OUT, sep = "\t"); print(out)
