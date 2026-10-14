#!/usr/bin/env Rscript
# de_threshold_objects.R — the deposit's threshold-test DE objects, from the deposited DESeq2 objects (ledger F032: they were
# carried over from an earlier analysis, and the Ost-0 ones no longer matched its rebuilt object).
#   DEGs/padj_<code>.RData                     padj_df_<code>: per gene (rows in the DESeq2 object's order) padj_<tp> for each
#                                              cold time point against the control, Wald test against a two-fold null
#                                              (lfcThreshold = 1), and MinPadj, the smallest of them (NA when any is NA)
#   per_timepoint_DEG_lists/<g>/DEGs_<code>_C_vs_<tp>_p0.05_l2fc1.txt     padj < 0.05 in that test
#   per_timepoint_DEG_lists/<g>/DEGs_<code>_C_vs_<tp>_p0.05_lfc0585.txt   padj < 0.05 against a 1.5-fold null (lfcThreshold =
#                                              log2(1.5)); the earlier lists of this name followed no recorded rule
# Lists: genes with a positive log2 fold change first, then negative, each in gene-id order (en_US collation, as the
# earlier lists).
# DESeq() with default settings; alpha = 0.05 for independent filtering.
# Usage: Rscript de_threshold_objects.R <deposit root> <out data dir>
suppressMessages({ library(DESeq2); library(data.table) })
invisible(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))   # gene-id order as in the earlier lists (case-insensitive collation)
a <- commandArgs(trailingOnly = TRUE); R <- a[1]; OUT <- a[2]
G <- list(col0 = c("dds_col0_soil.rda", "col0"), ost0 = c("dds_ost0_soil.rda", "ost0"), aspen = c("dds_Pt_new.rda", "pt"),
          birch = c("dds_Bp_new.rda", "bp"), spruce = c("dds_Pa.rda", "pa"), pine = c("dds_Ps.rda", "ps"))
TP <- c("6h", "24h", "3d", "10d")
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
dir.create(file.path(OUT, "DEGs"), recursive = TRUE, showWarnings = FALSE)
for (g in names(G)) {
  code <- G[[g]][2]; dds <- suppressMessages(DESeq(ld(file.path(R, "data/dds", G[[g]][1])), quiet = TRUE))
  v <- all.vars(design(dds))[1]; lv <- levels(colData(dds)[[v]])
  pd <- data.frame(row.names = rownames(dds)); dir.create(file.path(OUT, "DEGs/per_timepoint_DEG_lists", g), recursive = TRUE, showWarnings = FALSE)
  for (tp in TP) {
    hit <- lv[lv == tp | sub("_5C$", "", lv) == tp]; stopifnot(length(hit) == 1)
    for (thr in c(l2fc1 = 1, lfc0585 = log2(1.5))) {
      r <- results(dds, contrast = c(v, hit, "control"), lfcThreshold = thr, alpha = 0.05)
      if (thr == 1) pd[[paste0("padj_", tp)]] <- r$padj
      sel <- !is.na(r$padj) & r$padj < 0.05; up <- sort(rownames(r)[sel & r$log2FoldChange > 0]); dn <- sort(rownames(r)[sel & r$log2FoldChange < 0])
      tag <- if (thr == 1) "l2fc1" else "lfc0585"
      writeLines(c(up, dn), file.path(OUT, "DEGs/per_timepoint_DEG_lists", g, sprintf("DEGs_%s_C_vs_%s_p0.05_%s.txt", code, tp, tag)))
    }
  }
  pd$MinPadj <- do.call(pmin, unname(as.list(pd)))
  nm <- paste0("padj_df_", code); assign(nm, pd); save(list = nm, file = file.path(OUT, sprintf("DEGs/padj_%s.RData", code)))
  cat(sprintf("%-6s padj object and %d lists written\n", g, 2 * length(TP)))
}
