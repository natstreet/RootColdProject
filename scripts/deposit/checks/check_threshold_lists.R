#!/usr/bin/env Rscript
# check_threshold_lists.R — do the deposited per-timepoint lists *_p0.05_l2fc1.txt and *_p0.05_lfc0585.txt follow from
# the deposited DESeq2 objects? Candidate rule, tested on every genotype: DESeq() with default settings, Wald test of the
# time point against the control with lfcThreshold = 1 (l2fc1) or log2(1.5) (lfc0585), genes with padj < 0.05.
# A genotype that reproduces exactly confirms the rule; a difference elsewhere is then a stale list.
# Usage: Rscript check_threshold_lists.R <unpacked deposit root> <out.tsv>
suppressMessages({ library(DESeq2); library(data.table) })
a <- commandArgs(trailingOnly = TRUE); R <- a[1]; OUT <- a[2]
G <- list(col0 = c("dds_col0_soil.rda", "col0"), ost0 = c("dds_ost0_soil.rda", "ost0"), aspen = c("dds_Pt_new.rda", "pt"),
          birch = c("dds_Bp_new.rda", "bp"), spruce = c("dds_Pa.rda", "pa"), pine = c("dds_Ps.rda", "ps"))
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
rows <- list()
for (g in names(G)) {
  dds <- suppressMessages(DESeq(ld(file.path(R, "data/dds", G[[g]][1])), quiet = TRUE)); v <- all.vars(design(dds))[1]; lv <- levels(colData(dds)[[v]])
  for (f in list.files(file.path(R, "data/DEGs/per_timepoint_DEG_lists", g), pattern = "_(l2fc1|lfc0585)\\.txt$", full.names = TRUE)) {
    tp <- sub(".*_C_vs_([^_]+)_p0\\.05_.*", "\\1", basename(f)); thr <- if (grepl("l2fc1", f)) 1 else log2(1.5)
    hit <- lv[lv == tp | sub("_5C$", "", lv) == tp]
    r <- results(dds, contrast = c(v, hit, "control"), lfcThreshold = thr, alpha = 0.05)
    new <- rownames(r)[!is.na(r$padj) & r$padj < 0.05]; old <- readLines(f); old <- old[nzchar(old)]
    rows[[length(rows) + 1]] <- data.table(genotype = g, file = basename(f), old_genes = length(old), new_genes = length(new),
                                           only_old = length(setdiff(old, new)), only_new = length(setdiff(new, old)))
  }
}
out <- rbindlist(rows); fwrite(out, OUT, sep = "\t")
