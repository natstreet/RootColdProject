#!/usr/bin/env Rscript
# check_rebuilt_objects.R — consistency of the regenerated analysis objects of the deposit with its DESeq2 objects.
#   (a) data/expression/<g>_root_cold_stress_expression.txt = VST (blind) of data/dds, minus the matrix minimum (sample
#       names through data/quantification/sample_name_map.tsv)
#   (b) DEGs_<g>.RData = union of the per-timepoint *_lfc0 lists = genes with is_deg in data/open/differential_expression.tsv
#   (c) log2FoldChange_<code>.RData = DESeq2 log2 fold change, time point against control (lfcThreshold 0)
#   (d) the per-timepoint *_lfc0 lists = padj (lfcThreshold 0) <= 0.05 and |log2FC| >= 1
#   (e) open/differential_expression.tsv padj_lfc0 and log2FoldChange = (c), (d); its padj column against the padj objects
# Usage: Rscript check_rebuilt_objects.R <unpacked deposit root> <out.tsv>
suppressMessages({ library(DESeq2); library(data.table) })
a <- commandArgs(trailingOnly = TRUE); R <- a[1]; OUT <- a[2]
G <- list(col0 = c("dds_col0_soil.rda", "col0"), ost0 = c("dds_ost0_soil.rda", "ost0"), aspen = c("dds_Pt_new.rda", "pt"),
          birch = c("dds_Bp_new.rda", "bp"), spruce = c("dds_Pa.rda", "pa"), pine = c("dds_Ps.rda", "ps"))
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
de <- fread(file.path(R, "data/open/differential_expression.tsv"))
nm <- fread(file.path(R, "data/quantification/sample_name_map.tsv"), colClasses = "character")
rows <- list(); add <- function(...) rows[[length(rows) + 1]] <<- data.table(...)
for (g in names(G)) {
  dds0 <- ld(file.path(R, "data/dds", G[[g]][1])); code <- G[[g]][2]
  # (a)
  v <- assay(varianceStabilizingTransformation(dds0, blind = TRUE)); v <- v - min(v)
  ex <- fread(file.path(R, sprintf("data/expression/%s_root_cold_stress_expression.txt", g)), header = TRUE)
  m <- as.matrix(ex[, -1]); rownames(m) <- ex[[1]]
  mp <- nm[genotype == g]; colnames(v) <- setNames(mp$deposited_sample_id, mp$sequencing_name)[colnames(v)]   # publication codes
  same_dim <- identical(sort(rownames(m)), sort(rownames(v))) && identical(sort(colnames(m)), sort(colnames(v)))
  add(genotype = g, check = "a VST matrix", result = if (same_dim) sprintf("max |diff| %.3g", max(abs(m[rownames(v), colnames(v)] - v))) else
      sprintf("names differ: genes %d/%d, samples %d/%d", length(intersect(rownames(m), rownames(v))), nrow(v), length(intersect(colnames(m), colnames(v))), ncol(v)))
  # DESeq2 at lfcThreshold 0
  dds <- suppressMessages(DESeq(dds0, quiet = TRUE)); vv <- all.vars(design(dds))[1]; lv <- setdiff(levels(colData(dds)[[vv]]), "control")
  lf <- ld(file.path(R, sprintf("data/DEGs/log2FoldChange_%s.RData", code)))
  union_new <- character(0); lf_dev <- 0; list_mis <- 0; de_dev_l <- 0; de_dev_p <- 0; de_deg_mis <- 0
  for (tp in lv) {
    t0 <- sub("_5C$", "", tp)
    r <- results(dds, contrast = c(vv, tp, "control"), alpha = 0.05)
    lf_dev <- max(lf_dev, max(abs(lf[rownames(r), paste0("log2FoldChange_", t0)] - r$log2FoldChange), na.rm = TRUE))
    deg <- rownames(r)[!is.na(r$padj) & r$padj <= 0.05 & abs(r$log2FoldChange) >= 1]; union_new <- union(union_new, deg)
    f <- file.path(R, "data/DEGs/per_timepoint_DEG_lists", g, sprintf("DEGs_%s_C_vs_%s_p0.05_lfc0.txt", code, t0))
    old <- readLines(f); old <- old[nzchar(old)]; list_mis <- list_mis + length(setdiff(old, deg)) + length(setdiff(deg, old))
    d <- de[species == g & sub("_5C$", "", timepoint) == t0]; setkey(d, gene)
    de_dev_l <- max(de_dev_l, max(abs(d[rownames(r), log2FoldChange] - r$log2FoldChange), na.rm = TRUE))
    de_dev_p <- max(de_dev_p, max(abs(d[rownames(r), padj_lfc0] - r$padj), na.rm = TRUE))
    de_deg_mis <- de_deg_mis + sum(xor(d[rownames(r), is_deg %in% c(TRUE, "True", "TRUE")], rownames(r) %in% deg))
  }
  dg <- ld(file.path(R, sprintf("data/DEGs/DEGs_%s.RData", g)))
  add(genotype = g, check = "b DEG union", result = sprintf("object %d, recomputed %d, differing %d", length(dg), length(union_new), length(setdiff(dg, union_new)) + length(setdiff(union_new, dg))))
  add(genotype = g, check = "c log2FoldChange object", result = sprintf("max |diff| %.3g", lf_dev))
  add(genotype = g, check = "d lfc0 per-timepoint lists", result = sprintf("genes differing over the time points %d", list_mis))
  add(genotype = g, check = "e open DE table", result = sprintf("log2FC max |diff| %.3g; padj_lfc0 max |diff| %.3g; is_deg calls differing %d", de_dev_l, de_dev_p, de_deg_mis))
  po <- ld(file.path(R, sprintf("data/DEGs/padj_%s.RData", code))); d <- de[species == g]
  pc <- d[, .(gene, tp = sub("_5C$", "", timepoint), padj)]; pc[, obj := mapply(function(gg, t) po[gg, paste0("padj_", t)], gene, tp)]
  add(genotype = g, check = "e open DE padj = padj object", result = sprintf("max |diff| %.3g; NA pattern differing %d", max(abs(pc$padj - pc$obj), na.rm = TRUE), sum(is.na(pc$padj) != is.na(pc$obj))))
}
out <- rbindlist(rows); fwrite(out, OUT, sep = "\t")
