#!/usr/bin/env Rscript
# ga_biosynthesis_check.R — the values of the gibberellin (GA) paragraph (DECISIONS.md D049; ledger R037 item 1).
#   (a) for OG0000184 and OG0000925: the lowest Max.p.val of any gene pair in the col0-ost0, aspen-birch and spruce-pine
#       comparison tables (ComPlEx, abspearsonMR0.03)
#   (b) Norway spruce genes annotated GO:0009686 (gibberellin biosynthetic process) in data/annotation/Picab02_go_ids.csv:
#       annotated; expressed (in the spruce GO background, data/annotation/bg_spruce_PA.rda: nonzero VST in the cold time
#       course); cold DEGs (union of the spruce per-time-point DEG lists in data/DEGs/per_timepoint_DEG_lists/spruce/);
#       of those, down (sign of the largest |log2FC| over the four time points, TopL2FC in data/DEGs/log2FoldChange_pa.RData)
# Output: results/ga_biosynthesis.tsv (part, item, value). Usage (working tree): Rscript scripts/analysis/ga_biosynthesis_check.R
suppressPackageStartupMessages(library(data.table))
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
res <- list()
for (cmp in c("col0-ost0", "aspen-birch", "spruce-pine")) {
  x <- as.data.table(ld(sprintf("data/ComPlEx/PairwiseComparisons/comparison-%s-abspearsonMR0.03-table.RData", cmp)))
  for (og in c("OG0000184", "OG0000925")) {
    p <- suppressWarnings(min(as.numeric(x[OrthoGroup == og, Max.p.val])))
    res[[length(res) + 1]] <- data.table(part = "a", item = sprintf("%s lowest Max.p.val %s", og, cmp), value = if (is.finite(p)) format(signif(p, 3)) else "no pair tested") } }
go <- fread("data/annotation/Picab02_go_ids.csv")
ga <- unique(go[go_id == "GO:0009686", external_gene_name])
bg <- ld("data/annotation/bg_spruce_PA.rda")
degs <- unique(unlist(lapply(list.files("data/DEGs/per_timepoint_DEG_lists/spruce", "_p0\\.05_lfc0\\.txt$", full.names = TRUE), readLines)))
if (!setequal(degs, ld("data/DEGs/DEGs_spruce.RData"))) stop("the per-time-point lfc0 lists do not equal DEGs_spruce")
ex <- fread("data/expression/spruce_root_cold_stress_expression.txt"); m <- as.matrix(ex[, -1]); rownames(m) <- ex[[1]]
cpx <- rownames(m)[apply(m, 1, sd) > 0 & rowSums(m != 0) >= 3 & rowMeans(m) > 1]
lf <- ld("data/DEGs/log2FoldChange_pa.RData")
top <- vapply(lf$TopL2FC, function(v) if (length(v) && !is.na(v[1])) unname(v[1]) else NA_real_, numeric(1)); names(top) <- rownames(lf)
ga_deg <- intersect(intersect(ga, bg), degs)
res[[length(res) + 1]] <- data.table(part = "b", item = c("GO:0009686 spruce genes annotated", "expressed (GO background)", "expressed (ComPlEx filters)",
  "in the expression matrix", "cold DEGs", "DEGs down (largest |log2FC| negative)"),
  value = as.character(c(length(ga), length(intersect(ga, bg)), length(intersect(ga, cpx)), length(intersect(ga, rownames(m))), length(ga_deg), sum(top[ga_deg] < 0, na.rm = TRUE))))
out <- rbindlist(res); fwrite(out, "results/ga_biosynthesis.tsv", sep = "\t"); print(out, row.names = FALSE)
