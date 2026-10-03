#!/usr/bin/env Rscript
# tf_clique_directions.R — cold response of each member of the transcription-factor cliques of
# Fig. 6, as tested in the Methods, for the expression claims in the Discussion.
#
# For each Fig. 6 transcription-factor clique (scripts/analysis/tf_cliques.R), for each of its six members:
#   * per cold timepoint, DESeq2 log2 fold change and adjusted P from the same Wald test as
#     scripts/DE/differential_expression.R (lfcThreshold = 0, default independent filtering), and the
#     Methods DEG call (padj <= 0.05 and |log2FC| >= 1), with its direction;
#   * the member's expression level within its own species: the percentile of its mean VST across the
#     control samples among all genes of that species (VST levels are not comparable between species).
#
# Inputs (run from a directory containing the unpacked figshare deposit `data/`, merged with this repo's
# data/):
#   data/dds/dds_<species>.rda, data/expression/<species>_root_cold_stress_expression.txt,
#   data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS, data/annotation/gene_aliases_20140331.txt,
#   scripts/analysis/tf_cliques.R (the transcription-factor cliques)
# Output: tf_clique_directions.tsv
# Usage:  Rscript scripts/analysis/tf_clique_directions.R

suppressMessages(library(DESeq2))

source("scripts/analysis/tf_cliques.R")
al  <- read.delim("data/annotation/gene_aliases_20140331.txt", stringsAsFactors = FALSE)
SP  <- c(col0 = "dds_col0_soil", ost0 = "dds_ost0_soil", aspen = "dds_Pt_new", birch = "dds_Bp_new",
         spruce = "dds_Pa", pine = "dds_Ps")
getdds <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
cl <- tf_cliques()                     # the Fig. 6 transcription-factor cliques (tf = panel label)
symbols <- cl$tf
sym_of <- function(g) { i <- which(al$locus_name == g); if (length(i)) al$symbol[i[1]] else "" }

out <- list()
for (sp in names(SP)) {
  dds <- DESeq(getdds(Sys.glob(sprintf("data/dds/%s.rda", SP[[sp]]))[1]), quiet = TRUE)
  dv  <- all.vars(design(dds))[1]
  levs <- setdiff(levels(colData(dds)[[dv]]), "control")
  ex  <- read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", sp), row.names = 1,
                    check.names = FALSE)
  ctl <- grep("\\.01$", colnames(ex))                                   # control timepoint samples
  pct <- rank(rowMeans(ex[, ctl, drop = FALSE])) / nrow(ex) * 100
  for (lv in levs) {
    r <- results(dds, contrast = c(dv, lv, "control"), alpha = 0.05, lfcThreshold = 0)
    for (i in seq_len(nrow(cl))) {
      g <- cl[[sp]][i]; if (is.na(g)) next
      l2 <- if (g %in% rownames(r)) r[g, "log2FoldChange"] else NA
      pa <- if (g %in% rownames(r)) r[g, "padj"] else NA
      deg <- !is.na(pa) && pa <= 0.05 && abs(l2) >= 1
      out[[length(out) + 1]] <- data.frame(
        tf = cl$tf[i], OrthoGroup = cl$OrthoGroup[i], genotype = sp, gene = g,
        gene_symbol = if (sp %in% c("col0", "ost0")) sym_of(g) else "", timepoint = sub("_.*$", "", lv),   # 6h_5C -> 6h
        log2FC = round(l2, 3), padj = signif(pa, 3), deg = deg,
        direction = if (deg) ifelse(l2 > 0, "up", "down") else "",
        control_expression_percentile = round(unname(pct[g]), 1))
    }
  }
}
res <- do.call(rbind, out)
res$genotype <- factor(res$genotype, names(SP)); res$timepoint <- factor(res$timepoint, c("6h", "24h", "3d", "10d"))
res <- res[order(match(res$tf, symbols), res$genotype, res$timepoint), ]
write.table(res, "tf_clique_directions.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("cliques:", paste(sprintf("%s=%s", cl$tf, cl$OrthoGroup), collapse = ", "), "\n")
cat("Col-0 and Ost-0 members differ in:", paste(cl$tf[cl$col0 != cl$ost0], collapse = ", "), "\n")
