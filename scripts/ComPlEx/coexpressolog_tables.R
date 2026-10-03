#!/usr/bin/env Rscript
# coexpressolog_tables.R — the gene-pair tables built from the 15 ComPlEx pairwise comparisons:
#   all_expressed_genes.RData  every tested orthologue gene pair of every comparison (expr_genes)
#   co_expressologs.RData      the pairs with Max.p.val < 0.1, the co-expressologs (co_expressologs)
# Columns: OrthologGroup, Species1, GeneSpecies1, Species2, GeneSpecies2, [Species1pVal, Species2pVal,] MaxpVal,
# effectSize (= Min.effect.size). Replaces the corresponding chunks of Genelists.Rmd, which no longer runs as committed.
# Usage (repository root): Rscript scripts/ComPlEx/coexpressolog_tables.R [comparisons dir] [out dir]
suppressPackageStartupMessages(library(dplyr))
a <- commandArgs(trailingOnly = TRUE)
IN <- if (length(a) >= 1) a[1] else "data/ComPlEx/PairwiseComparisons"
OUT <- if (length(a) >= 2) a[2] else "data/ComPlEx"
PAIRS <- c("col0-ost0", "col0-aspen", "col0-birch", "col0-spruce", "col0-pine", "ost0-aspen", "ost0-birch", "ost0-spruce",
           "ost0-pine", "aspen-birch", "aspen-spruce", "aspen-pine", "birch-spruce", "birch-pine", "spruce-pine")
PVAL <- 0.1
rows <- lapply(PAIRS, function(p) {
  e <- new.env(); load(file.path(IN, sprintf("comparison-%s-abspearsonMR0.03-table.RData", p)), envir = e)
  ct <- as.data.frame(get("comparison_table", e)); sp <- strsplit(p, "-")[[1]]
  data.frame(OrthologGroup = ct$OrthoGroup, Species1 = sp[1], GeneSpecies1 = ct$Species1, Species2 = sp[2],
             GeneSpecies2 = ct$Species2, Species1pVal = ct$Species1.p.val, Species2pVal = ct$Species2.p.val,
             MaxpVal = ct$Max.p.val, effectSize = ct$Min.effect.size)
})
expr_genes <- do.call(rbind, rows)
expr_genes$MaxpVal <- as.numeric(expr_genes$MaxpVal)
# row order as in the deposited objects: all pairs by orthogroup, then pair (PAIRS order), then Max.p.val (ties keep
# the comparison table's order); co-expressologs by pair (alphabetical), then Max.p.val, keeping those row names
pair_idx <- match(paste(expr_genes$Species1, expr_genes$Species2, sep = "-"), PAIRS)
expr_genes <- expr_genes[order(expr_genes$OrthologGroup, pair_idx, expr_genes$MaxpVal, seq_len(nrow(expr_genes))), ]
rownames(expr_genes) <- NULL
co_expressologs <- expr_genes[expr_genes$MaxpVal < PVAL, c("OrthologGroup", "Species1", "GeneSpecies1",
                                                           "Species2", "GeneSpecies2", "MaxpVal", "effectSize")]
co_expressologs <- co_expressologs[order(paste(co_expressologs$Species1, co_expressologs$Species2, sep = "-"),
                                         co_expressologs$MaxpVal, as.integer(rownames(co_expressologs))), ]
save(expr_genes, file = file.path(OUT, "all_expressed_genes.RData"))
save(co_expressologs, file = file.path(OUT, "co_expressologs.RData"))
cat(sprintf("%d tested gene pairs; %d co-expressologs (Max.p.val < %s) -> %s\n", nrow(expr_genes), nrow(co_expressologs), PVAL, OUT))
