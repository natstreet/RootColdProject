#!/usr/bin/env Rscript
# make_TableS1.R — Table S1: GO (biological-process) enrichment of the commonly differentially
# expressed orthogroups (those with a DEG in all six genotypes), using the Arabidopsis (Col-0)
# DEGs as the study set.
#
# Enrichment: topGO, **classic** Fisher's exact test, Benjamini-Hochberg adjustment, padj < 0.05.
# (All GO enrichment in the paper uses the classic algorithm.)
#
# Inputs (unpack data.tar.gz from figshare into ./data):
#   data/*/Orthogroups_20240613.tsv
#   data/DEGs/og_summary.RData            (orthogroup x genotype DEG summary; the 270 common
#                                          orthogroups are those with n_species_with_DEGs == 6.
#                                          This is the definition of the common set and is
#                                          provided in the figshare archive.)
#   data/DEGs/DEGs_col0.RData
#   data/annotation/Athal_go_ids.tsv      (gene -> GO mapping, Arabidopsis)
#   data/annotation/bg_col0_ATC.rda       (expressed-gene background, Col-0)
# Output:
#   tables/TableS1_commonDEG_GO_classic.tsv
#
# Usage: Rscript scripts/GO/make_TableS1.R

suppressMessages({ library(topGO); library(GO.db) })

ld <- function(f, o = NULL) { e <- new.env(); load(f, envir = e); if (is.null(o)) get(ls(e)[1], envir = e) else get(o, e) }
ogfile <- Sys.glob("data/*/Orthogroups_20240613.tsv")[1]
og <- read.delim(ogfile, header = TRUE, sep = "\t", quote = "\"", check.names = FALSE, stringsAsFactors = FALSE)

# --- common-DEG orthogroups (DEG in all six genotypes) ---
ogsf <- Sys.glob("data/DEGs/og_summary.RData")
if (length(ogsf)) {
  osum <- ld(ogsf[1])
  common <- osum$OrthoGroup[osum$n_species_with_DEGs == 6]
} else {
  stop("data/DEGs/og_summary.RData not found. It defines the 270 common-DEG orthogroups ",
       "(n_species_with_DEGs == 6) and is included in the figshare data archive; unpack it into ",
       "data/DEGs/.")
}
sub <- og[og[["Ortholog_Group"]] %in% common, ]
at <- unique(trimws(unlist(strsplit(paste(sub[["Arabidopsis_thaliana"]], collapse = ", "), "[,;[:space:]]+"))))
at <- at[grepl("^AT[0-9MC]G", at)]

studyset <- intersect(at, ld("data/DEGs/DEGs_col0.RData"))
geneID2GO <- readMappings("data/annotation/Athal_go_ids.tsv")[-1]
g2 <- geneID2GO[names(geneID2GO) %in% ld("data/annotation/bg_col0_ATC.rda")]
gl <- factor(as.integer(names(g2) %in% studyset)); names(gl) <- names(g2)

GOdata <- new("topGOdata", ontology = "BP", allGenes = gl, annot = annFUN.gene2GO, gene2GO = g2)
res <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
tab <- GenTable(GOdata, classicFisher = res, orderBy = "classicFisher", topNodes = length(usedGO(GOdata)))
tab$classicFisher <- as.numeric(tab$classicFisher)
tab$padj <- round(p.adjust(tab$classicFisher, "BH"), 4)
tab <- tab[order(tab$padj), ]
sig <- tab[tab$padj < 0.05, ]

dir.create("tables", showWarnings = FALSE)
write.table(sig, "tables/TableS1_commonDEG_GO_classic.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote tables/TableS1_commonDEG_GO_classic.tsv (%d significant BP terms; padj < 0.05)\n", nrow(sig)))
cat("Top 5 terms:\n"); print(head(sig[, c("GO.ID", "Term", "Annotated", "Significant", "padj")], 5), row.names = FALSE)
