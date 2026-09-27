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
#   data/DEGs/og_summary.RData            (per-orthogroup DEG summary; the common set is the
#                                          orthogroups with n_species_with_DEGs == 6)
#   data/DEGs/DEGs_col0.RData
#   data/annotation/Athal_go_ids.tsv      (gene -> GO mapping, Arabidopsis)
#   data/annotation/bg_col0_ATC.rda       (expressed-gene background, Col-0)
# Output:
#   tables/TableS1_commonDEG_GO_classic.tsv
#   tables/Table_S1.xlsx   (the supplementary table: Legend sheet + "Table S1" sheet; padj to 4 s.f., p to 3)
#
# Usage: Rscript scripts/GO/make_TableS1.R

source("scripts/GO/go_table_utils.R")

ld <- function(f, o = NULL) { e <- new.env(); load(f, envir = e); if (is.null(o)) get(ls(e)[1], envir = e) else get(o, e) }
ogfile <- Sys.glob("data/*/Orthogroups_20240613.tsv")[1]
og <- read.delim(ogfile, header = TRUE, sep = "\t", quote = "\"", check.names = FALSE, stringsAsFactors = FALSE)

# --- common-DEG orthogroups (DEG in all six genotypes) ---
ogsf <- Sys.glob("data/DEGs/og_summary.RData")
if (length(ogsf)) {
  osum <- ld(ogsf[1])
  common <- osum$OrthoGroup[osum$n_species_with_DEGs == 6]
} else {
  stop("data/DEGs/og_summary.RData not found; it defines the common-DEG orthogroups ",
       "(n_species_with_DEGs == 6) and is provided in this repository under data/DEGs/.")
}
sub <- og[og[["Ortholog_Group"]] %in% common, ]
at <- unique(trimws(unlist(strsplit(paste(sub[["Arabidopsis_thaliana"]], collapse = ", "), "[,;[:space:]]+"))))
at <- at[grepl("^AT[0-9MC]G", at)]

studyset <- intersect(at, ld("data/DEGs/DEGs_col0.RData"))
geneID2GO <- read_go_mapping("data/annotation/Athal_go_ids.tsv")
g2 <- geneID2GO[names(geneID2GO) %in% ld("data/annotation/bg_col0_ATC.rda")]
gl <- factor(as.integer(names(g2) %in% studyset)); names(gl) <- names(g2)

GOdata <- new("topGOdata", ontology = "BP", allGenes = gl, annot = annFUN.gene2GO, gene2GO = g2)
res <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
sig <- go_result(GOdata, res)

dir.create("tables", showWarnings = FALSE)
write.table(sig, "tables/TableS1_commonDEG_GO_classic.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
xl <- sig[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "classicFisher", "padj")]
names(xl)[names(xl) == "classicFisher"] <- "classicFisher P"
write_supp_xlsx(xl, "tables/Table_S1.xlsx", "Table S1", sig_cols = c("classicFisher P", "padj"),
  legend = sprintf(paste("Table S1 GO enrichment of the %d commonly differentially expressed orthogroups (topGO",
                         "classic Fisher's exact test, biological process; padj < 0.05)."), length(common)))
cat(sprintf("Wrote tables/TableS1_commonDEG_GO_classic.tsv and tables/Table_S1.xlsx (%d significant BP terms; padj < 0.05)\n", nrow(sig)))
cat("Top 5 terms:\n"); print(head(sig[, c("GO.ID", "Term", "Annotated", "Significant", "padj")], 5), row.names = FALSE)
