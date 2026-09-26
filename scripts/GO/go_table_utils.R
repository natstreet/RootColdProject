# go_table_utils.R — shared by make_TableS1.R, make_TableS2.R and make_TableS3.R.
#
# go_result():      the significant terms of one topGO classic Fisher test, with full GO term names (GenTable
#                   cuts names at 40 characters), the exact p-value and the Benjamini-Hochberg padj. Terms are
#                   selected on the unrounded padj (< 0.05) and ordered by padj, then p.
# write_supp_xlsx(): a supplementary table in the journal layout — a "Legend" sheet holding the legend and a
#                   "Table Sn" sheet holding the table — with padj given to four significant figures and
#                   p-values to three.

suppressMessages({ library(topGO); library(GO.db); library(openxlsx) })

go_result <- function(GOdata, res, alpha = 0.05) {
  tab <- GenTable(GOdata, classicFisher = res, orderBy = "classicFisher", topNodes = length(usedGO(GOdata)))
  tab$Term <- unname(AnnotationDbi::Term(GO.db::GOTERM[tab$GO.ID]))
  tab$classicFisher <- unname(score(res)[tab$GO.ID])     # exact p (GenTable prints p < 1e-30 as text)
  tab$padj <- p.adjust(tab$classicFisher, "BH")
  tab <- tab[order(tab$padj, tab$classicFisher, tab$GO.ID), ]
  tab[tab$padj < alpha, ]
}

write_supp_xlsx <- function(tab, file, sheet, legend, sig_cols) {
  # padj to four significant figures, so that a value just below the 0.05 cut-off (e.g. 0.04995) does not
  # display as 0.05; the other statistics to three
  for (k in sig_cols) tab[[k]] <- signif(tab[[k]], if (k == "padj") 4 else 3)
  wb <- createWorkbook()
  addWorksheet(wb, "Legend"); writeData(wb, "Legend", legend, colNames = FALSE)
  addWorksheet(wb, sheet); writeData(wb, sheet, tab)
  setColWidths(wb, sheet, cols = seq_along(tab), widths = "auto")
  saveWorkbook(wb, file, overwrite = TRUE)
}
