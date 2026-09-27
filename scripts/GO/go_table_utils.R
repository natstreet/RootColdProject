# go_table_utils.R — shared by make_TableS1.R, make_TableS3.R and make_TableS4.R.
#
# go_result():      the significant terms of one topGO classic Fisher test, with full GO term names (GenTable
#                   cuts names at 40 characters), the exact p-value and the Benjamini-Hochberg padj. Terms are
#                   selected on the unrounded padj (< 0.05) and ordered by padj, then p.
# write_supp_xlsx(): a supplementary table in the journal layout — a "Legend" sheet holding the legend and a
#                   "Table Sn" sheet holding the table — with padj given to four significant figures and
#                   p-values to three.

suppressMessages({ library(topGO); library(GO.db) })

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
  # writexl writes only the parts a workbook uses. The Legend sheet is a one-cell sheet: a zero-row table whose
  # only header is the legend text.
  leg <- stats::setNames(data.frame(matrix(nrow = 0, ncol = 1)), legend)
  sheets <- list(Legend = leg); sheets[[sheet]] <- tab
  writexl::write_xlsx(sheets, file, format_headers = FALSE)
  check_xlsx_parts(file)
}

# Stop if any relationship target in the workbook's .rels files is missing from the archive (a dangling reference
# makes Excel offer to repair the file and openpyxl refuse to open it).
check_xlsx_parts <- function(file) {
  parts <- utils::unzip(file, list = TRUE)$Name
  rels <- grep("(^|/)_rels/[^/]*[.]rels$", parts, value = TRUE)
  ex <- tempfile(); on.exit(unlink(ex, recursive = TRUE)); utils::unzip(file, files = rels, exdir = ex)
  resolve <- function(base, t) {                                # join and collapse "." and ".." segments
    seg <- c(if (startsWith(t, "/") || base == ".") character(0) else strsplit(base, "/")[[1]], strsplit(t, "/")[[1]])
    out <- character(0)
    for (x in seg[seg != "" & seg != "."]) out <- if (x == "..") head(out, -1) else c(out, x)
    paste(out, collapse = "/")
  }
  bad <- character(0)
  for (r in rels) {
    x <- paste(readLines(file.path(ex, r), warn = FALSE), collapse = "")
    rel <- regmatches(x, gregexpr("<Relationship [^>]*>", x))[[1]]
    rel <- rel[!grepl('TargetMode="External"', rel)]
    tg <- sub('^.*Target="([^"]+)".*$', "\\1", rel)
    base <- dirname(dirname(r))                                  # the folder the .rels file describes
    for (t in tg) if (!resolve(base, t) %in% parts) bad <- c(bad, sprintf("%s -> %s", r, t))
  }
  if (length(bad)) stop("xlsx has dangling part references: ", paste(bad, collapse = "; "))
  invisible(TRUE)
}
