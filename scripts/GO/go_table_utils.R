# go_table_utils.R — shared by make_TableS1.R, make_TableS2.R and make_TableS3.R.
#
# go_result():      the significant terms of one topGO classic Fisher test, with full GO term names (GenTable
#                   cuts names at 40 characters), the exact p-value and the Benjamini-Hochberg padj. Terms are
#                   selected on the unrounded padj (< 0.05) and ordered by padj, then p.
# read_go_mapping(): a gene -> GO mapping file of the form "gene<tab>GO:a,GO:b", with or without a header line
#                   (Potra02_go_ids.tsv has one; Athal_go_ids.tsv and Betpen_go_ids.tsv do not), checked so that
#                   every gene in the file is in the mapping.
# write_supp_xlsx(): a supplementary table in the journal layout — a "Legend" sheet holding the legend and a
#                   "Table Sn" sheet holding the table — with padj given to four significant figures and
#                   p-values to three.

suppressMessages({ library(topGO); library(GO.db) })

read_go_mapping <- function(file) {
  first <- strsplit(readLines(file, n = 1), "\t")[[1]]
  has_header <- length(first) < 2 || !grepl("^GO:", trimws(first[2]))
  m <- topGO::readMappings(file)
  if (has_header) m <- m[-1]                                # readMappings reads a header line as a gene
  genes <- unique(sub("\t.*$", "", readLines(file)))
  if (has_header) genes <- genes[-1]
  genes <- genes[nzchar(genes)]
  if (!setequal(names(m), genes))
    stop(file, ": the mapping holds ", length(m), " genes, the file ", length(genes))
  m
}

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
  # tab may be one table (written to `sheet`) or a named list of tables (one sheet each; `sheet` is then ignored)
  tabs <- if (is.data.frame(tab)) stats::setNames(list(tab), sheet) else tab
  for (s in names(tabs)) for (k in intersect(sig_cols, names(tabs[[s]])))
    tabs[[s]][[k]] <- signif(tabs[[s]][[k]], if (startsWith(k, "padj")) 4 else 3)
  # writexl writes only the parts a workbook uses. The Legend sheet is a one-cell sheet: a zero-row table whose
  # only header is the legend text.
  leg <- stats::setNames(data.frame(matrix(nrow = 0, ncol = 1)), legend)
  sheets <- c(list(Legend = leg), tabs)
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

# size_matched_p(): exact P of observing at least `obs` term-carrying units among a study set drawn at random within size
#                   strata (stratum s: n[s] units drawn without replacement from N[s], of which K[t, s] carry term t), i.e.
#                   the distribution of a sum of hypergeometric counts, by convolution. The exact form of a size-matched
#                   random-draw null, without Monte Carlo error. Returns P and the expected count per term.
size_matched_p <- function(obs, K, N, n) {
  stopifnot(length(obs) == nrow(K), ncol(K) == length(N), length(N) == length(n), all(n <= N))
  res <- t(vapply(seq_along(obs), function(t) {
    pmf <- 1
    for (s in seq_along(N)) if (n[s] > 0) pmf <- convolve(pmf, rev(dhyper(0:n[s], K[t, s], N[s] - K[t, s], n[s])), type = "open")
    pmf <- pmax(pmf, 0); pmf <- pmf / sum(pmf)
    c(p = min(1, sum(pmf[seq_along(pmf) - 1 >= obs[t]])), expected = sum(n * K[t, ] / N))
  }, c(p = 0, expected = 0)))
  data.frame(p = res[, "p"], expected = res[, "expected"])
}
