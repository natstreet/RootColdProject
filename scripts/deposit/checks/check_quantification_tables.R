#!/usr/bin/env Rscript
# check_quantification_tables.R — the gene counts of data/quantification/ must equal the deposited DESeq2 counts, for every
# gene and every sample of every genotype (exactly; quantification_tables.py follows the counts' rounding convention), and
# lane_map.tsv must list a run file for every sample. Stops on any difference. Used by the deposit build as a check.
# Usage: Rscript check_quantification_tables.R <deposit root> [quantification dir, default <root>/data/quantification]
suppressMessages({ library(DESeq2); library(data.table) })
a <- commandArgs(trailingOnly = TRUE); R <- a[1]; Q <- if (length(a) > 1) a[2] else file.path(R, "data/quantification")
G <- c(col0 = "dds_col0_soil.rda", ost0 = "dds_ost0_soil.rda", aspen = "dds_Pt_new.rda", birch = "dds_Bp_new.rda",
       spruce = "dds_Pa.rda", pine = "dds_Ps.rda")
lm <- fread(file.path(Q, "lane_map.tsv"), colClasses = "character"); fail <- FALSE
for (g in names(G)) {
  e <- new.env(); load(file.path(R, "data/dds", G[[g]]), envir = e); cn <- counts(get(ls(e)[1], e))
  q <- fread(cmd = paste("gzip -dc", shQuote(file.path(Q, sprintf("%s_gene_counts.tsv.gz", g)))), header = TRUE)
  x <- as.matrix(q[, -1]); rownames(x) <- q[[1]]
  ok_names <- setequal(rownames(x), rownames(cn)) && setequal(colnames(x), colnames(cn))
  d <- if (ok_names) abs(x[rownames(cn), colnames(cn)] - cn) else NA
  runs <- lm[genotype == g & used == "yes", .N, by = deposited_sample_id]
  no_run <- setdiff(colnames(cn), runs$deposited_sample_id)
  ok <- ok_names && max(d) == 0 && !length(no_run)
  cat(sprintf("%-6s genes %d, samples %d; names match %s; max |count diff| %s; samples without a run file %d -> %s\n", g, nrow(cn), ncol(cn),
              ok_names, if (ok_names) max(d) else "n/a", length(no_run), if (ok) "OK" else "FAIL"))
  if (!ok) fail <- TRUE
}
if (fail) stop("quantification tables do not reproduce the deposited counts")
