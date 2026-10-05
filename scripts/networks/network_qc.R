#!/usr/bin/env Rscript
# network_qc.R — checks for a network from build_coexpression_network.R, to run on every new network before it is
# hosted. Recomputes the filter and normalisation from the same counts and groups, then reports:
#   genes, edges, degree        as in the summary file, recounted from the edge file
#   rank range                  every edge has 1 <= HRR <= HRR_MAX
#   sign                        every edge's jackknifed correlation has the network's sign (none zero)
#   filter                      every gene in the network passes the expression filter
#   single-sample hubs          among the 30 genes with most links (ties: lower summed HRR), those with more than half
#                               of their variance (sum of squared deviations of the normalised values) in one sample
# and stops (exit 1) when a check fails. Single-sample hubs are reported, not failed: the benchmark found 0 to 3 per
# network (median 0) after the jackknife; more than 5 is flagged for a look at the gene's expression.
# Usage: Rscript network_qc.R <counts.tsv> <groups.tsv> <network.tsv> [positive|negative] [HRR_MAX=100] [exclude=a;b]
suppressMessages({ library(data.table); library(edgeR) })
a <- commandArgs(trailingOnly = TRUE)
SIGN <- if (length(a) >= 4) a[4] else "positive"; HRR_MAX <- if (length(a) >= 5) as.numeric(sub("HRR_MAX=", "", a[5])) else 100
excl <- if (length(a) >= 6) strsplit(sub("^exclude=", "", a[6]), ";")[[1]] else character(0)
t <- fread(a[1], header = TRUE); g <- t[[1]]; t[[1]] <- NULL
cnt <- rowsum(as.matrix(t[, vapply(t, is.numeric, TRUE), with = FALSE]), g); storage.mode(cnt) <- "double"
gt <- fread(a[2], colClasses = "character"); grp <- setNames(gt$group, gt$sample)
cnt <- cnt[, setdiff(intersect(colnames(cnt), names(grp)), excl), drop = FALSE]; cnt <- cnt[, colSums(cnt) > 0]; grp <- grp[colnames(cnt)]
gs <- table(grp); k <- max(3, min(gs[gs >= 2]))
x <- cnt[rowSums(cpm(cnt) >= 1) >= k, , drop = FALSE]
f <- normLibSizes(DGEList(x), method = "TMM")$samples; eff <- f$lib.size * f$norm.factors
x <- log2(t(t(x) / eff * mean(eff)) + 1)
e <- fread(a[3])
deg <- table(c(e$source, e$target)); hs <- tapply(c(e$hrr, e$hrr), c(e$source, e$target), sum)
hub <- names(deg)[order(-as.numeric(deg), hs[names(deg)])][seq_len(min(30, length(deg)))]
lev <- apply(x[hub, , drop = FALSE], 1, function(v) { d <- (v - mean(v))^2; max(d) / sum(d) })
fails <- c(
  rank = sum(e$hrr < 1 | e$hrr > HRR_MAX),
  sign = sum(if (SIGN == "positive") e$jackknife_r <= 0 else e$jackknife_r >= 0),
  filter = sum(!unique(c(e$source, e$target)) %in% rownames(x)))
cat(sprintf("%s (%s): %d samples, k = %d, %d genes pass the filter; %d edges among %d genes; degree median %g, max %d\n",
            basename(a[3]), SIGN, ncol(x), k, nrow(x), nrow(e), length(deg), median(as.numeric(deg)), max(deg)))
cat(sprintf("  edges outside HRR 1-%g: %d; edges of the wrong sign: %d; genes failing the filter: %d\n",
            HRR_MAX, fails[["rank"]], fails[["sign"]], fails[["filter"]]))
cat(sprintf("  single-sample genes among the 30 most linked: %d%s\n", sum(lev > 0.5),
            if (sum(lev > 0.5) > 5) " (more than 5: look at these genes' expression)" else ""))
if (sum(lev > 0.5)) cat("   ", paste(names(lev)[lev > 0.5], collapse = " "), "\n")
if (any(fails > 0)) { cat("FAILED\n"); quit(status = 1) }
cat("passed\n")
