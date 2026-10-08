#!/usr/bin/env Rscript
# design_stability.R — how sampling design sets the leave-one-sample-out stability of co-expression neighbourhoods
# (NRS, 2026-10-08: is the root data's sample sensitivity a matter of replicates per sampling point, of the number of
# sampling points, or of the total; does the wood data hold up better?). Informs future designs; ComPlEx is not re-run.
#
# Samples are named <replicate>.<point>: root <genotype><replicate>.<time point 01-05>, wood <tree>.<section>.
# For each design (r replicates per point, k points) and draw (10, seeds 1-10):
#   replicates  r drawn at random among the dataset's replicates (trees, or plants per time point)
#   points      root: control plus k-1 other time points at random; wood: k sections evenly spaced over the sections
#               present in all r chosen trees, with a random start (sections are positions along the developmental gradient)
#   samples     every chosen replicate at every chosen point that exists (n <= r x k)
# Stability, as stability_matched_n.R: for 300 random genes (seed 1) of the dataset's filtered universe (SD > 0, at least
# 3 non-zero values, mean > 1 on the full data; fixed for all designs), the top 3% of genes by |Pearson r|; Jaccard between
# the neighbourhood from the design's samples and from those samples minus one; per gene the worst (minimum) and the mean
# over left-out samples; per draw the median over genes.
# Input values: root, the ComPlEx expression inputs as they are; wood, log2(CPM + 1) of the raw salmon counts, without
# the samples not sequenced or excluded as outliers in the paper (README of the counts folder).
# Usage: Rscript design_stability.R <name> <matrix> <none|log2cpm> <root|wood> <r values, comma> <k values, comma> <out tsv> [exclude, comma]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
NAME <- a[1]; MAT <- a[2]; TR <- a[3]; KIND <- a[4]; RS <- as.integer(strsplit(a[5], ",")[[1]]); KS <- as.integer(strsplit(a[6], ",")[[1]]); OUT <- a[7]
EXC <- if (length(a) > 7) strsplit(a[8], ",")[[1]] else character(0)
DENS <- 0.03; NDRAW <- 10; NG <- 300
x <- fread(cmd = if (grepl("\\.gz$", MAT)) paste("gzip -dc", shQuote(MAT)) else paste("cat", shQuote(MAT)))
m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m <- m[, !colnames(m) %in% EXC & colSums(is.na(m)) == 0, drop = FALSE]
if (TR == "log2cpm") m <- log2(t(t(m) / colSums(m)) * 1e6 + 1)
m <- m[apply(m, 1, sd) > 0 & rowSums(m != 0) >= 3 & rowMeans(m) > 1, ]
set.seed(1); q <- sample(rownames(m), NG); qi <- match(q, rownames(m)); K3 <- ceiling(DENS * (nrow(m) - 1))
rep_of <- sub("\\..*", "", colnames(m)); pt_of <- as.integer(sub(".*\\.", "", colnames(m)))
topsets <- function(mm) { z <- t(scale(t(mm))); z[!is.finite(z)] <- 0
  r <- abs(tcrossprod(z[qi, , drop = FALSE], z)) / (ncol(mm) - 1)
  lapply(seq_along(qi), function(i) { v <- r[i, ]; v[qi[i]] <- -Inf; thr <- -sort(-v, partial = K3)[K3]; which(v >= thr) }) }
jac <- function(x, y) length(intersect(x, y)) / length(union(x, y))
res <- list()
for (R in RS) for (K in KS) for (d in seq_len(NDRAW)) {
  set.seed(d); reps <- sort(unique(rep_of)); if (R > length(reps)) next
  ch <- sample(reps, R)
  common <- sort(Reduce(intersect, lapply(ch, function(r) unique(pt_of[rep_of == r]))))
  if (KIND == "root") { if (K > length(common)) next; pts <- c(common[1], sort(sample(common[-1], K - 1)))
  } else { if (K > length(common)) next; step <- length(common) / K; pts <- common[floor((0:(K - 1)) * step + runif(1, 0, step)) + 1] }
  cols <- which(rep_of %in% ch & pt_of %in% pts); if (length(cols) < 4) next
  mm <- m[, cols]; full <- topsets(mm)
  J <- sapply(seq_len(ncol(mm)), function(s) mapply(jac, full, topsets(mm[, -s])))
  res[[length(res) + 1]] <- data.table(dataset = NAME, r = R, k = K, n = length(cols), draw = d,
    median_worst_jaccard = median(apply(J, 1, min)), median_mean_jaccard = median(rowMeans(J)))
  cat(NAME, R, K, d, length(cols), "\n")
}
fwrite(rbindlist(res), OUT, sep = "\t")
