#!/usr/bin/env Rscript
# build_coexpression_network.R — one whole-genome co-expression network from one experiment's raw counts.
#
# Method (chosen by the PlantGenIE network benchmark, 20 experiments; see docs/coexpression_networks.md):
#   1. Samples: those listed in the groups table, less any excluded; samples without reads are dropped.
#   2. Gene filter: keep genes with >= 1 count per million (edgeR::cpm on the raw counts) in at least k samples,
#      k = max(3, the size of the smallest group with at least two samples).
#   3. Normalisation: TMM scaling factors on the kept genes (edgeR::normLibSizes); log2(scaled count + 1), the counts
#      scaled to the mean effective library size.
#   4. Association: jackknifed Pearson correlation (Heyer et al. 1999). For each gene pair the correlation is
#      recomputed leaving out each sample in turn and the value closest to zero is kept (0 when the leave-one-out
#      values differ in sign), so that no single sample can create a link. Computed exactly from running sums.
#      Values are rounded to 10 decimal places, so equal values tie exactly.
#   5. Ranks: within each gene's row, partners ranked by the jackknifed correlation (ties averaged; the gene itself
#      excluded). With sign=positive (default) the strongest positive partner ranks 1; with sign=negative the strongest
#      negative partner ranks 1. Each network keeps only edges of its own sign, so the positive and the negative
#      network are two separate layers, each with its own ranks. With sign=absolute, partners are ranked by the
#      absolute jackknifed correlation and both signs share one network (edges with a jackknifed correlation of zero
#      are left out).
#   6. Edges: a pair is linked when its highest reciprocal rank, HRR = max(rank i->j, rank j->i), is at most HRR_MAX
#      (default 100) and its jackknifed correlation has the network's sign. Each edge carries both ranks, HRR, the mutual rank sqrt(rank i->j * rank j->i), the jackknifed
#      correlation and the ordinary Pearson correlation of the same normalised values.
#
# Inputs: counts TSV (first column gene id, one column per sample; rows with the same gene id are summed, so a
#         transcript-level table can be given with gene ids in the first column), groups TSV (sample, group)
# Output: <out>.tsv  source, target, rank_source_to_target, rank_target_to_source, hrr, mr, jackknife_r, pearson_r
#         <out>.summary.tsv  one line: samples, k, genes before and after the filter, edges, degree summary
# Usage:  Rscript build_coexpression_network.R <counts.tsv> <groups.tsv> <out prefix> [HRR_MAX=100] [exclude=a;b]
#                                              [sign=positive|negative|absolute]
suppressMessages({ library(data.table); library(edgeR); library(matrixStats) })
a <- commandArgs(trailingOnly = TRUE)
COUNTS <- a[1]; GROUPS <- a[2]; OUT <- a[3]
opt <- list(HRR_MAX = "100", exclude = "", block = "500", sign = "positive")
for (x in a[-(1:3)]) { kv <- strsplit(x, "=", fixed = TRUE)[[1]]; if (!kv[1] %in% names(opt)) stop("unknown option ", x); opt[[kv[1]]] <- kv[2] }
HRR_MAX <- as.integer(opt$HRR_MAX); BLOCK <- as.integer(opt$block)
if (!opt$sign %in% c("positive", "negative", "absolute")) stop("sign must be positive, negative or absolute")
SGN <- if (opt$sign == "negative") -1 else 1
ABS <- opt$sign == "absolute"
CPM_MIN <- 1

## ---- 1. counts and samples ----
t <- fread(COUNTS, header = TRUE)                     # sample names can be all digits: never guess the header
g <- t[[1]]; t[[1]] <- NULL
t <- t[, vapply(t, is.numeric, logical(1)), with = FALSE]   # drop text columns such as gene names
cnt <- rowsum(as.matrix(t), g); storage.mode(cnt) <- "double"
gt <- fread(GROUPS, colClasses = "character", header = TRUE); grp <- setNames(gt$group, gt$sample)
excl <- if (nzchar(opt$exclude)) strsplit(opt$exclude, ";")[[1]] else character(0)
keep_s <- setdiff(intersect(colnames(cnt), names(grp)), excl)
if (length(setdiff(names(grp), c(colnames(cnt), excl))))
  stop("samples in the groups table but not in the counts: ", paste(setdiff(names(grp), c(colnames(cnt), excl)), collapse = ", "))
cnt <- cnt[, keep_s, drop = FALSE]; cnt <- cnt[, colSums(cnt) > 0, drop = FALSE]; grp <- grp[colnames(cnt)]

## ---- 2. gene filter ----
gs <- table(grp); k <- max(3L, max(2L, min(gs[gs >= 2])))
pass <- rowSums(cpm(cnt) >= CPM_MIN) >= k
x <- cnt[pass, , drop = FALSE]

## ---- 3. TMM normalisation ----
f <- normLibSizes(DGEList(x), method = "TMM")$samples
eff <- f$lib.size * f$norm.factors
if (any(!is.finite(eff) | eff <= 0)) stop("TMM scaling factor zero or not finite")
x <- log2(t(t(x) / eff * mean(eff)) + 1)
n <- nrow(x); S <- ncol(x); genes <- rownames(x)
cat(sprintf("%s: %d samples, %d groups (smallest replicated %d; k = %d); %d of %d genes pass the filter\n",
            basename(OUT), S, length(gs), min(gs[gs >= 2]), k, n, nrow(cnt)))

## ---- 4-5. jackknifed correlation and within-row ranks, in row blocks ----
X <- x - rowMeans(x); Cxx <- rowSums(X^2); cf <- S / (S - 1)   # leaving sample s out changes each centred
Z <- X / sqrt(Cxx); Z[!is.finite(Z)] <- 0                        # cross-product by -S/(S-1) times its term at s
top_i <- top_j <- integer(0); top_r <- top_jk <- top_p <- numeric(0)
for (from in seq(1L, n, by = BLOCK)) {
  to <- min(from + BLOCK - 1L, n); b <- from:to
  Cxy <- tcrossprod(X[b, , drop = FALSE], X)
  lo <- matrix(Inf, length(b), n); hi <- matrix(-Inf, length(b), n)
  for (s in seq_len(S)) {
    xb <- X[b, s]; xa <- X[, s]
    vb <- Cxx[b] - cf * xb^2; va <- Cxx - cf * xa^2
    vb[vb <= 1e-10 * Cxx[b]] <- 0; va[va <= 1e-10 * Cxx] <- 0    # constant without sample s: no correlation
    r <- (Cxy - cf * outer(xb, xa)) / sqrt(outer(vb, va)); r[!is.finite(r)] <- 0
    lo <- pmin(lo, r); hi <- pmax(hi, r)
  }
  jk <- round(ifelse(lo > 0, lo, ifelse(hi < 0, hi, 0)), 10); rm(lo, hi, Cxy, r)
  jk[cbind(seq_along(b), b)] <- NA
  sc <- if (ABS) abs(jk) else SGN * jk                    # higher = stronger partner for this network
  R <- rowRanks(-sc, ties.method = "average")            # 1 = strongest partner
  w <- which(R <= HRR_MAX & sc > 0, arr.ind = TRUE)
  top_i <- c(top_i, b[w[, 1]]); top_j <- c(top_j, w[, 2]); top_r <- c(top_r, R[w]); top_jk <- c(top_jk, jk[w])
  top_p <- c(top_p, rowSums(Z[b[w[, 1]], , drop = FALSE] * Z[w[, 2], , drop = FALSE]))
  rm(jk, R, w)
  if ((to %/% BLOCK) %% 10 == 0) cat(sprintf("  %d of %d genes\n", to, n))
}

## ---- 6. reciprocal edges ----
d <- data.table(i = top_i, j = top_j, r = top_r, jk = top_jk, p = top_p)
e <- merge(d[i < j], d[i > j, .(i = j, j = i, r_back = r)], by = c("i", "j"))   # both directions within HRR_MAX
out <- e[, .(source = genes[i], target = genes[j], rank_source_to_target = r, rank_target_to_source = r_back,
             hrr = pmax(r, r_back), mr = round(sqrt(r * r_back), 2), jackknife_r = round(jk, 4), pearson_r = round(p, 4))]
setorder(out, source, hrr)
fwrite(out, paste0(OUT, ".tsv"), sep = "\t")
deg <- table(c(out$source, out$target))
summ <- data.table(network = basename(OUT), sign = opt$sign, samples = S, groups = length(gs), k = k, cpm_min = CPM_MIN,
                   genes_in_counts = nrow(cnt), genes_after_filter = n, hrr_max = HRR_MAX, edges = nrow(out),
                   genes_with_edges = length(deg), mean_degree = round(mean(as.numeric(deg)), 1),
                   median_degree = median(as.numeric(deg)), max_degree = max(as.numeric(deg)),
                   edges_jackknife_r_zero = sum(out$jackknife_r == 0))
fwrite(summ, paste0(OUT, ".summary.tsv"), sep = "\t")
cat(sprintf("%d edges among %d genes -> %s.tsv\n", nrow(out), length(deg), OUT))
