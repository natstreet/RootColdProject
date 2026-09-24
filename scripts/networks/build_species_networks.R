#!/usr/bin/env Rscript
# build_species_networks.R — one whole-genome co-expression network per species.
#
# These are the per-species networks of the resource (data/networks/network_<species>_hrr100.tsv):
# for any gene they give its co-expression neighbourhood across the whole genome. They are distinct
# from the networks behind the pairwise co-expressolog comparisons, which ComPlEx.R builds for each
# species pair on the genes with an orthologue in that pair (top 3% of edges by mutual rank).
#
# Method, per species:
#   1. Drop genes with zero variance across samples.
#   2. Pearson correlation between every gene pair; absolute value taken.
#   3. Rank |r| within each gene's row (base R rank, ties averaged): R[i,j].
#   4. Keep the edge i-j when its highest reciprocal rank, HRR(i,j) = max(R[i,j], R[j,i]), is at most
#      HRR_MAX (default 100). Each edge also carries the mutual rank MR = sqrt(R[i,j] * R[j,i]) and the
#      Pearson r.
#
# Run from the directory containing data/ :
#     Rscript build_species_networks.R [HRR_MAX] [outdir]
# Defaults: HRR_MAX = 100, outdir = data/networks/
# Input:  data/expression/<species>_root_cold_stress_expression.txt
# Output: data/networks/network_<species>_hrr<HRR_MAX>.tsv, and network_thresholds.tsv (per-species
#         summary, including the mutual-rank values at the 3% and 1% edge densities for comparison)
#
# Memory: the n x n correlation matrix is never held; correlations are computed in row blocks and
# reduced to ranks, and only the rank matrix persists (n^2 doubles; about 16 GB for pine).

suppressPackageStartupMessages({
  library(matrixStats)
})

args   <- commandArgs(trailingOnly = TRUE)
HRR_MAX <- if (length(args) >= 1) as.integer(args[1]) else 100L
OUTDIR  <- if (length(args) >= 2) args[2] else "data/networks"
BLOCK   <- 2000L
DENSITIES <- c(0.03, 0.01)   # the published network cut, and its neighbourhood cut

# Sparsification by HRR <= HRR_MAX rather than a fixed number of neighbours per gene: a pair is kept
# only if each gene ranks the other within its own top HRR_MAX, so highly connected genes keep many
# neighbours and peripheral genes few. Every retained edge carries both ranks, so a user can filter
# further (a lower HRR, a mutual-rank cut or an absolute Pearson r) without recomputing. The 3% density
# used for the pairwise comparisons is not used here; at whole-genome scale it implies a mean degree of
# 0.03 x (genes - 1), several hundred neighbours per gene.

SPECIES <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

read_expr <- function(species) {
  f <- sprintf("data/expression/%s_root_cold_stress_expression.txt", species)
  d <- read.delim(f, sep = "\t", header = TRUE, check.names = FALSE)
  E <- as.matrix(d[, -1, drop = FALSE])
  rownames(E) <- as.character(d[[1]])          # keep the identifiers
  E[rowSds(E) > 0, , drop = FALSE]
}

# Rank matrix of |Pearson r|, computed in row blocks so the correlation matrix is
# never materialised in full. Row i of the result holds the ranks of |r(i, .)|
# within row i, exactly as base `rank` would give them.
#
# COR_DP: |r| is rounded to this many decimal places BEFORE ranking. This is not
# cosmetic and must not be removed.
#
# Between 6.2% (birch) and 10.4% (pine) of the variable genes in these matrices share
# an identical expression profile with at least one other gene, so |r| = 1 exactly for
# those pairs and large exact-tie groups exist throughout. Meanwhile stats::cor(x) and
# stats::cor(x, y) take different BLAS paths and disagree by up to 1.33e-15. A
# difference that small decides whether two values are *exactly* equal, which moves a
# whole tie group's rank boundary -- so unrounded ranks depend on which cor() form the
# caller happened to use. Measured on a 3000-gene subset, comparing the two cor paths:
#
#     unrounded   50861 of 9000000 rank cells differ
#     12 dp          17 cells differ
#     10 dp           0
#      8 dp           0
#
# 10 dp is the finest precision that gives exact agreement, so it is used here. Without
# it, a deposited network is not reproducible on a machine with a different BLAS.
COR_DP <- 10L

rank_matrix <- function(E, block = BLOCK, dp = COR_DP) {
  n  <- nrow(E)
  tE <- t(E)
  R  <- matrix(0, nrow = n, ncol = n)
  for (from in seq(1L, n, by = block)) {
    to  <- min(from + block - 1L, n)
    C   <- round(abs(stats::cor(tE[, from:to, drop = FALSE], tE, method = "pearson")), dp)
    C[is.na(C)] <- 0
    R[from:to, ] <- matrixStats::rowRanks(C, ties.method = "average")
    rm(C)
  }
  R
}

# Genes with an identical expression profile are indistinguishable in a co-expression
# network: they sit at |r| = 1 with each other and are therefore mutual rank-1
# neighbours. They are retained rather than collapsed, because which member of such a
# group a user cares about is not ours to decide, but the count is reported so the
# resulting cliques are not mistaken for biology.
duplicate_profile_count <- function(E) {
  sum(duplicated(E)) + sum(duplicated(E, fromLast = TRUE))
}

# Walk the rank matrix in row blocks. Ranks in R are ASCENDING over |r|, so the
# strongest partner of gene i has the LARGEST value in row i. Convert to a
# co-expression rank (1 = strongest) before applying HRR.
#
# Returns the HRR-retained edges, each carrying both ranks, MR and HRR, plus the MR
# value at each requested density (for documentation, not for filtering).
hrr_edges_and_thresholds <- function(R, genes, hrr_max, densities, block = BLOCK) {
  n <- nrow(R)
  # co-expression rank: 1 for the strongest partner. R holds ascending ranks of |r|
  # including the self-comparison (|r| = 1, hence rank n), so rank_coexp = n - R.
  set.seed(1L)
  keep_frac <- min(1, 3e7 / (n * (n - 1) / 2))
  pool <- vector("list", 0L)
  si <- integer(0); ti <- integer(0)
  rs <- integer(0); rt <- integer(0); mrv <- numeric(0)
  for (from in seq(1L, n, by = block)) {
    to  <- min(from + block - 1L, n)
    nb  <- to - from + 1L
    Rb  <- R[from:to, , drop = FALSE]          # ascending ranks, rows of this block
    Rbt <- t(R[, from:to, drop = FALSE])       # the reciprocal ranks
    Cs  <- n - Rb                               # 1 = strongest, for the block's rows
    Ct  <- n - Rbt                              # reciprocal co-expression rank
    HRR <- pmax(Cs, Ct)
    # remove self-edges: gene i's own column within its row
    for (i in seq_len(nb)) HRR[i, from + i - 1L] <- NA_integer_
    sel <- which(HRR <= hrr_max, arr.ind = TRUE)
    if (nrow(sel)) {
      si  <- c(si, from + sel[, 1] - 1L)
      ti  <- c(ti, sel[, 2])
      rs  <- c(rs, Cs[sel])
      rt  <- c(rt, Ct[sel])
      mrv <- c(mrv, sqrt(Rb[sel] * Rbt[sel]))
    }
    # sample of the strictly-upper triangle of MR, for the density thresholds only
    rowid <- matrix(rep(from:to, times = n), nrow = nb)
    colid <- matrix(rep(seq_len(n), each = nb), nrow = nb)
    idx <- which(rowid < colid)
    v <- sqrt(Rb[idx] * Rbt[idx])
    if (keep_frac < 1) v <- v[stats::runif(length(v)) < keep_frac]
    pool[[length(pool) + 1L]] <- v
    rm(Rb, Rbt, Cs, Ct, HRR, sel, rowid, colid, idx, v)
  }
  allv <- unlist(pool, use.names = FALSE)
  thr  <- vapply(densities, function(dn)
                 unname(stats::quantile(allv, 1 - dn, names = FALSE)), numeric(1))
  # each undirected edge appears twice (i->j and j->i); keep one copy
  keep <- si < ti
  list(edges = data.frame(source = genes[si[keep]], target = genes[ti[keep]],
                          rank_source_to_target = rs[keep],
                          rank_target_to_source = rt[keep],
                          hrr = pmax(rs[keep], rt[keep]),
                          mr = mrv[keep], stringsAsFactors = FALSE),
       thresholds = setNames(thr, paste0("mr_at_density_", densities)),
       n_sampled = length(allv))
}

if (isTRUE(getOption("bsn.functions_only"))) {
  message("build_species_networks.R: functions loaded, main loop skipped")
} else {

summary_rows <- list()
for (s in SPECIES) {
  t0 <- Sys.time()
  E  <- read_expr(s)
  cat(sprintf("%-7s %6d variable genes, %d samples\n", s, nrow(E), ncol(E)))
  R  <- rank_matrix(E)
  out <- hrr_edges_and_thresholds(R, rownames(E), HRR_MAX, DENSITIES)
  rm(R); gc(verbose = FALSE)

  # attach the underlying Pearson r for every retained edge, so a tool can threshold
  # on a quantity users understand rather than only on mutual rank
  tE <- t(E)
  ei <- match(out$edges$source, rownames(E)); ej <- match(out$edges$target, rownames(E))
  rr <- numeric(length(ei))
  for (from in seq(1L, length(ei), by = 200000L)) {
    to <- min(from + 200000L - 1L, length(ei))
    rr[from:to] <- vapply(from:to, function(q)
      stats::cor(tE[, ei[q]], tE[, ej[q]], method = "pearson"), numeric(1))
  }
  out$edges$pearson_r <- round(rr, 4)
  out$edges$mr <- round(out$edges$mr, 2)
  out$edges <- out$edges[order(out$edges$source, out$edges$hrr), ]

  f <- file.path(OUTDIR, sprintf("network_%s_hrr%d.tsv", s, HRR_MAX))
  write.table(out$edges, f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  deg <- table(c(out$edges$source, out$edges$target))
  summary_rows[[s]] <- data.frame(
    species = s, genes = nrow(E), samples = ncol(E), hrr_max = HRR_MAX,
    cor_decimal_places = COR_DP,
    genes_with_duplicate_profile = duplicate_profile_count(E),
    edges = nrow(out$edges),
    genes_with_edges = length(deg),
    mean_degree = round(mean(as.numeric(deg)), 1),
    median_degree = stats::median(as.numeric(deg)),
    max_degree = max(as.numeric(deg)),
    mean_degree_at_density_0.03 = round(0.03 * (nrow(E) - 1)),
    mr_at_density_0.03 = round(out$thresholds[["mr_at_density_0.03"]], 1),
    mr_at_density_0.01 = round(out$thresholds[["mr_at_density_0.01"]], 1),
    min_abs_r_retained = round(min(abs(out$edges$pearson_r)), 3),
    median_abs_r_retained = round(stats::median(abs(out$edges$pearson_r)), 3),
    stringsAsFactors = FALSE)
  cat(sprintf("        %s edges, %.1f min\n", format(nrow(out$edges), big.mark = ","),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(E, out, tE, rr); gc(verbose = FALSE)
}

sm <- do.call(rbind, summary_rows)
write.table(sm, file.path(OUTDIR, "network_thresholds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(sm, row.names = FALSE)
cat("\nwrote", nrow(sm), "networks to", normalizePath(OUTDIR), "\n")

}
