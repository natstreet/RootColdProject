#!/usr/bin/env Rscript
# pseudo_family_complex.R — the ComPlEx co-expressolog test for one species pair, applied to the pseudo-families of pseudo_families.R
# (ledger R071, D097). The networks and the test are those of scripts/ComPlEx/ComPlEx.R, unchanged: expression filter (SD > 0, at least
# 3 non-zero values, mean > 1), genes restricted to those with an expressed orthologue in the pair (data/ComPlEx/orthologs-<s1>-<s2>.RData),
# absolute Pearson correlation, mutual rank sqrt(R * t(R)), neighbours = the top 3% of the upper triangle, and per pair a hypergeometric
# test in each direction of the overlap between one gene's neighbours and the orthologues (through the real orthology table) of the other
# gene's neighbours; pairs without overlap in either direction are dropped, BH in each direction, Max.p.val = the larger, co-expressolog if
# Max.p.val < 0.1. The test is vectorised over precomputed neighbour sets (same arithmetic as the loop in ComPlEx.R).
# Known-answer check: the real orthologue pairs are scored first and must reproduce the deposited comparison table (the same pairs kept, and
# Max.p.val within 0.5% relative, the precision of its stored three significant figures); the script stops otherwise.
# Pseudo pairs: in each replicate, every pseudo-family's genes of <s1> x its genes of <s2> (genes in the pair's networks only); BH over the
# replicate's pseudo pairs, as ComPlEx does over a comparison's pairs.
# Outputs (<out dir>): rate_<s1>-<s2>.tsv (pairs tested and co-expressologs, real and per replicate), coexpressologs_<s1>-<s2>.tsv.gz (the
# pseudo co-expressolog pairs: replicate, OrthoGroup, gene1, gene2), kat_<s1>-<s2>.txt
# Usage (working tree): Rscript scripts/validation/pseudo_family_complex.R <s1> <s2> <pseudo families tsv> <out dir>
suppressPackageStartupMessages({ library(data.table); library(matrixStats) })
a <- commandArgs(trailingOnly = TRUE); s1 <- a[1]; s2 <- a[2]; PF <- a[3]; OUT <- a[4]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
expr <- function(s) { x <- read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", s), check.names = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m <- m[rowSds(m) > 0 & rowCounts(m != 0) >= 3, , drop = FALSE]; m[rowMeans(m) > 1, , drop = FALSE] }
e1 <- expr(s1); e2 <- expr(s2)
ortho <- as.data.table(ld(sprintf("data/ComPlEx/orthologs-%s-%s.RData", s1, s2)))[, .(Species1, Species2, OrthoGroup)]
ortho <- ortho[Species1 %in% rownames(e1) & Species2 %in% rownames(e2)]
e1 <- e1[rownames(e1) %in% ortho$Species1, , drop = FALSE]; e2 <- e2[rownames(e2) %in% ortho$Species2, , drop = FALSE]
net <- function(m) { C <- abs(cor(t(m), method = "pearson")); R <- t(apply(C, 1, rank)); rm(C); Nt <- sqrt(R * t(R)); rm(R); diag(Nt) <- 0
  thr <- sort(Nt[upper.tri(Nt, diag = FALSE)], decreasing = TRUE); thr <- thr[round(0.03 * length(thr))]
  nb <- lapply(seq_len(nrow(Nt)), function(i) which(Nt[i, ] >= thr)); names(nb) <- rownames(Nt); nb }
t0 <- proc.time(); nb1 <- net(e1); gc(); nb2 <- net(e2); gc(); cat("networks:", (proc.time() - t0)[["elapsed"]], "s\n")
g1 <- names(nb1); g2 <- names(nb2); N1 <- length(g1); N2 <- length(g2)
# orthologues (through the real table) of each gene's neighbours, as gene indices in the other species
i1 <- match(ortho$Species1, g1); i2 <- match(ortho$Species2, g2)
o1_of_2 <- split(i1, i2); o2_of_1 <- split(i2, i1)                              # species-1 orthologues of each species-2 gene, and back
map_set <- function(nb, o) lapply(nb, function(v) unique(unlist(o[as.character(v)], use.names = FALSE)))
on2 <- map_set(nb2, o1_of_2); on1 <- map_set(nb1, o2_of_1); gc()                 # on2[[b]]: species-1 genes orthologous to b's neighbours
score <- function(a1, b2) {                                                     # a1, b2: gene indices of the pairs
  m1 <- lengths(nb1)[a1]; k1 <- lengths(on2)[b2]; x1 <- mapply(function(p, q) sum(nb1[[p]] %in% on2[[q]]), a1, b2)
  m2 <- lengths(nb2)[b2]; k2 <- lengths(on1)[a1]; x2 <- mapply(function(p, q) sum(nb2[[q]] %in% on1[[p]]), a1, b2)
  p1 <- ifelse(x1 > 1, phyper(x1 - 1, m1, N1 - m1, k1, lower.tail = FALSE), 1)
  p2 <- ifelse(x2 > 1, phyper(x2 - 1, m2, N2 - m2, k2, lower.tail = FALSE), 1)
  keep <- x1 > 0 & x2 > 0
  data.table(keep, p1, p2) }
bh_max <- function(sc) { sc <- copy(sc); sc[keep == TRUE, `:=`(q1 = p.adjust(p1, "fdr"), q2 = p.adjust(p2, "fdr"))]; sc[, maxq := pmax(q1, q2)]; sc }
# known-answer check on the real pairs
t0 <- proc.time(); rs <- bh_max(score(match(ortho$Species1, g1), match(ortho$Species2, g2)))
real <- cbind(ortho, rs)[keep == TRUE]
dep <- as.data.table(ld(sprintf("data/ComPlEx/PairwiseComparisons/comparison-%s-%s-abspearsonMR0.03-table.RData", s1, s2)))[, .(Species1, Species2, Max.p.val)]
cmp <- merge(real[, .(Species1, Species2, maxq)], dep, by = c("Species1", "Species2"), all = TRUE)
cmp[, dep_p := as.numeric(Max.p.val)]                                          # stored as text to three significant figures
mism <- cmp[is.na(maxq) | is.na(dep_p) | !(abs(maxq - dep_p) <= 5e-3 * pmax(maxq, dep_p) | pmax(maxq, dep_p) < 1e-300)]
kat <- sprintf("%s-%s: real pairs kept %d, deposited %d, Max.p.val differing by more than 0.5%% (relative) %d; %.0f s", s1, s2, nrow(real), nrow(dep),
               nrow(mism), (proc.time() - t0)[["elapsed"]])
writeLines(kat, file.path(OUT, sprintf("kat_%s-%s.txt", s1, s2))); cat(kat, "\n")
if (nrow(mism)) { fwrite(head(mism, 50), file.path(OUT, sprintf("kat_mismatch_%s-%s.tsv", s1, s2)), sep = "\t"); stop("known-answer check failed") }
# pseudo pairs
pf <- fread(PF)
rate <- list(data.table(set = "real", replicate = 0L, pairs_tested = nrow(ortho), pairs_kept = nrow(real), coexpressologs = real[maxq < 0.1, .N]))
coex <- list()
for (r in sort(unique(pf$replicate))) {
  p <- pf[replicate == r]; A <- p[genotype == s1 & gene %in% g1, .(OrthoGroup, gene1 = gene)]; B <- p[genotype == s2 & gene %in% g2, .(OrthoGroup, gene2 = gene)]
  pr <- merge(A, B, by = "OrthoGroup", allow.cartesian = TRUE)
  sc <- bh_max(score(match(pr$gene1, g1), match(pr$gene2, g2))); pr <- cbind(pr, sc)
  rate[[length(rate) + 1]] <- data.table(set = "pseudo", replicate = r, pairs_tested = nrow(pr), pairs_kept = pr[keep == TRUE, .N], coexpressologs = pr[keep == TRUE & maxq < 0.1, .N])
  coex[[length(coex) + 1]] <- pr[keep == TRUE & maxq < 0.1, .(replicate = r, OrthoGroup, gene1, gene2)]
  message(s1, "-", s2, " replicate ", r, ": ", nrow(pr), " pairs") }
fwrite(rbindlist(rate), file.path(OUT, sprintf("rate_%s-%s.tsv", s1, s2)), sep = "\t")
fwrite(rbindlist(coex), file.path(OUT, sprintf("coexpressologs_%s-%s.tsv.gz", s1, s2)), sep = "\t")
