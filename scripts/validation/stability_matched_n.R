#!/usr/bin/env Rscript
# stability_matched_n.R — is spruce's sample sensitivity explained by its design (3 replicates per time point, 15
# samples) or by its data? (ledger R041.) Within each genotype, the stability of co-expression neighbourhoods under
# leave-one-sample-out, with the angiosperms also subsampled to spruce's design.
#   neighbourhood  for a query gene, its top 3% of filtered genes by |Pearson r| (a per-gene proxy for ComPlEx's network
#                  of 3% density on mutual ranks; ComPlEx itself is not re-run)
#   stability      Jaccard similarity between the neighbourhood from the full sample set and from the set minus one
#                  sample; per query gene, the minimum over left-out samples (the worst case, as in D052) and the mean
#   query genes    members of that genotype in the representative cliques of the 133 core orthogroups (kept and lost
#                  marked) and 300 random filtered genes (seed 1)
#   designs        full: all samples. matched15: for 5 x 5 genotypes, 3 replicates drawn per time point (15 samples, as
#                  spruce), 10 draws (seeds 1-10). Spruce and pine are run as they are.
# Filters as ComPlEx (SD > 0, at least 3 non-zero values, mean > 1), applied to each design's sample set; the leave-one-out
# runs keep that gene universe (a gene constant without one sample gets zero correlations).
# Inputs (working dir): expression/<genotype>_root_cold_stress_expression.txt, kept_lost_orthogroups.tsv, wider.RDS
# Output: matched_n_stability.tsv.gz (genotype, design, draw, gene, set, min_jaccard, mean_jaccard)
# Usage: Rscript stability_matched_n.R <expression dir> <kept_lost_orthogroups.tsv> <wider.RDS> <out file> [genotype]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); EXD <- a[1]; KL <- a[2]; WID <- a[3]; OUT <- a[4]
G <- if (length(a) > 4) a[5] else c("col0", "ost0", "aspen", "birch", "spruce", "pine")
DENS <- 0.03; NDRAW <- 10; NRAND <- 300
kl <- fread(KL); w <- as.data.table(readRDS(WID))
filt <- function(m) m[apply(m, 1, sd) > 0 & rowSums(m != 0) >= 3 & rowMeans(m) > 1, , drop = FALSE]
topset <- function(m, q) {   # rows q of m: indices of the top DENS of genes by |r| (self excluded)
  z <- t(scale(t(m))); z[!is.finite(z)] <- 0; r <- abs(tcrossprod(z[q, , drop = FALSE], z)) / (ncol(m) - 1)
  k <- ceiling(DENS * (nrow(m) - 1)); qi <- match(q, rownames(m))
  lapply(seq_along(q), function(i) { v <- r[i, ]; v[qi[i]] <- -Inf; rownames(m)[order(v, decreasing = TRUE)[seq_len(k)]] }) }
jac <- function(x, y) length(intersect(x, y)) / length(union(x, y))
run <- function(m, q, setlab) {
  m <- filt(m); q <- intersect(q, rownames(m)); full <- topset(m, q)
  J <- sapply(seq_len(ncol(m)), function(s) { sub <- topset(m[, -s], q); mapply(jac, full, sub) })   # same gene universe as the full set
  data.table(gene = q, set = setlab[q], min_jaccard = apply(J, 1, min), mean_jaccard = rowMeans(J)) }
res <- list()
for (g in G) {
  ex <- fread(file.path(EXD, sprintf("%s_root_cold_stress_expression.txt", g))); m <- as.matrix(ex[, -1]); rownames(m) <- ex[[1]]
  kept <- w[cliqueID %in% kl[core_status == "core", clique]][[g]]; lost <- w[cliqueID %in% kl[core_status == "sample-sensitive", clique]][[g]]
  set.seed(1); rnd <- sample(setdiff(rownames(filt(m)), c(kept, lost)), NRAND)
  setlab <- c(setNames(rep("random", length(rnd)), rnd), setNames(rep("lost", length(lost)), lost), setNames(rep("kept", length(kept)), kept))
  setlab <- setlab[!duplicated(names(setlab), fromLast = TRUE)]; q <- names(setlab)
  res[[paste(g, "full")]] <- run(m, q, setlab)[, `:=`(genotype = g, design = "full", draw = 0L, samples = ncol(m))]
  tp <- sub(".*\\.", "", colnames(m))
  if (all(table(tp) == 5)) for (d in seq_len(NDRAW)) { set.seed(d)
    cols <- unlist(lapply(split(seq_along(tp), tp), function(ix) sample(ix, 3)))
    res[[paste(g, d)]] <- run(m[, cols], q, setlab)[, `:=`(genotype = g, design = "matched15", draw = d, samples = length(cols))] }
  cat(g, "done\n")
}
fwrite(rbindlist(res), OUT, sep = "\t")
