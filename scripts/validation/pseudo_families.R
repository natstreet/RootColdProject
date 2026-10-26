#!/usr/bin/env Rscript
# pseudo_families.R — unrelated-partner null for the conserved core (ledger R071, D097). For every orthogroup eligible for the core
# (results/orthogroup_size_rule.tsv, 6,889; D080) and every genotype, a pseudo-family takes the genes of a DIFFERENT orthogroup, drawn at
# random among eligible orthogroups with the same number of genes passing the expression filter in that genotype (the ComPlEx filter:
# SD > 0, at least 3 non-zero values, mean > 1); if none has the same number, among those in the same decile of that number. The donor
# is never the real orthogroup and differs between the six genotypes of a pseudo-family. Col-0 and Ost-0 draw separately.
# Output: <out tsv> with replicate, OrthoGroup (the real orthogroup the pseudo-family stands for), genotype, donor, gene, and how the
# donor was matched (exact or decile).
# Usage (working tree): Rscript scripts/validation/pseudo_families.R <replicates> <seed> <out tsv>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); NR <- as.integer(a[1]); set.seed(as.integer(a[2])); OUT <- a[3]
GENO <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
COL <- c(col0 = "Arabidopsis_thaliana", ost0 = "Arabidopsis_thaliana", aspen = "Populus_tremula", birch = "Betula_pendula",
         spruce = "Picea_abies", pine = "Pinus_sylvestris")
elig <- fread("results/orthogroup_size_rule.tsv")$OrthoGroup
og <- read.delim("data/ComPlEx/Orthogroups_20240613.tsv", check.names = FALSE, colClasses = "character"); names(og)[1] <- "OrthoGroup"
og <- og[og$OrthoGroup %in% elig, ]; stopifnot(nrow(og) == length(elig))
pass <- lapply(setNames(GENO, GENO), function(s) { x <- read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", s))
  m <- as.matrix(x[, -1]); x[[1]][matrixStats::rowSds(m) > 0 & rowSums(m != 0) >= 3 & rowMeans(m) > 1] })
G <- rbindlist(lapply(GENO, function(s) { v <- strsplit(og[[COL[[s]]]], ",", fixed = TRUE)
  d <- data.table(OrthoGroup = rep(og$OrthoGroup, lengths(v)), gene = trimws(unlist(v)))[gene %in% pass[[s]]]; d[, genotype := s] }))
n <- dcast(G[, .N, by = .(OrthoGroup, genotype)], OrthoGroup ~ genotype, value.var = "N", fill = 0)
n <- merge(data.table(OrthoGroup = elig), n, by = "OrthoGroup", all.x = TRUE); for (s in GENO) n[is.na(get(s)), (s) := 0L]
genes_of <- split(G$gene, paste(G$OrthoGroup, G$genotype))
res <- vector("list", NR)
for (r in seq_len(NR)) {
  out <- vector("list", length(GENO))
  used <- matrix("", nrow(n), length(GENO), dimnames = list(n$OrthoGroup, GENO))
  for (j in seq_along(GENO)) { s <- GENO[j]; cnt <- n[[s]]; dec <- cut(cnt, unique(quantile(cnt, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)
    pool_exact <- split(n$OrthoGroup, cnt); pool_dec <- split(n$OrthoGroup, dec)
    don <- character(nrow(n)); how <- character(nrow(n))
    for (i in seq_len(nrow(n))) {
      excl <- c(n$OrthoGroup[i], used[i, seq_len(j - 1)])
      pool <- pool_exact[[as.character(cnt[i])]]; how[i] <- "exact"
      if (sum(!pool %in% excl) == 0) { pool <- pool_dec[[as.character(dec[i])]]; how[i] <- "decile" }
      stopifnot(sum(!pool %in% excl) > 0)
      repeat { d <- pool[sample.int(length(pool), 1)]; if (!d %in% excl) break }   # uniform over the pool minus the excluded
      don[i] <- d }
    used[, j] <- don
    gl <- genes_of[paste(don, s)]
    out[[j]] <- data.table(replicate = r, OrthoGroup = rep(n$OrthoGroup, lengths(gl)), genotype = s, donor = rep(don, lengths(gl)),
                           gene = unlist(gl, use.names = FALSE), matched = rep(how, lengths(gl))) }
  res[[r]] <- rbindlist(out); message("replicate ", r) }
res <- rbindlist(res)
stopifnot(res[, all(OrthoGroup != donor)])
chk <- unique(res[, .(replicate, OrthoGroup, genotype, donor)])[, .(d = uniqueN(donor), n = .N), by = .(replicate, OrthoGroup)]
stopifnot(chk[, all(d == n)])                                                   # a different donor for every genotype of a pseudo-family
fwrite(res, OUT, sep = "\t")
print(unique(res[, .(replicate, OrthoGroup, genotype, matched)])[, .N, by = matched])
