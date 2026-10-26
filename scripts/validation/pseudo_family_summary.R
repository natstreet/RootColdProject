#!/usr/bin/env Rscript
# pseudo_family_summary.R — the unrelated-partner null (ledger R071, D097) summarised from the 15 outputs of pseudo_family_complex.R.
# Per replicate: (a) pseudo-families with a co-expressolog pair in all 15 comparisons; (b) pseudo-families with a complete six-genotype
# clique (one gene per genotype, all 15 pairs co-expressologs), and (b') those the size rule counts (the real orthogroup's fitted
# probability < 0.5 in results/orthogroup_size_rule.tsv, i.e. the pseudo-family takes the size class of the orthogroup it stands for).
# Observed: the same counts from the real comparison tables (264 / 137 / 130). Also by size bin (deciles of the size rule's expressed
# genes) and the co-expressolog rate per comparison, real vs pseudo.
# The clique search is has_clique() of pair_level_null.R (pruning plus forward checking; known-answer test there and here).
# Outputs (<out dir>): null_summary.tsv, null_by_bin.tsv, rate_by_comparison.tsv, per_replicate.tsv
# Usage (working tree): Rscript scripts/validation/pseudo_family_summary.R <dir with rate_*.tsv and coexpressologs_*.tsv.gz> <out dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); IN <- a[1]; OUT <- a[2]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
GENO <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
src <- readLines("scripts/validation/pair_level_null.R")                       # reuse has_clique() unchanged
eval(parse(text = src[grep("^has_clique <- function", src):(grep("^# known-answer test", src) - 1)]))
toy <- CJ(i = 1:6, j = 1:6)[i < j][, .(g1 = paste0(GENO[i], ":x"), g2 = paste0(GENO[j], ":x"))]
stopifnot(has_clique(toy), !has_clique(toy[-3]))
files <- list.files(IN, "^coexpressologs_.*[.]tsv[.]gz$", full.names = TRUE); stopifnot(length(files) == 15)
pairs <- rbindlist(lapply(files, function(f) { p <- strsplit(sub("^coexpressologs_(.*)[.]tsv[.]gz$", "\\1", basename(f)), "-")[[1]]
  x <- fread(f); if (!nrow(x)) return(NULL); x[, `:=`(g1 = paste0(p[1], ":", gene1), g2 = paste0(p[2], ":", gene2), cmp = paste(p, collapse = "-"))] }))
rates <- rbindlist(lapply(list.files(IN, "^rate_.*[.]tsv$", full.names = TRUE), function(f) fread(f)[, cmp := sub("^rate_(.*)[.]tsv$", "\\1", basename(f))]))
reps <- sort(unique(rates[set == "pseudo", replicate])); stopifnot(rates[set == "pseudo", .N, by = cmp][, all(N == length(reps))])
sr <- fread("results/orthogroup_size_rule.tsv")
bins <- setNames(cut(sr$expressed_genes, unique(quantile(sr$expressed_genes, 0:10 / 10)), include.lowest = TRUE, labels = FALSE), sr$OrthoGroup)
rule <- setNames(sr$p_core_given_size < 0.5, sr$OrthoGroup)
count <- function(pp) {                                                         # pp: co-expressolog pairs with OrthoGroup, cmp, g1, g2
  a15 <- pp[, .(n = uniqueN(cmp)), by = OrthoGroup][n == 15, OrthoGroup]
  cl <- if (length(a15)) a15[vapply(split(pp[OrthoGroup %in% a15, .(g1, g2, OrthoGroup)], by = "OrthoGroup")[a15], has_clique, TRUE)] else character(0)
  cr <- cl[rule[cl] %in% TRUE]
  list(all15 = a15, clique = cl, clique_rule = cr) }
# observed, from the real comparison tables (the families themselves)
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
obs_pairs <- rbindlist(lapply(list.files("data/ComPlEx/PairwiseComparisons", "table[.]RData$", full.names = TRUE), function(f) {
  p <- strsplit(sub("^comparison-(.*)-abspearsonMR0[.]03-table[.]RData$", "\\1", basename(f)), "-")[[1]]
  x <- as.data.table(ld(f))[as.numeric(Max.p.val) < 0.1]
  x[, .(OrthoGroup, cmp = paste(p, collapse = "-"), g1 = paste0(p[1], ":", Species1), g2 = paste0(p[2], ":", Species2))] }))
obs <- count(obs_pairs); cat("observed:", lengths(obs), "\n")
per <- rbindlist(lapply(reps, function(r) { cc <- count(pairs[replicate == r])
  rbindlist(lapply(names(cc), function(k) data.table(replicate = r, statistic = k, bin = c(seq_len(max(bins, na.rm = TRUE)), NA),
    n = as.integer(table(factor(bins[cc[[k]]], seq_len(max(bins, na.rm = TRUE))), useNA = "always"))))) }))
fwrite(per, file.path(OUT, "per_replicate.tsv"), sep = "\t")
tot <- per[, .(n = sum(n)), by = .(replicate, statistic)]
s <- tot[, .(observed = length(obs[[statistic[1]]]), replicates = .N, null_mean = round(mean(n), 2), null_min = min(n), null_max = max(n),
             null_lo95 = quantile(n, 0.025), null_hi95 = quantile(n, 0.975)), by = statistic]
s[, P := mapply(function(k, o) (sum(tot[statistic == k, n] >= o) + 1) / (sum(tot$statistic == k) + 1), statistic, observed)]
s[, fold := round(observed / pmax(null_mean, 1 / replicates), 1)]
fwrite(s, file.path(OUT, "null_summary.tsv"), sep = "\t"); print(s)
ob <- rbindlist(lapply(names(obs), function(k) data.table(statistic = k, bin = c(seq_len(max(bins, na.rm = TRUE)), NA),
  observed = as.integer(table(factor(bins[obs[[k]]], seq_len(max(bins, na.rm = TRUE))), useNA = "always")))))
bb <- merge(ob, per[, .(null_mean = round(mean(n), 2), null_max = max(n)), by = .(statistic, bin)], by = c("statistic", "bin"), all = TRUE)
edges <- sr[, .(min_genes = min(expressed_genes), max_genes = max(expressed_genes), orthogroups = .N), by = .(bin = bins[OrthoGroup])]
bb <- merge(bb, edges, by = "bin", all.x = TRUE); fwrite(bb, file.path(OUT, "null_by_bin.tsv"), sep = "\t"); print(bb)
rr <- rates[, .(real_rate = round(coexpressologs[set == "real"] / pairs_kept[set == "real"], 4),
                pseudo_rate_mean = round(mean(coexpressologs[set == "pseudo"] / pairs_kept[set == "pseudo"]), 4),
                pseudo_pairs_tested_mean = round(mean(pairs_tested[set == "pseudo"])), real_pairs_tested = pairs_tested[set == "real"]), by = cmp]
fwrite(rr, file.path(OUT, "rate_by_comparison.tsv"), sep = "\t"); print(rr)
