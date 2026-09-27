#!/usr/bin/env Rscript
# conservation_sizeaware_null.R — family-size-aware null for the conserved co-expression core (companion to the
# label-permutation null in conservation_label_permutation.R).
#
# The label-permutation null replaces each comparison's co-expressolog orthogroups by a uniform random draw of the same
# size, so every orthogroup is equally likely to be drawn. Large orthogroups contain more gene pairs and are therefore
# more often co-expressologs in any one comparison. This null keeps that dependence: within each of the 15
# comparisons, an orthogroup is included with the observed inclusion rate of orthogroups of its size class in that
# comparison, and comparisons are independent. The statistic is the number of orthogroups included in all 15.
#
# Universe: the orthogroups with at least one gene in each of the five species columns of Orthogroups_20240613.tsv
# (Arabidopsis, aspen, birch, Norway spruce, Scots pine), 7,318 at the time of writing. This is not the 6,643
# "orthogroups with genes in all six genotypes" of og_summary.RData (scripts/validation/common_orthogroups.R).
# Orthogroup size is the number of genes over those five columns; size classes are bounded by EDGES.
#
# Inputs (deposit v4 or later, unpacked into ./data): data/open/co_expressologs.tsv, data/*/Orthogroups_20240613.tsv
# Outputs: results/conservation_sizeaware_null.csv (null_core_size), results/conservation_sizeaware_summary.tsv
# Usage:   Rscript scripts/validation/conservation_sizeaware_null.R

P <- 0.1; N_SIM <- 5000; set.seed(1)
EDGES <- c(2, 6, 10, 15, 20, 30, 50, 75, 100, 150, 250, 500, 1e6)
SP <- c("Arabidopsis_thaliana", "Populus_tremula", "Betula_pendula", "Picea_abies", "Pinus_sylvestris")

co <- read.delim("data/open/co_expressologs.tsv", stringsAsFactors = FALSE)
co <- co[co$MaxpVal < P, ]
cmp <- paste(pmin(co$Species1, co$Species2), pmax(co$Species1, co$Species2))
per <- lapply(split(co$OrthologGroup, cmp), unique)
stopifnot(length(per) == 15)

og <- read.delim(Sys.glob("data/*/Orthogroups_20240613.tsv")[1], check.names = FALSE, stringsAsFactors = FALSE,
                 colClasses = "character")
n <- sapply(SP, function(s) { x <- og[[s]]; x[is.na(x)] <- ""; ifelse(nzchar(x), lengths(strsplit(x, ",")), 0L) })
keep <- rowSums(n > 0) == length(SP)
U <- og$Ortholog_Group[keep]; size <- rowSums(n[keep, , drop = FALSE])

M <- t(vapply(per, function(s) U %in% s, logical(length(U))))       # comparisons x universe
obs <- sum(colSums(M) == nrow(M))
b <- findInterval(size, EDGES, left.open = FALSE)                  # size class, as numpy.digitize(size, EDGES)
R <- M * 0
for (k in unique(b)) R[, b == k] <- rowMeans(M[, b == k, drop = FALSE])
expected <- sum(apply(R, 2, prod))

null <- vapply(seq_len(N_SIM), function(i) sum(colSums(matrix(runif(length(R)), nrow(R)) < R) == nrow(R)), integer(1))

dir.create("results", showWarnings = FALSE)
write.csv(data.frame(null_core_size = null), "results/conservation_sizeaware_null.csv", row.names = FALSE)
out <- data.frame(quantity = c("observed", "universe_orthogroups_with_genes_in_all_five_species", "expected", "null_mean", "null_max", "fold", "empirical_P", "n_sim"),
                  value = c(obs, length(U), expected, mean(null), max(null), obs / expected,
                            (sum(null >= obs) + 1) / (N_SIM + 1), N_SIM))
write.table(out, "results/conservation_sizeaware_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
