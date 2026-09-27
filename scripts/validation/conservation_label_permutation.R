#!/usr/bin/env Rscript
# conservation_label_permutation.R — permutation test for the conserved co-expression core
# (Fig. S6a).
#
# Observed statistic: the number of orthogroups that carry a co-expressolog (MaxpVal < 0.1) in all
# 15 pairwise species comparisons. Null: in each permutation, each comparison's set of co-expressolog
# orthogroups is replaced by a random draw of the same size from the eligible universe (orthogroups
# with genes in all five species, so present in all six genotypes), and the all-15 count is taken.
#
# Inputs (run from a directory containing the unpacked figshare deposit `data/`):
#   data/ComPlEx/co_expressologs.RData        co-expressolog pairs with MaxpVal
#   data/ComPlEx/Orthogroups_20240613.tsv     orthogroup membership (one column per species)
# Output: conservation_null_null.csv (column null_core_size), conservation_null_obs.txt
# Usage:  Rscript scripts/validation/conservation_label_permutation.R

ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
P <- 0.1; N_PERM <- 5000; set.seed(1)

co <- ld("data/ComPlEx/co_expressologs.RData")
co <- co[co$MaxpVal < P, ]
per <- tapply(co$OrthologGroup, paste(co$Species1, co$Species2), function(x) unique(as.character(x)))
stopifnot(length(per) == 15)
obs <- sum(table(unlist(per)) == 15)

og <- read.delim(Sys.glob("data/*/Orthogroups_20240613.tsv")[1], check.names = FALSE, stringsAsFactors = FALSE)
spcols <- c("Arabidopsis_thaliana", "Populus_tremula", "Betula_pendula", "Picea_abies", "Pinus_sylvestris")
universe <- og$Ortholog_Group[rowSums(sapply(spcols, function(s) nzchar(og[[s]]) & !is.na(og[[s]]))) == 5]
stopifnot(all(unique(unlist(per)) %in% og$Ortholog_Group))
sizes <- vapply(per, function(x) sum(x %in% universe), integer(1))

null <- vapply(seq_len(N_PERM), function(i) {
  draws <- unlist(lapply(sizes, function(k) sample(universe, k)))
  sum(tabulate(match(draws, universe), nbins = length(universe)) == 15)
}, integer(1))

write.csv(data.frame(null_core_size = null), "conservation_null_null.csv", row.names = FALSE)
writeLines(as.character(obs), "conservation_null_obs.txt")
cat(sprintf("observed = %d | universe = %d orthogroups | null mean = %.3f, max = %d | P = %.2g (%d permutations)\n",
            obs, length(universe), mean(null), max(null), (sum(null >= obs) + 1) / (N_PERM + 1), N_PERM))
