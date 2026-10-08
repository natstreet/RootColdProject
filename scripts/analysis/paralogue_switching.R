#!/usr/bin/env Rscript
# paralogue_switching.R — do Col-0 and Ost-0 keep conserved co-expression through different paralogues? (DECISIONS.md D062
# item 2, D063; design R040 as amended in A035; on the final core, D066.) Descriptive, no test.
# For each core orthogroup whose representative clique (lowest MeanCliqueSum) has different Col-0 and Ost-0 members:
#   (a) whether a complete clique of the orthogroup has the same locus as its Col-0 and Ost-0 member
#   (b) if so, the MeanCliqueSum difference between the best such clique and the representative, and its rank (percentile)
#       among all pairwise MeanCliqueSum differences between the orthogroup's complete cliques ("near-equivalence" scale)
#   category  representative choice ((a) yes) | candidate genuine switch ((a) no; reviewed for quantification artefact by
#             paralogue_switching_artefact.R: CDS identity >= 95% -> possible quantification artefact)
# For every case: whether the representative's Col-0/Ost-0 locus pair is a col0-ost0 co-expressolog (Max.p.val < 0.1), and
#   the two loci's mean expression in each accession.
# Inputs (working tree): data/ComPlEx/cliques/conserved_genes_COMPLETE_{filtered,wider}.RDS, data/ComPlEx/PairwiseComparisons/
#   comparison-col0-ost0-abspearsonMR0.03-table.RData, data/expression/{col0,ost0}_root_cold_stress_expression.txt
# Output: results/paralogue_switching.tsv. Usage (working tree): Rscript scripts/analysis/paralogue_switching.R
suppressPackageStartupMessages(library(data.table))
ld <- function(f) { e <- new.env(); load(f, envir = e); as.data.table(get(ls(e)[1], e)) }
rep <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
w <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS"))[OrthoGroup %in% rep$OrthoGroup]
co <- ld("data/ComPlEx/PairwiseComparisons/comparison-col0-ost0-abspearsonMR0.03-table.RData")[, .(Species1, Species2, p = as.numeric(Max.p.val))]
mexp <- function(g) { x <- fread(sprintf("data/expression/%s_root_cold_stress_expression.txt", g)); setNames(rowMeans(as.matrix(x[, -1])), x[[1]]) }
ec <- mexp("col0"); eo <- mexp("ost0")
sw <- rep[col0 != ost0]
res <- rbindlist(lapply(seq_len(nrow(sw)), function(i) {
  r <- sw[i]; cl <- w[OrthoGroup == r$OrthoGroup]
  same <- cl[col0 == ost0]; best <- if (nrow(same)) same[which.min(MeanCliqueSum)] else NULL
  d <- if (!is.null(best)) best$MeanCliqueSum - r$MeanCliqueSum else NA_real_
  pd <- if (nrow(cl) > 1) as.vector(dist(cl$MeanCliqueSum)) else numeric(0)
  pr <- if (!is.na(d) && length(pd)) mean(pd <= abs(d)) else NA_real_
  pc <- co[(Species1 == r$col0 & Species2 == r$ost0) | (Species1 == r$ost0 & Species2 == r$col0), p]
  data.table(OrthoGroup = r$OrthoGroup, representative = r$cliqueID, col0_member = r$col0, ost0_member = r$ost0,
             complete_cliques = nrow(cl), same_locus_cliques = nrow(same), best_same_locus_clique = if (is.null(best)) "" else best$cliqueID,
             same_locus = if (is.null(best)) "" else best$col0, meancliquesum_difference = signif(d, 3),
             difference_percentile = round(pr, 3), category = if (is.null(best)) "candidate genuine switch" else "representative choice",
             pair_coexpressolog_p = if (length(pc)) signif(min(pc), 3) else NA_real_,
             col0_member_mean_in_col0 = round(ec[r$col0], 2), col0_member_mean_in_ost0 = round(eo[r$col0], 2),
             ost0_member_mean_in_col0 = round(ec[r$ost0], 2), ost0_member_mean_in_ost0 = round(eo[r$ost0], 2)) }))
fwrite(res, "results/paralogue_switching.tsv", sep = "\t", na = "")
cat("core representatives with different Col-0 and Ost-0 members:", nrow(res), "of", nrow(rep), "\n"); print(res[, .N, by = category])
print(res[category == "representative choice", .(median_difference_percentile = median(difference_percentile, na.rm = TRUE))])
