#!/usr/bin/env Rscript
# stability_why_spruce.R — design and cold-response signal per genotype, to explain why spruce members fail the sample-
# stability rule most often (ledger R041). Sample IDs are <genotype><replicate>.<time point>; time points 01-05 are
# control, 6 h, 24 h, 3 d and 10 d at 5 C (spruce_cold_roots.csv).
#   design        samples, time points, replicates per time point
#   time_adjR2    per gene, adjusted R2 of a one-way model of expression on time point (adjusted, because the unadjusted R2
#                 has a null expectation of (k-1)/(n-1), which is larger with fewer samples); genes passing ComPlEx's filters
#                 (SD > 0, at least 3 non-zero values, mean > 1); reported for all such genes and for the members of the
#                 representative cliques of orthogroups that stay core or are sample-sensitive
#   replicate_r   mean Pearson correlation between replicates of the same time point, and between samples of different
#                 time points (all filtered genes)
# Inputs (working tree): data/expression/<genotype>_root_cold_stress_expression.txt,
#   <interpretation dir>/kept_lost_orthogroups.tsv (stability_kept_lost.R), data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS
# Outputs (<interpretation dir>): why_spruce_design.tsv, why_spruce_signal.tsv
# Usage (working tree): Rscript stability_why_spruce.R <interpretation dir>
suppressPackageStartupMessages(library(data.table))
D <- commandArgs(trailingOnly = TRUE)[1]
G <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
kl <- fread(file.path(D, "kept_lost_orthogroups.tsv"))
w <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS"))
des <- list(); sig <- list()
for (g in G) {
  ex <- fread(sprintf("data/expression/%s_root_cold_stress_expression.txt", g)); m <- as.matrix(ex[, -1]); rownames(m) <- ex[[1]]
  keep <- apply(m, 1, sd) > 0 & rowSums(m != 0) >= 3 & rowMeans(m) > 1; m <- m[keep, ]
  tp <- factor(sub(".*\\.", "", colnames(m))); n <- ncol(m); k <- nlevels(tp)
  des[[g]] <- data.table(genotype = g, samples = n, time_points = k, replicates_per_time_point = paste(table(tp), collapse = ","),
                         genes_after_filter = nrow(m))
  X <- model.matrix(~ tp); fit <- m %*% X %*% solve(crossprod(X)) %*% t(X)
  sst <- rowSums((m - rowMeans(m))^2); sse <- rowSums((m - fit)^2)
  adj <- 1 - (sse / (n - k)) / (sst / (n - 1))
  r <- cor(m); same <- outer(tp, tp, "==") & upper.tri(r); diff <- outer(tp, tp, "!=") & upper.tri(r)
  sets <- list(all_filtered = rownames(m),
               kept_members = intersect(w[cliqueID %in% kl[core_status == "core", clique]][[g]], rownames(m)),
               lost_members = intersect(w[cliqueID %in% kl[core_status == "sample-sensitive", clique]][[g]], rownames(m)))
  for (s in names(sets)) sig[[paste(g, s)]] <- data.table(genotype = g, gene_set = s, genes = length(sets[[s]]),
    median_time_adjR2 = median(adj[sets[[s]]]), share_adjR2_above_0.5 = mean(adj[sets[[s]]] > 0.5),
    replicate_r = mean(r[same]), between_time_r = mean(r[diff]))
}
des <- rbindlist(des); sig <- rbindlist(sig)
fwrite(des, file.path(D, "why_spruce_design.tsv"), sep = "\t"); fwrite(sig, file.path(D, "why_spruce_signal.tsv"), sep = "\t")
print(des, row.names = FALSE); print(sig, row.names = FALSE)
