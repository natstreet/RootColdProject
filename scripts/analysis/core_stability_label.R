#!/usr/bin/env Rscript
# core_stability_label.R — the sample-stability label of each conserved-core orthogroup (DECISIONS.md D066): the core is the
# 133 (D043); the D052/D059 test (scripts/validation/loo_stability_score.R) labels, it does not filter.
#   core_status            stable (at least one complete clique has all six members passing) or sample-sensitive
#   complete_cliques, passing_cliques, share_passing (passing / complete)
#   representative         the core's representative clique (lowest MeanCliqueSum; unchanged by the label)
#   representative_passes  whether that clique passes; representative_failing_members (genotypes that fail)
#   stable_representative  the lowest-MeanCliqueSum passing clique (extra column; empty when none passes)
# Input (working tree): results/loo_stability/loo_orthogroups.tsv (the scorer's output), data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS
# Output: results/core_stability_label.tsv. Usage (working tree): Rscript scripts/analysis/core_stability_label.R
suppressPackageStartupMessages(library(data.table))
og <- fread("results/loo_stability/loo_orthogroups.tsv")
core <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
if (!setequal(og$OrthoGroup, core$OrthoGroup)) stop("the stability scoring does not cover exactly the core orthogroups")
x <- og[, .(OrthoGroup, core_status = ifelse(core_status == "core", "stable", "sample-sensitive"), complete_cliques, passing_cliques,
            share_passing = round(passing_cliques / complete_cliques, 4), representative = core$cliqueID[match(OrthoGroup, core$OrthoGroup)],
            representative_passes = old_representative_passes, representative_failing_members = old_representative_failing,
            stable_representative = representative)]
fwrite(x[order(OrthoGroup)], "results/core_stability_label.tsv", sep = "\t")
print(x[, .N, by = core_status])
