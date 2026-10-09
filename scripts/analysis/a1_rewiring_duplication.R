#!/usr/bin/env Rscript
# a1_rewiring_duplication.R — A1 (ledger N159; criteria fixed before SCI's run; R045 asks for a committed producer): are opposed cold
# responses between lineages more frequent in paralogue-carried than in orthologue-carried core orthogroups?
#   opposed      direction "Arabidopsis vs trees" or "angiosperms vs conifers"; not opposed: "same" or "mixed"; no direction (too
#                little DE) excluded (tables/core_response_patterns.tsv, core_response_patterns.R)
#   class        orthology_class (results/orthology_class.tsv, orthology_class.R; 137-tree N1)
#   primary      Fisher's exact test, paralogue- vs orthologue-carried x opposed
#   adjusted     logistic regression opposed ~ paralogue-carried + cluster (factor) + log10(expressed genes)
#   supported    primary OR > 1 with P < 0.05 AND adjusted paralogue coefficient > 0 with P < 0.05
# Output: results/a1_rewiring_duplication.tsv. Usage (working tree): Rscript scripts/analysis/a1_rewiring_duplication.R
suppressPackageStartupMessages(library(data.table))
rp <- fread("tables/core_response_patterns.tsv"); oc <- fread("results/orthology_class.tsv")
x <- merge(rp[, .(OrthoGroup, cluster, direction)], oc[, .(OrthoGroup, orthology_class, expressed_genes)], by = "OrthoGroup")
stopifnot(nrow(x) == nrow(rp))
x <- x[direction %in% c("Arabidopsis vs trees", "angiosperms vs conifers", "same", "mixed")]
x[, opposed := direction %in% c("Arabidopsis vs trees", "angiosperms vs conifers")][, paralogue := orthology_class == "paralogue-carried"]
t <- table(factor(x$paralogue, c(TRUE, FALSE)), factor(x$opposed, c(TRUE, FALSE))); f <- fisher.test(t)
g <- glm(opposed ~ paralogue + factor(cluster) + log10(expressed_genes), family = binomial, data = x)
co <- summary(g)$coefficients["paralogueTRUE", ]
res <- data.table(orthogroups_with_direction = nrow(x),
  paralogue_opposed = t["TRUE", "TRUE"], paralogue_total = sum(t["TRUE", ]), orthologue_opposed = t["FALSE", "TRUE"], orthologue_total = sum(t["FALSE", ]),
  fisher_or = unname(f$estimate), fisher_ci_low = f$conf.int[1], fisher_ci_high = f$conf.int[2], fisher_p = f$p.value,
  adjusted_coef = co[["Estimate"]], adjusted_p = co[["Pr(>|z|)"]],
  cluster1_opposed = x[cluster == 1, sum(opposed)], cluster1_with_direction = x[cluster == 1, .N])
res[, verdict := ifelse(fisher_or > 1 & fisher_p < 0.05 & adjusted_coef > 0 & adjusted_p < 0.05, "supported", "not supported")]
fwrite(res, "results/a1_rewiring_duplication.tsv", sep = "\t"); print(t(res))
