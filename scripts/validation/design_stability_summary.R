#!/usr/bin/env Rscript
# design_stability_summary.R — summary of design_stability.R: per dataset and design (r replicates per point, k points),
# the median over draws of the per-draw median worst-case leave-one-out neighbourhood Jaccard, with n; and, to separate
# replicates from points, a linear model per data type (root, wood) of the draw values on log2 r and log2 k with a dataset
# term (worst-case and mean measures); and designs of equal n with different allocations.
# Outputs (<dir>): design_summary.tsv, design_model.tsv, design_equal_n.tsv. Usage: Rscript design_stability_summary.R <dir with *.tsv>
suppressPackageStartupMessages(library(data.table))
D <- commandArgs(trailingOnly = TRUE)[1]
x <- rbindlist(lapply(list.files(D, "^(root|wood)_.*\\.tsv$", full.names = TRUE), fread))
x[, type := sub("_.*", "", dataset)]
s <- x[, .(n = as.numeric(median(n)), draws = .N, worst_j = round(median(median_worst_jaccard), 3), mean_j = round(median(median_mean_jaccard), 3)), by = .(type, dataset, r, k)]
setorder(s, type, dataset, r, k); fwrite(s, file.path(D, "design_summary.tsv"), sep = "\t")
mods <- rbindlist(lapply(c("median_worst_jaccard", "median_mean_jaccard"), function(y) rbindlist(lapply(c("root", "wood"), function(t) {
  f <- lm(as.formula(paste(y, "~ log2(r) + log2(k) + dataset")), data = x[type == t]); co <- summary(f)$coefficients
  data.table(measure = y, type = t, term = rownames(co), estimate = co[, 1], se = co[, 2], p = co[, 4], r2 = summary(f)$r.squared) }))))
# same total n, different allocation: per dataset, median over draws at each n with both r and k varying
eq <- s[, .(designs = paste0("r", r, "k", k, "=", worst_j, "/", mean_j, collapse = "  ")), by = .(type, dataset, n)][grepl("  ", designs)]
fwrite(eq, file.path(D, "design_equal_n.tsv"), sep = "\t")
fwrite(mods[grepl("log2", term)], file.path(D, "design_model.tsv"), sep = "\t")
print(dcast(s, dataset ~ paste0("r", r, "k", k), value.var = "worst_j")); print(mods[grepl("log2", term)])
