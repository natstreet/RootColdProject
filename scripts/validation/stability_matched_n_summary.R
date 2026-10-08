#!/usr/bin/env Rscript
# stability_matched_n_summary.R — summary of stability_matched_n.R: per genotype, design and gene set, the median over query
# genes of the worst-case (min_jaccard) and mean leave-one-out neighbourhood Jaccard; matched15 values are the median over
# the 10 draws of each draw's median. Output: matched_n_summary.tsv. Usage: Rscript stability_matched_n_summary.R <dir with matched_*.tsv.gz> <out file>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
x <- rbindlist(lapply(list.files(a[1], "^matched_.*\\.tsv\\.gz$", full.names = TRUE), function(f) fread(cmd = paste("gzip -dc", shQuote(f)))))
d <- x[, .(min_j = median(min_jaccard), mean_j = median(mean_jaccard), genes = .N), by = .(genotype, design, draw, samples, set)]
s <- d[, .(samples = samples[1], draws = .N, genes = genes[1], median_worst_jaccard = round(median(min_j), 3), median_mean_jaccard = round(median(mean_j), 3),
           range_worst_over_draws = paste(round(range(min_j), 3), collapse = "-")), by = .(genotype, design, set)]
setorder(s, set, design, genotype); fwrite(s, a[2], sep = "\t"); print(s, row.names = FALSE)
