#!/usr/bin/env Rscript
# network_link_summary.R — per genotype, the positive and negative networks of the resource: number of links, the
# weakest and the median absolute jackknifed correlation among the links (the correlation the links start at), and the
# share of links that are negative.
# Usage: Rscript network_link_summary.R <networks dir> <out tsv> <name>=<file stem> ...
#   e.g. col0=network_col0  (reads <stem>_hrr100.tsv and <stem>_neg_hrr100.tsv)
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); D <- a[1]; OUT <- a[2]
rows <- rbindlist(lapply(strsplit(a[-(1:2)], "=", fixed = TRUE), function(x) {
  p <- fread(file.path(D, paste0(x[2], "_hrr100.tsv")), select = "jackknife_r")$jackknife_r
  n <- fread(file.path(D, paste0(x[2], "_neg_hrr100.tsv")), select = "jackknife_r")$jackknife_r
  stopifnot(all(p > 0), all(n < 0))
  data.table(genotype = x[1], positive_links = length(p), positive_weakest_r = round(min(p), 3),
             positive_median_r = round(median(p), 3), negative_links = length(n),
             negative_weakest_r = round(max(n), 3), negative_median_r = round(median(n), 3),
             negative_share = round(length(n) / (length(p) + length(n)), 3))
}))
fwrite(rows, OUT, sep = "\t"); print(rows)
