#!/usr/bin/env Rscript
# compare_sign_variants.R — one species pair's ComPlEx comparison table built three ways (ComPlEx.R with cor_sign abs,
# the published method; pos, positive correlation only; neg, negative correlation only): how many orthologue pairs are
# co-expressologs (Max.p.val < 0.1, the threshold of the conserved-core cliques) under each, and how the sets overlap.
# A control checks that the abs table equals the reference table (the one the analysis used).
# Usage: Rscript compare_sign_variants.R <abs table> <pos table> <neg table> <reference table> <out tsv>
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
ld <- function(f) { e <- new.env(); load(f, envir = e); as.data.table(get(ls(e)[1], e)) }
P <- 0.1
key <- function(t) paste(t$Species1, t$Species2)
tabs <- lapply(setNames(a[1:4], c("abs", "pos", "neg", "reference")), ld)
ctl <- all.equal(as.data.frame(tabs$abs), as.data.frame(tabs$reference))
cat("control: abs table equals the reference table:", isTRUE(ctl), "\n")
sig <- lapply(tabs[c("abs", "pos", "neg")], function(t) key(t)[as.numeric(t$Max.p.val) < P])
tested <- lapply(tabs[c("abs", "pos", "neg")], function(t) nrow(t))
pn <- union(sig$pos, sig$neg)
res <- data.table(
  quantity = c("pairs in table (>= 1 shared neighbour each way)", "co-expressologs (Max.p.val < 0.1)",
               "co-expressologs shared with abs", "abs co-expressologs found by pos or neg",
               "found by pos or neg but not abs", "shared by pos and neg"),
  abs = c(tested$abs, length(sig$abs), length(sig$abs), length(intersect(sig$abs, pn)), NA, NA),
  pos = c(tested$pos, length(sig$pos), length(intersect(sig$pos, sig$abs)), NA, NA, length(intersect(sig$pos, sig$neg))),
  neg = c(tested$neg, length(sig$neg), length(intersect(sig$neg, sig$abs)), NA, NA, NA),
  pos_or_neg = c(NA, length(pn), length(intersect(pn, sig$abs)), NA, length(setdiff(pn, sig$abs)), NA))
print(res)
fwrite(res, a[5], sep = "\t")
