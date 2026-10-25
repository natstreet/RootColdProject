#!/usr/bin/env Rscript
# fig4_captions_family_level.R — Fig. 4 block captions from the family-level Table S3 (D090), by a stated rule:
# per cluster, the N most significant family-level terms (lowest padj) that pass Fisher, the size-matched test and the
# majority rule; if fewer pass all three, the next from terms passing Fisher and size-matched, then Fisher and majority,
# then Fisher alone. A term whose
# set of families equals that of a term already chosen for the cluster is skipped (parent/child terms carried by the same
# families). The basis of each caption is recorded. Captions are provisional until NRS/SCI approve them.
# Input: the family-level Table S3 tsv (tableS3_arabidopsis.R). Output: a caption table (cluster, GO.ID, Term, basis, padj).
# Usage: Rscript scripts/figures/fig4_captions_family_level.R <TableS3_family_level.tsv> <out tsv> [N, default 2]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); N <- if (length(a) >= 3) as.integer(a[3]) else 2L
f <- fread(a[1]); setnames(f, "Cluster (manuscript)", "cluster")
f[, basis := fifelse(size_matched & majority, "Fisher, size-matched and majority",
              fifelse(size_matched, "Fisher and size-matched", fifelse(majority, "Fisher and majority", "Fisher only")))]
f[, rank_basis := match(basis, c("Fisher, size-matched and majority", "Fisher and size-matched", "Fisher and majority", "Fisher only"))]
out <- rbindlist(lapply(sort(unique(f$cluster)), function(k) {
  d <- f[cluster == k][order(rank_basis, padj, GO.ID)]; keep <- integer(0); seen <- character(0)
  for (i in seq_len(nrow(d))) { if (length(keep) == N) break
    if (d$family_ids[i] %in% seen) next
    keep <- c(keep, i); seen <- c(seen, d$family_ids[i]) }
  d[keep, .(cluster, GO.ID, Term, basis, padj)]
}))
fwrite(out, a[2], sep = "\t"); print(out)
