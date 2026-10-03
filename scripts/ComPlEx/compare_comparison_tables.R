#!/usr/bin/env Rscript
# compare_comparison_tables.R — how much a rerun of the ComPlEx pairwise comparisons differs from a reference set
# (e.g. the published, deposited tables): per species pair, the orthologue pairs tested, the co-expressologs at
# Max.p.val < 0.1 (the threshold used for the cliques) and < 0.05, the overlap of the two co-expressolog sets
# (shared, only old, only new, Jaccard), and whether every column of the two tables is identical.
# Usage: Rscript compare_comparison_tables.R <reference dir> <new dir> <out.tsv>
#   both dirs hold comparison-<sp1>-<sp2>-abspearsonMR0.03-table.RData (new dir may hold them one folder deeper)
a <- commandArgs(trailingOnly = TRUE)
ld <- function(f) { e <- new.env(); load(f, envir = e); as.data.frame(get(ls(e)[1], e)) }
ref <- list.files(a[1], "^comparison-.*-table\\.RData$", full.names = TRUE)
new <- list.files(a[2], "^comparison-.*-table\\.RData$", full.names = TRUE, recursive = TRUE)
rows <- list()
for (f in ref) {
  g <- new[basename(new) == basename(f)]
  pair <- sub("^comparison-(.*)-abspearsonMR0.03-table\\.RData$", "\\1", basename(f))
  if (!length(g)) { rows[[pair]] <- data.frame(pair = pair, status = "no new table"); next }
  o <- ld(f); n <- ld(g[1])
  key <- function(x) paste(x$Species1, x$Species2)
  sig <- function(x, t) key(x)[as.numeric(x$Max.p.val) < t]
  s0 <- sig(o, 0.1); s1 <- sig(n, 0.1)
  rows[[pair]] <- data.frame(pair = pair, status = "compared",
    tested_old = nrow(o), tested_new = nrow(n),
    coexp_p0.1_old = length(s0), coexp_p0.1_new = length(s1),
    coexp_p0.05_old = length(sig(o, 0.05)), coexp_p0.05_new = length(sig(n, 0.05)),
    shared = length(intersect(s0, s1)), only_old = length(setdiff(s0, s1)), only_new = length(setdiff(s1, s0)),
    jaccard = round(length(intersect(s0, s1)) / length(union(s0, s1)), 3),
    identical_tables = isTRUE(all.equal(o[order(key(o)), ], n[order(key(n)), ], check.attributes = FALSE)))
}
res <- do.call(rbind, lapply(rows, function(r) { for (c in setdiff(names(rows[[which.max(sapply(rows, ncol))]]), names(r))) r[[c]] <- NA; r }))
write.table(res, a[3], sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
