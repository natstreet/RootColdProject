#!/usr/bin/env Rscript
# core_network_flag.R — flag the conserved-core clique members that have no edge in either network layer of their
# genotype (data/networks/network_<sp>_hrr100.tsv and _neg_hrr100.tsv), D047. Such a member is at the floor of the
# expression filters: it passed ComPlEx's filter but not the networks' (at least 1 CPM in at least k samples), so its
# co-expression evidence rests on very low expression. Appends, as the last column of conserved_core.tsv and
# conserved_core_flagged.tsv:
#   members_without_network_edges   the genotypes whose member has no network edge, ";"-separated; empty if none
# Every other column is left as it was. Stops if a table already has the column.
# Usage: Rscript scripts/deposit/core_network_flag.R <deposit data/open dir> <deposit data/networks dir>
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); O <- a[1]; N <- a[2]
G <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
in_net <- setNames(lapply(G, function(s) {
  x <- rbind(fread(file.path(N, sprintf("network_%s_hrr100.tsv", s)), select = 1:2),
             fread(file.path(N, sprintf("network_%s_neg_hrr100.tsv", s)), select = 1:2))
  unique(c(x[[1]], x[[2]])) }), G)
for (f in c("conserved_core.tsv", "conserved_core_flagged.tsv")) {
  p <- file.path(O, f)
  x <- read.delim(p, colClasses = "character", na.strings = character(0), check.names = FALSE, quote = "")
  if ("members_without_network_edges" %in% names(x)) stop(f, " already has members_without_network_edges")
  x$members_without_network_edges <- vapply(seq_len(nrow(x)), function(i)
    paste(G[vapply(G, function(s) !(x[[s]][i] %in% in_net[[s]]), TRUE)], collapse = ";"), "")
  write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  fl <- x[x$members_without_network_edges != "", c("OrthoGroup", "members_without_network_edges")]
  cat(sprintf("%s: %d of %d rows flagged%s\n", f, nrow(fl), nrow(x),
              if (nrow(fl)) paste0(" (", paste(sprintf("%s: %s", fl$OrthoGroup, fl$members_without_network_edges), collapse = ", "), ")") else ""))
}
