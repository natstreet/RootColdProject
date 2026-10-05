#!/usr/bin/env Rscript
# value_changes.R — the tracked numbers of the article (Figure 1 panels, Figure 4 Venn regions, complete-clique
# summary, common orthogroups, ecotype DEG overlap, Figure S5) compared between two producer runs. Each tree is a
# working directory of scripts/reproduce_figures_tables.sh plus complete_cliques_per_orthogroup.R,
# common_orthogroups.R and ecotype_deg_overlap.R; a value file is read from the tree or its results/.
# Columns: source, item, published (first tree), rebuilt (second tree), changed.
# Usage: Rscript scripts/analysis/value_changes.R <first tree> <second tree> <out.tsv>
P <- commandArgs(TRUE)[1]; T <- commandArgs(TRUE)[2]; out <- commandArgs(TRUE)[3]
f <- function(d, n) { p <- c(file.path(d, n), file.path(d, "results", n)); read.delim(p[file.exists(p)][1], check.names = FALSE, colClasses = "character") }
rows <- list()
add <- function(src, key, old, new) rows[[length(rows) + 1]] <<- data.frame(source = src, item = key, published = old, rebuilt = new)
a <- f(P, "Figure1_panel_values.tsv"); b <- f(T, "Figure1_panel_values.tsv")
k <- function(x) paste(x$panel, x$species, x$category); m <- match(k(a), k(b))
for (i in seq_len(nrow(a))) { add("Figure 1", paste(a$panel[i], a$species[i], a$category[i], "(count)"), a$count[i], b$count[m[i]])
  if (nzchar(a$percent[i])) add("Figure 1", paste(a$panel[i], a$species[i], a$category[i], "(%)"), a$percent[i], b$percent[m[i]]) }
for (n in c("Figure4_venn_regions.tsv", "complete_cliques_summary.tsv", "common_orthogroups.tsv", "ecotype_deg_overlap.tsv")) {
  a <- f(P, n); b <- f(T, n); for (i in seq_len(nrow(a))) add(sub(".tsv$", "", n), a[[1]][i], a[[2]][i], b[[2]][match(a[[1]][i], b[[1]])]) }
a <- f(P, "Figure_S5_values.tsv"); b <- f(T, "Figure_S5_values.tsv")
for (i in seq_len(nrow(a))) add("Figure S5", paste(a$comparison[i], a$quantity[i]), a$value[i], b$value[k2 <- match(paste(a$comparison[i], a$quantity[i]), paste(b$comparison, b$quantity))])
r <- do.call(rbind, rows); r$changed <- r$published != r$rebuilt
write.table(r, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(r), "values compared;", sum(r$changed, na.rm = TRUE), "changed\n")
print(r[r$source == "Figure 1" & r$changed, c("item", "published", "rebuilt")], row.names = FALSE)
