#!/usr/bin/env Rscript
# core_size_columns.R — add the size rule (D026) and the orthology class (D043) to the deposit's conserved-core table,
# and write the table of the orthogroups the size rule does not count (reported, not hidden).
#
# Columns appended to conserved_core.tsv (every other column left byte for byte as it was):
#   expressed_genes, complete_cliques, p_core_given_size     results/orthogroup_size_rule.tsv (scripts/ComPlEx/size_rule.R)
#   orthology_class, n_orthology_consistent_cliques, representative_orthology_consistent, best_consistent_clique,
#   best_col0 ... best_pine, orthology_basis                 results/orthology_class.tsv (scripts/analysis/orthology_class.R)
# conserved_core_flagged.tsv: the orthogroups with a complete clique that the size rule does not count, with their
# representative clique (data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered_flagged.RDS, core_cliques.R) and the
# same appended columns. They have no cluster: the clusters are of the core only.
# Stops if conserved_core.tsv already has the columns, if an orthogroup lacks a value, or if the core table and the size
# rule disagree on which orthogroups are counted.
#
# Usage (working directory with results/ and data/): Rscript scripts/deposit/core_size_columns.R <deposit data/open dir>
open <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(open), dir.exists(open))
p <- file.path(open, "conserved_core.tsv")
x <- read.delim(p, colClasses = "character", na.strings = character(0), check.names = FALSE, quote = "")
rule <- read.delim("results/orthogroup_size_rule.tsv", colClasses = "character")
oc <- read.delim("results/orthology_class.tsv", colClasses = "character", na.strings = character(0))
G <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
oc_cols <- c("orthology_class", "n_orthology_consistent_cliques", "representative_orthology_consistent",
             "best_consistent_clique", paste0("best_", G), "orthology_basis")
size_cols <- c("expressed_genes", "complete_cliques", "p_core_given_size")
if (any(c(size_cols, oc_cols) %in% names(x))) stop("conserved_core.tsv already has the size or orthology columns")
if (!setequal(x$OrthoGroup, rule$OrthoGroup[rule$counted_as_conserved == "TRUE"]))
  stop("conserved_core.tsv and the size rule disagree on the counted orthogroups")
add <- function(t) {
  m <- match(t$OrthoGroup, rule$OrthoGroup); k <- match(t$OrthoGroup, oc$OrthoGroup)
  if (anyNA(m) || anyNA(k)) stop("no size rule or orthology class for ", paste(t$OrthoGroup[is.na(m) | is.na(k)], collapse = ", "))
  t[size_cols] <- rule[m, size_cols]; t[oc_cols] <- oc[k, oc_cols]; t
}
x <- add(x)
write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

fl <- as.data.frame(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered_flagged.RDS"))
fl[] <- lapply(fl, as.character)
if (!setequal(fl$OrthoGroup, rule$OrthoGroup[rule$has_complete_clique == "TRUE" & rule$counted_as_conserved == "FALSE"]))
  stop("the flagged object and the size rule disagree")
fl <- add(fl[order(-as.numeric(rule$expressed_genes[match(fl$OrthoGroup, rule$OrthoGroup)])), ])
write.table(fl, file.path(open, "conserved_core_flagged.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("conserved_core.tsv: %d orthogroups, %d columns added; conserved_core_flagged.tsv: %d orthogroups (%s)\n",
            nrow(x), length(c(size_cols, oc_cols)), nrow(fl), paste(fl$OrthoGroup, collapse = ", ")))
