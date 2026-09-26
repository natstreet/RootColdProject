#!/usr/bin/env Rscript
# orthogroup_hog_split.R — how many orthogroups the phylogenetic hierarchical orthogroups (HOGs) subdivide
# (Methods: the two orthology levels differ for the orthogroups that split into more than one HOG).
#
# The HOG table is N0.tsv of the same OrthoFinder v2.5.2 run as Orthogroups_20240613.tsv (27 species); a copy is
# committed as data/ComPlEx/Phylogenetic_Hierarchical_Orthogroups_N0.tsv.gz (from
# https://github.com/tallgran/conifer-JasonHill, data/N0.tsv.gz). The script first checks that the two tables
# are the same run: the same orthogroup ids with the same number of genes per species in each, and stops if not.
# Gene identifiers are not compared directly because the two files spell some of them differently.
#
# Inputs (run from a directory containing the unpacked figshare deposit `data/`, merged with this repo's data/):
#   data/ComPlEx/Orthogroups_20240613.tsv, data/ComPlEx/Phylogenetic_Hierarchical_Orthogroups_N0.tsv.gz
# Output: orthogroup_hog_split.tsv
# Usage:  Rscript scripts/analysis/orthogroup_hog_split.R

og  <- read.delim("data/ComPlEx/Orthogroups_20240613.tsv", check.names = FALSE, stringsAsFactors = FALSE)
hog <- read.delim(gzfile("data/ComPlEx/Phylogenetic_Hierarchical_Orthogroups_N0.tsv.gz"), check.names = FALSE,
                  stringsAsFactors = FALSE)
sp  <- setdiff(names(og), "Ortholog_Group")
stopifnot(identical(sp, names(hog)[-(1:3)]))
ngenes <- function(x) ifelse(is.na(x) | x == "", 0L, lengths(strsplit(x, ", ", fixed = TRUE)))

# same run: per-orthogroup, per-species gene counts identical
a <- sapply(sp, function(s) ngenes(og[[s]])); rownames(a) <- og$Ortholog_Group
b <- rowsum(sapply(sp, function(s) ngenes(hog[[s]])), hog$OG)
if (!setequal(rownames(a), rownames(b)) || !identical(a[rownames(b), ], b))
  stop("the HOG table is not from the same OrthoFinder run as Orthogroups_20240613.tsv")

n_hog <- table(hog$OG)
res <- data.frame(orthogroups = nrow(og), genes = sum(a),
                  orthogroups_with_more_than_one_HOG = sum(n_hog > 1),
                  percent = round(100 * sum(n_hog > 1) / nrow(og), 2))
write.table(res, "orthogroup_hog_split.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(res)
