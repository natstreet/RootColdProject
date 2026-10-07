#!/usr/bin/env Rscript
# test_rebuilt_objects.R — runs the producers of four analysis objects into a temporary folder and compares each
# result with the copy under data/ (the deposited object). Exits 1 if any differs.
#
#   ones_and_zeros_all.RData               scripts/ComPlEx/ones_and_zeros.R      identical()
#   conserved_genes_COMPLETE_filtered.RDS  scripts/ComPlEx/size_rule.R, then core_cliques.R with its table (D026):
#                                          identical() after ordering by orthogroup, and the same for the
#                                          _flagged.RDS object of the orthogroups the rule does not count
#   clique_HMcluster_list.rda              scripts/ComPlEx/clique_clusters.R     identical()
#   og_summary.RData                       scripts/DE/og_summary.R               identical() after ordering by orthogroup
# Usage (repository root, deposit unpacked under data/): Rscript scripts/tests/test_rebuilt_objects.R
suppressPackageStartupMessages(library(dplyr))
ld  <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
tmp <- tempfile("rebuilt_"); dir.create(tmp)
run <- function(script, out, ...) {
  st <- system2("Rscript", c(script, out, ...), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(st, "status")) || !file.exists(out)) stop(script, " failed:\n", paste(st, collapse = "\n"))
  cat(tail(st, 1), "\n")
  out
}
fails <- 0
check <- function(ok, what) { cat(if (ok) "ok  " else "FAIL", what, "\n"); if (!ok) fails <<- fails + 1 }

x <- ld(run("scripts/ComPlEx/ones_and_zeros.R", file.path(tmp, "ones_and_zeros_all.RData")))
y <- ld("data/ComPlEx/ones_and_zeros_all.RData")
check(identical(x, y), sprintf("ones_and_zeros_all: %d x %d, cells differing %s", nrow(x), ncol(x),
                               if (identical(dim(x), dim(y))) sum(as.matrix(x) != as.matrix(y)) else "dimensions differ"))

rule <- file.path(tmp, "orthogroup_size_rule")
st <- system2("Rscript", c("scripts/ComPlEx/size_rule.R", rule), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(st, "status")) || !file.exists(paste0(rule, ".tsv"))) stop("size_rule.R failed:\n", paste(st, collapse = "\n"))
x <- readRDS(run("scripts/ComPlEx/core_cliques.R", file.path(tmp, "conserved_genes_COMPLETE_filtered.RDS"), paste0(rule, ".tsv")))
y <- readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS")
check(identical(arrange(x, OrthoGroup)[, names(y)], arrange(y, OrthoGroup)),
      sprintf("conserved_genes_COMPLETE_filtered: %d rebuilt vs %d deposited orthogroups", nrow(x), nrow(y)))
x <- readRDS(file.path(tmp, "conserved_genes_COMPLETE_filtered_flagged.RDS"))
y <- readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered_flagged.RDS")
check(identical(arrange(x, OrthoGroup)[, names(y)], arrange(y, OrthoGroup)),
      sprintf("conserved_genes_COMPLETE_filtered_flagged: %d rebuilt vs %d deposited orthogroups", nrow(x), nrow(y)))

x <- ld(run("scripts/ComPlEx/clique_clusters.R", file.path(tmp, "clique_HMcluster_list.rda")))
y <- ld("data/ComPlEx/cliques/clique_HMcluster_list.rda")
check(identical(x, y), sprintf("clique_HMcluster_list: %d cliques, cluster sizes %s", nrow(x), paste(table(x$Cluster), collapse = "/")))

x <- ld(run("scripts/DE/og_summary.R", file.path(tmp, "og_summary.RData")))
y <- ld("data/DEGs/og_summary.RData")
check(identical(arrange(x, OrthoGroup), arrange(y, OrthoGroup)),
      sprintf("og_summary: %d orthogroups, %d with a DEG in all six genotypes", nrow(x), sum(x$n_species_with_DEGs == 6)))

unlink(tmp, recursive = TRUE)
cat(fails, "of 4 differ\n")
quit(status = if (fails) 1 else 0)
