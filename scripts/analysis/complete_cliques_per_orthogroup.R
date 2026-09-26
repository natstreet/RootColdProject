#!/usr/bin/env Rscript
# complete_cliques_per_orthogroup.R — how many complete cliques each conserved-core orthogroup yields (Methods).
#
# A complete clique has one gene per genotype with every pair a co-expressolog. The conserved core is one
# representative complete clique per orthogroup (116 orthogroups); this counts, for each of those orthogroups, all
# its complete cliques.
#
# Inputs (run from a directory containing the unpacked figshare deposit `data/`, merged with this repo's data/):
#   data/ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS   every complete clique (cliqueID per gene pair)
#   data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS  the conserved core (one clique per orthogroup)
# Output: complete_cliques_per_orthogroup.tsv (per core orthogroup: number of complete cliques) and
#         complete_cliques_summary.tsv (core orthogroups, those with more than one complete clique, median count)
# Usage:  Rscript scripts/analysis/complete_cliques_per_orthogroup.R

cg   <- as.data.frame(readRDS("data/ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS"))
core <- as.data.frame(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
stopifnot(!anyDuplicated(core$OrthoGroup), all(core$cliqueID %in% cg$cliqueID))

n <- tapply(cg$cliqueID, cg$OrthoGroup, function(x) length(unique(x)))
per <- data.frame(OrthoGroup = core$OrthoGroup, complete_cliques = as.integer(n[core$OrthoGroup]))
per <- per[order(per$OrthoGroup), ]
write.table(per, "complete_cliques_per_orthogroup.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

summ <- data.frame(
  quantity = c("conserved-core orthogroups", "core orthogroups with more than one complete clique",
               "median complete cliques per core orthogroup", "minimum", "maximum"),
  value = c(nrow(per), sum(per$complete_cliques > 1), stats::median(per$complete_cliques),
            min(per$complete_cliques), max(per$complete_cliques)))
write.table(summ, "complete_cliques_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE)
