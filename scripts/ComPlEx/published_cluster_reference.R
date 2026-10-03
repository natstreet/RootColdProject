#!/usr/bin/env Rscript
# published_cluster_reference.R — a record of the manuscript cluster (1-4) of each conserved-core orthogroup as
# published (deposit v2: clique_HMcluster_list.rda and the manuscript -> analysis map of that version), so a rebuilt
# clustering can be numbered consistently with the paper (manuscript_cluster_map.R).
# Usage (repository root, published deposit unpacked under <dir>/data): Rscript scripts/ComPlEx/published_cluster_reference.R <dir>
# Output: scripts/ComPlEx/reference/published_manuscript_clusters.tsv (orthogroup, manuscript_cluster)
d <- commandArgs(trailingOnly = TRUE)[1]
e <- new.env(); load(file.path(d, "data/ComPlEx/cliques/clique_HMcluster_list.rda"), envir = e)
x <- get("clique_cluster_order", e)
PUBLISHED_MAP <- c(`1` = 2, `2` = 4, `3` = 3, `4` = 1)      # manuscript -> analysis cluster, as published
out <- data.frame(orthogroup = sub("_[0-9]+$", "", x$ID),
                  manuscript_cluster = as.integer(names(PUBLISHED_MAP))[match(x$Cluster, PUBLISHED_MAP)])
write.table(out[order(out$manuscript_cluster, out$orthogroup), ], "scripts/ComPlEx/reference/published_manuscript_clusters.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("published manuscript clusters:", paste(names(table(out$manuscript_cluster)), table(out$manuscript_cluster), sep = "=", collapse = " "), "\n")
