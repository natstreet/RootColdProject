#!/usr/bin/env Rscript
# orthology_class_tables.R — descriptive tables of the orthology class over the conserved core (D043; no tests).
# The class is a gene-tree call that probably overstates paralogue-carried cases, and it depends strongly on family size,
# so every table is to be read with the size table beside it.
#
# Inputs (working directory): results/orthology_class.tsv (scripts/analysis/orthology_class.R),
#   tables/core_response_patterns.tsv (scripts/analysis/core_response_patterns.R: manuscript cluster, lineage pattern,
#   direction), data/annotation/gene_aliases_20140331.txt
# Outputs: results/orthology_class_by_size.tsv       class x expressed-gene band
#          results/orthology_class_by_cluster.tsv    class x manuscript cluster
#          results/orthology_class_by_response.tsv   class x lineage pattern x direction
#          results/paralogue_carried_small.tsv       paralogue-carried core orthogroups with <= 20 expressed genes: the
#                                                    representative clique's members and the Col-0 member's symbol
# Usage: Rscript scripts/analysis/orthology_class_tables.R
suppressPackageStartupMessages(library(data.table))
BANDS <- c(0, 10, 20, 50, 100, 200, Inf); SMALL <- 20
G <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
oc <- fread("results/orthology_class.tsv")[counted_as_conserved == TRUE]
rp <- fread("tables/core_response_patterns.tsv")
if (!setequal(oc$OrthoGroup, rp$OrthoGroup)) stop("orthology classes and response patterns cover different orthogroups")
d <- merge(oc, rp[, .(OrthoGroup, cluster, lineage, direction)], by = "OrthoGroup")
d[, size_band := cut(expressed_genes, BANDS, labels = c("<=10", "11-20", "21-50", "51-100", "101-200", ">200"))]

wide <- function(x, by) {
  t <- dcast(x[, .N, by = c(by, "orthology_class")], as.formula(paste(paste(by, collapse = " + "), "~ orthology_class")),
             value.var = "N", fill = 0)
  for (k in c("orthologue-carried", "paralogue-carried")) if (!k %in% names(t)) t[[k]] <- 0L
  t[, total := `orthologue-carried` + `paralogue-carried`]
  t[, pct_orthologue_carried := round(100 * `orthologue-carried` / total, 1)]
  t
}
fwrite(wide(d, "size_band")[order(size_band)], "results/orthology_class_by_size.tsv", sep = "\t")
fwrite(wide(d, "cluster")[order(cluster)], "results/orthology_class_by_cluster.tsv", sep = "\t")
fwrite(wide(d, c("lineage", "direction"))[order(-total)], "results/orthology_class_by_response.tsv", sep = "\t")

core <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
al <- fread("data/annotation/gene_aliases_20140331.txt", colClasses = "character")
strip <- function(x) toupper(sub("[.][0-9]+$", "", x))
small <- merge(d[orthology_class == "paralogue-carried" & expressed_genes <= SMALL,
                 .(OrthoGroup, expressed_genes, n_complete_cliques, cluster, lineage, direction)],
               core[, c("OrthoGroup", G), with = FALSE], by = "OrthoGroup")
small[, col0_symbol := vapply(col0, function(g) paste(unique(al$symbol[toupper(al$locus_name) == strip(g)]), collapse = "/"), "")]
setorder(small, expressed_genes, OrthoGroup)
fwrite(small, "results/paralogue_carried_small.tsv", sep = "\t")

for (f in c("by_size", "by_cluster", "by_response")) { cat("==", f, "\n"); print(fread(sprintf("results/orthology_class_%s.tsv", f)), row.names = FALSE) }
cat("== paralogue-carried with <=", SMALL, "expressed genes:", nrow(small), "\n")
print(small[, .(OrthoGroup, expressed_genes, cluster, col0, col0_symbol)], row.names = FALSE)
