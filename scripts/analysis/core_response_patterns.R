#!/usr/bin/env Rscript
# core_response_patterns.R — how strongly, and in which genotypes, the members of each conserved-core clique respond
# to cold.
#
# For each core clique and genotype, the member's DEG calls over the four time points (DE/de_results.tsv: adjusted
# P <= 0.05 and |log2FC| >= 1, as in the Methods) give its response: "up" or "down" when every DEG call has that sign,
# "mixed" when they differ, "none" when the member is not differentially expressed at any time point. The clique's
# lineage pattern is the smallest of these groups that contains every responding genotype:
#   Arabidopsis (Col-0, Ost-0), broadleaf trees (aspen, birch), conifers (spruce, pine),
#   angiosperms (Col-0, Ost-0, aspen, birch), otherwise "across lineages";
# "none" or "single genotype" when fewer than two genotypes respond. Direction is "same" when every responding member
# responds in the same single direction; "angiosperms vs conifers" or "Arabidopsis vs trees" (aspen, birch, spruce,
# pine) when the responding members on one side of that division all respond in one direction and those on the other
# side all in the opposite one; otherwise "mixed" (including a member whose own DEG calls differ in sign).
#
# Inputs (repository root): data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS, DE/de_results.tsv,
#   data/ComPlEx/cliques/clique_HMcluster_list.rda, data/ComPlEx/cliques/manuscript_cluster_map.tsv,
#   data/annotation/gene_aliases_20140331.txt
# Output: tables/core_response_patterns.tsv (one row per clique, ordered by number of responding genotypes, then by
#   the median over genotypes of the member's largest |log2FC|), and a summary on standard output.
# Usage:  Rscript scripts/analysis/core_response_patterns.R
suppressPackageStartupMessages(library(data.table))
SP <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
GROUPS <- list(Arabidopsis = c("col0", "ost0"), `broadleaf trees` = c("aspen", "birch"),
               conifers = c("spruce", "pine"), angiosperms = c("col0", "ost0", "aspen", "birch"))
ANGIO <- GROUPS$angiosperms; CONIF <- GROUPS$conifers
ARAB <- GROUPS$Arabidopsis; TREES <- c("aspen", "birch", "spruce", "pine")
strip <- function(x) sub("[.][0-9]+$", "", as.character(x))
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

core <- as.data.frame(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
de <- fread("DE/de_results.tsv")
hm <- ld("data/ComPlEx/cliques/clique_HMcluster_list.rda")
mm <- read.delim("data/ComPlEx/cliques/manuscript_cluster_map.tsv")
al <- read.delim("data/annotation/gene_aliases_20140331.txt", colClasses = "character")
sym <- function(g) { v <- al$symbol[toupper(al$locus_name) == strip(g)]; if (length(v) && nzchar(v[1])) v[1] else "" }

resp <- function(s, g) {
  x <- de[species == s & gene == g]
  if (nrow(x) != 4) stop("expected four time points for ", s, " ", g)
  d <- sign(x$log2FoldChange[x$deg])
  list(call = if (!length(d)) "none" else if (all(d > 0)) "up" else if (all(d < 0)) "down" else "mixed",
       maxabs = max(abs(x$log2FoldChange), na.rm = TRUE), n_tp = sum(x$deg))
}
rows <- lapply(seq_len(nrow(core)), function(i) {
  r <- lapply(SP, function(s) resp(s, core[[s]][i])); names(r) <- SP
  call <- vapply(r, `[[`, "", "call"); on <- SP[call != "none"]
  lineage <- if (length(on) == 0) "none" else if (length(on) == 1) "single genotype" else {
    hit <- names(GROUPS)[vapply(GROUPS, function(g) all(on %in% g), TRUE)]
    if (length(hit)) hit[which.min(lengths(GROUPS[hit]))] else "across lineages" }
  dirs <- call[on]
  direction <- if (length(on) < 2) "" else if (length(unique(dirs)) == 1 && dirs[1] != "mixed") "same" else {
    split <- function(A, B) { a <- unique(dirs[on %in% A]); b <- unique(dirs[on %in% B])
      length(a) == 1 && length(b) == 1 && setequal(c(a, b), c("up", "down")) }
    if (split(ANGIO, CONIF)) "angiosperms vs conifers" else if (split(ARAB, TREES)) "Arabidopsis vs trees" else "mixed" }
  cid <- core$cliqueID[i]
  as.data.table(c(list(OrthoGroup = core$OrthoGroup[i], cliqueID = cid,
             cluster = mm$manuscript_cluster[match(hm$Cluster[match(cid, hm$ID)], mm$analysis_cluster)],
             col0_gene = strip(core$col0[i]), col0_symbol = sym(core$col0[i]),
             n_genotypes_DE = length(on), DE_genotypes = paste(on, collapse = ","), lineage = lineage,
             direction = direction,
             median_max_abs_lfc = round(median(vapply(r, `[[`, 0, "maxabs")), 3)),
             setNames(as.list(call), paste0("response_", SP)),
             setNames(as.list(vapply(r, `[[`, 0L, "n_tp")), paste0("timepoints_DE_", SP))))
})
out <- rbindlist(rows)
setorder(out, -n_genotypes_DE, -median_max_abs_lfc)
dir.create("tables", showWarnings = FALSE)
fwrite(out, "tables/core_response_patterns.tsv", sep = "\t")

cat(sprintf("%d core cliques -> tables/core_response_patterns.tsv\n\nGenotypes in which the member is DE:\n", nrow(out)))
print(out[, .N, by = n_genotypes_DE][order(n_genotypes_DE)])
cat("\nLineage pattern (two or more genotypes DE):\n")
print(out[n_genotypes_DE >= 2, .N, by = .(lineage)][order(-N)])
cat("\nDirection (two or more genotypes DE):\n")
print(out[n_genotypes_DE >= 2, .N, by = .(direction)][order(-N)])
