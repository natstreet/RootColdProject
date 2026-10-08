#!/usr/bin/env Rscript
# motif_pairs.R — orthologue pairs and gene lists for the promoter motif track (DECISIONS.md D062 item 4, D063; R040 as
# amended in A035).
#
# Pairs: every gene pair ComPlEx tested within an orthogroup that has a complete clique, in each of the 10 comparisons
# among col0, aspen, birch, spruce and pine (Ost-0 shares the Arabidopsis promoters, so its comparisons are left out).
#   co-expressolog      Max.p.val < 0.1 (the threshold used for complete cliques and the D052 test)
#   non-co-expressolog  tested, Max.p.val >= 0.1; pairs ComPlEx did not test (a gene filtered out) are not used
# The orthogroup set is every orthogroup with a complete clique, a superset of the final core; the test script restricts
# to the final core, so promoters and motif scans do not depend on R039.
# Inputs (working tree): data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS,
#   data/ComPlEx/PairwiseComparisons/comparison-<s1>-<s2>-abspearsonMR0.03-table.RData
# Outputs (<out dir>): motif_pairs.tsv.gz (s1, s2, OrthoGroup, gene1, gene2, max_p, coexpressolog);
#   genes_<genotype>.txt (one gene per line, the genes of those pairs)
# Usage (working tree): Rscript scripts/motifs/motif_pairs.R <out dir>
suppressPackageStartupMessages(library(data.table))
OUT <- commandArgs(trailingOnly = TRUE)[1]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P_THR <- 0.1
G <- c("col0", "aspen", "birch", "spruce", "pine")
ogs <- unique(as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS"))$OrthoGroup)
ld <- function(f) { e <- new.env(); load(f, envir = e); as.data.table(get(ls(e)[1], e)) }
cmp <- CJ(i = seq_along(G), j = seq_along(G))[i < j]
pairs <- rbindlist(lapply(seq_len(nrow(cmp)), function(k) {
  a <- G[cmp$i[k]]; b <- G[cmp$j[k]]
  f <- sprintf("data/ComPlEx/PairwiseComparisons/comparison-%s-%s-abspearsonMR0.03-table.RData", a, b)
  if (!file.exists(f)) { f <- sprintf("data/ComPlEx/PairwiseComparisons/comparison-%s-%s-abspearsonMR0.03-table.RData", b, a); s <- c(b, a) } else s <- c(a, b)
  x <- ld(f)[OrthoGroup %in% ogs]
  x[, .(s1 = s[1], s2 = s[2], OrthoGroup, gene1 = Species1, gene2 = Species2, max_p = as.numeric(Max.p.val))] }))
if (pairs[, uniqueN(paste(s1, s2))] != 10) stop("expected 10 comparisons")
pairs[, coexpressolog := max_p < P_THR]
fwrite(pairs, file.path(OUT, "motif_pairs.tsv.gz"), sep = "\t")
for (g in G) writeLines(sort(unique(c(pairs[s1 == g, gene1], pairs[s2 == g, gene2]))), file.path(OUT, sprintf("genes_%s.txt", g)))
print(pairs[, .(pairs = .N, coexpressologs = sum(coexpressolog), orthogroups = uniqueN(OrthoGroup)), by = .(s1, s2)], row.names = FALSE)
cat("genes per genotype:", sapply(G, function(g) length(readLines(file.path(OUT, sprintf("genes_%s.txt", g))))), "\n")
