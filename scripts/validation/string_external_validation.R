#!/usr/bin/env Rscript
# string_external_validation.R — independent test that the conserved-core clusters are coherent co-regulated modules
# (Fig. S4b).
#
# Tests whether the Arabidopsis members of the conserved core (the clique clusters of Fig. 4) are more co-expressed
# within clusters than between clusters in STRING v12 (an external A. thaliana co-expression compendium). This
# validates the cluster structure as co-regulation, complementing the permutation test that establishes the
# cross-species conservation itself.
#
# Inputs (download once from STRING, https://string-db.org, into the working directory):
#   3702.protein.links.detailed.v12.0.txt.gz   (co-expression channel = the "coexpression" column)
#   3702.protein.aliases.v12.0.txt.gz          (STRING protein -> TAIR locus)
#   core_clusters.tsv                          (Col0 gene + Cluster of each conserved-core clique; written by
#                                              scripts/validation/make_core_clusters.R from the deposit)
# Method: the co-expression network at subscore >= THRESHOLD among proteins mappable to TAIR loci; among the
# conserved-core Arabidopsis genes (with cluster labels), count co-expression links that fall within a cluster; compare
# with N_PERM permutations that re-assign the same genes to clusters at random (cluster sizes kept). STRING lists each
# pair in both directions; each gene pair is counted once.
# Outputs: string_null.csv (column null), string_obs.txt — the inputs of make_FigS4.R
# Usage:   Rscript scripts/validation/string_external_validation.R [clusters tsv] [output prefix]
#          (defaults core_clusters.tsv and no prefix; with a prefix, the outputs are <prefix>_string_null.csv and
#          <prefix>_string_obs.txt, e.g. for the sample-stable and sample-sensitive orthogroups separately, D066)

suppressMessages(library(data.table))
THRESHOLD <- 400; N_PERM <- 5000; set.seed(1)
LINKS <- "3702.protein.links.detailed.v12.0.txt.gz"; ALIASES <- "3702.protein.aliases.v12.0.txt.gz"

string_tables <- function() {
  al <- fread(cmd = paste("gzip -dc", shQuote(ALIASES)), sep = "\t", header = FALSE, skip = 1, select = 1:2, col.names = c("prot", "alias"), quote = "")
  al[, alias := toupper(alias)]
  al <- al[grepl("^AT[1-5MC]G[0-9]{5}$", alias)]
  at2s <- al[!duplicated(alias)]                                   # first STRING protein listed for each locus
  lk <- fread(cmd = paste("gzip -dc", shQuote(LINKS)), sep = " ", select = c("protein1", "protein2", "coexpression"))
  lk <- lk[protein1 < protein2 & coexpression >= THRESHOLD &
           protein1 %in% at2s$prot & protein2 %in% at2s$prot, .(a = protein1, b = protein2)]
  list(at2s = setNames(at2s$prot, at2s$alias), links = lk)
}

st <- string_tables()
ARGS <- commandArgs(trailingOnly = TRUE); CLF <- if (length(ARGS) > 0) ARGS[1] else "core_clusters.tsv"; PFX <- if (length(ARGS) > 1) paste0(ARGS[2], "_") else ""
cl <- fread(CLF)[!is.na(Col0) & !is.na(Cluster)]
prot <- st$at2s[toupper(cl$Col0)]
nodes_in_net <- unique(c(st$links$a, st$links$b))
ok <- !is.na(prot) & prot %in% nodes_in_net
g2c <- setNames(as.integer(cl$Cluster[ok]), prot[ok])
g2c <- g2c[!duplicated(names(g2c))]
e <- st$links[a %in% names(g2c) & b %in% names(g2c)]
nodes <- sort(names(g2c)); lab <- unname(g2c[nodes])
ei <- match(e$a, nodes); ej <- match(e$b, nodes)
obs <- sum(lab[ei] == lab[ej])
null <- vapply(seq_len(N_PERM), function(i) { p <- sample(lab); sum(p[ei] == p[ej]) }, integer(1))
pval <- (sum(null >= obs) + 1) / (N_PERM + 1)

write.csv(data.frame(null = null), paste0(PFX, "string_null.csv"), row.names = FALSE, quote = FALSE)
writeLines(as.character(obs), paste0(PFX, "string_obs.txt"))
within <- lab[ei][lab[ei] == lab[ej]]
cat(sprintf("core genes with cluster + STRING: %d; core-core gene pairs: %d\n", length(nodes), nrow(e)))
cat(sprintf("observed within-cluster gene pairs: %d\n", obs))
cat("within-cluster pairs by cluster:", paste(sprintf("%s: %d", names(table(within)), table(within)), collapse = ", "), "\n")
cat(sprintf("null mean=%.1f sd=%.1f max=%d\n", mean(null), sqrt(mean((null - mean(null))^2)), max(null)))
cat(sprintf("fold=%.1fx ; empirical P=%.4g\n", obs / max(mean(null), 1e-9), pval))
