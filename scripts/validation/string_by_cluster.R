#!/usr/bin/env Rscript
# string_by_cluster.R — the Fig. S6b STRING test broken down by clique cluster, and repeated without manuscript
# cluster 2 (the translation / ribosome-biogenesis cluster), to ask whether cluster coherence in STRING holds beyond
# that cluster.
#
# Same data and settings as string_external_validation.R: STRING v12 A. thaliana co-expression channel, subscore
# >= 400, each gene pair counted once, core genes that have at least one such link, 5,000 permutations of cluster
# labels (cluster sizes kept), seed 1.
#   1. per cluster, from the full permutation: observed and null-mean within-cluster pairs, empirical P
#   2. manuscript cluster 2 excluded: the same test on the remaining genes
#   3. the within-cluster pairs of every cluster, with Arabidopsis gene symbols
#
# Run from a directory holding data/ (the deposit plus this repository's data/), core_clusters.tsv (from
# make_core_clusters.R) and the two STRING files; manuscript cluster numbers are read from MAP in
# scripts/GO/make_TableS4.R.
# Outputs: string_by_cluster.tsv (tests 1 and 2), string_within_cluster_pairs.tsv (3)
# Usage:   Rscript scripts/validation/string_by_cluster.R

suppressMessages(library(data.table))
THRESHOLD <- 400; N_PERM <- 5000; SEED <- 1
LINKS <- "3702.protein.links.detailed.v12.0.txt.gz"; ALIASES <- "3702.protein.aliases.v12.0.txt.gz"

src <- readLines("scripts/GO/make_TableS4.R")
MAP <- eval(parse(text = sub("^MAP <- ", "", regmatches(src, regexpr("MAP <- c\\(.*\\)", src)))))  # manuscript -> analysis
an2ms <- setNames(as.integer(names(MAP)), MAP)

al <- fread(cmd = paste("gzip -dc", shQuote(ALIASES)), sep = "\t", header = FALSE, skip = 1, select = 1:2,
            col.names = c("prot", "alias"), quote = "")
al[, alias := toupper(alias)]
al <- al[grepl("^AT[1-5MC]G[0-9]{5}$", alias)][!duplicated(alias)]  # first STRING protein listed for each locus
at2s <- setNames(al$prot, al$alias)
lk <- fread(cmd = paste("gzip -dc", shQuote(LINKS)), sep = " ", select = c("protein1", "protein2", "coexpression"))
lk <- unique(lk[protein1 < protein2 & coexpression >= THRESHOLD & protein1 %in% al$prot & protein2 %in% al$prot,
                .(a = protein1, b = protein2)])
linked <- unique(c(lk$a, lk$b))

cl <- fread("core_clusters.tsv")
s <- at2s[toupper(cl$Col0)]; ok <- !is.na(s) & s %in% linked
g2c <- setNames(unname(an2ms[as.character(cl$Cluster[ok])]), s[ok])              # manuscript cluster number
g2c <- g2c[!duplicated(names(g2c))]
# a STRING protein can carry several Arabidopsis loci (paralogues); label each protein with the core gene it stands for
s2at <- setNames(toupper(cl$Col0[ok]), s[ok]); s2at <- s2at[!duplicated(names(s2at))]
core <- lk[a %in% names(g2c) & b %in% names(g2c)]

perm_test <- function(g2c, pairs) {
  set.seed(SEED)
  nodes <- sort(names(g2c)); lab <- unname(g2c[nodes]); ei <- match(pairs$a, nodes); ej <- match(pairs$b, nodes)
  ks <- sort(unique(lab))
  counts <- function(l) { same <- l[ei] == l[ej]; c(all = sum(same), vapply(ks, function(k) sum(same & l[ei] == k), 0L)) }
  obs <- counts(lab)
  null <- vapply(seq_len(N_PERM), function(i) counts(sample(lab)), numeric(length(ks) + 1))
  list(obs = obs, null = null, ks = ks, size = vapply(ks, function(k) sum(lab == k), 0L))
}
p <- function(o, n) (sum(n >= o) + 1) / (N_PERM + 1)
sd0 <- function(x) sqrt(mean((x - mean(x))^2))                     # population sd, as numpy's std
row <- function(test, cluster, genes, pairs, o, n, fold = TRUE)
  data.frame(test = test, cluster = cluster, genes = genes, core_pairs = pairs, observed = o,
             null_mean = round(mean(n), 2), null_sd = round(sd0(n), 2),
             fold = if (fold && mean(n) > 0) as.character(round(o / mean(n), 2)) else "", empirical_P = p(o, n))

r <- perm_test(g2c, core)
rows <- list(row("all clusters", "all", length(g2c), nrow(core), r$obs[1], r$null[1, ]))
for (i in seq_along(r$ks))
  rows[[length(rows) + 1]] <- row("all clusters, per cluster", r$ks[i], r$size[i], "", r$obs[i + 1], r$null[i + 1, ])
g2c_x <- g2c[g2c != 2]; core_x <- core[a %in% names(g2c_x) & b %in% names(g2c_x)]
r <- perm_test(g2c_x, core_x)
rows[[length(rows) + 1]] <- row("cluster 2 excluded", "1,3,4", length(g2c_x), nrow(core_x), r$obs[1], r$null[1, ])
for (i in seq_along(r$ks))
  rows[[length(rows) + 1]] <- row("cluster 2 excluded, genes per cluster", r$ks[i], r$size[i], "", r$obs[i + 1],
                                  r$null[i + 1, ], fold = FALSE)
out <- do.call(rbind, rows)
write.table(out, "string_by_cluster.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

al2 <- fread("data/annotation/gene_aliases_20140331.txt", sep = "\t", quote = "")
al2 <- al2[!duplicated(toupper(locus_name))]
sym <- setNames(al2$symbol, toupper(al2$locus_name))
w <- core[g2c[a] == g2c[b]]
w <- data.frame(cluster = unname(g2c[w$a]), gene1 = unname(s2at[w$a]), gene2 = unname(s2at[w$b]))
w <- w[order(w$cluster, w$gene1, w$gene2), ]
w$symbol1 <- ifelse(is.na(sym[w$gene1]), "", sym[w$gene1]); w$symbol2 <- ifelse(is.na(sym[w$gene2]), "", sym[w$gene2])
write.table(w[, c("cluster", "gene1", "symbol1", "gene2", "symbol2")], "string_within_cluster_pairs.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
