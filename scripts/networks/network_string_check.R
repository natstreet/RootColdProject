#!/usr/bin/env Rscript
# network_string_check.R — functional support in STRING v12 for A. thaliana co-expression networks built with
# build_coexpression_network.R (e.g. the positive, negative and absolute networks of one experiment).
#
# Two STRING evidence channels, each at score >= 400:
#   independent:  max(experimental, database), evidence that does not come from expression
#   coexpression: STRING's co-expression channel (an external A. thaliana compendium)
# Measures, per network:
#   edges         the fraction of the network's edges whose two genes are linked in STRING, over the fraction among
#                 random pairs of the same genes (enrichment); for the absolute network also split by the sign of
#                 the edge's jackknifed correlation
#   neighbourhood for a sample of genes (the same genes in every network: present in all of them with at least
#                 MIN_DEG partners), the density of STRING links among the gene's partners, over its mean in random sets
#                 of the same size drawn from the network's genes (N_NULL sets); summarised as the mean observed density
#                 over the mean random density across the sampled genes, and the share of genes with any link among
#                 their partners
# Usage: Rscript network_string_check.R <STRING dir> <out prefix> <label>=<network.tsv> [<label>=<network.tsv> ...]
suppressMessages({ library(data.table); library(Matrix) })
a <- commandArgs(trailingOnly = TRUE); SD <- a[1]; OUT <- a[2]
nets <- lapply(strsplit(a[-(1:2)], "=", fixed = TRUE), function(x) list(label = x[1], file = x[2]))
THR <- 400; N_GENES <- 2000L; N_NULL <- 20L; MIN_DEG <- 10L; N_RANDOM_PAIRS <- 2e6
set.seed(1)

E <- lapply(nets, function(n) fread(n$file)); names(E) <- vapply(nets, `[[`, "", "label")
genes <- sort(unique(unlist(lapply(E, function(e) c(e$source, e$target)))))
al <- fread(cmd = paste("gzip -dc", shQuote(file.path(SD, "3702.protein.aliases.v12.0.txt.gz"))), sep = "\t", header = FALSE,
            skip = 1, select = 1:2, col.names = c("prot", "alias"), quote = "")
al[, alias := toupper(alias)]; al <- al[grepl("^AT[1-5MC]G[0-9]{5}$", alias)][!duplicated(alias)]
g2p <- setNames(al$prot, al$alias)
sg <- genes[toupper(sub("[.][0-9]+$", "", genes)) %in% names(g2p)]       # network genes present in STRING
pid <- g2p[toupper(sub("[.][0-9]+$", "", sg))]
lk <- fread(cmd = paste("gzip -dc", shQuote(file.path(SD, "3702.protein.links.detailed.v12.0.txt.gz"))), sep = " ",
            select = c("protein1", "protein2", "coexpression", "experimental", "database"))
lk <- lk[protein1 %in% pid & protein2 %in% pid]
gi <- function(p) match(p, pid)
adj <- function(keep) { i <- gi(lk$protein1[keep]); j <- gi(lk$protein2[keep])
  m <- sparseMatrix(i = c(i, j), j = c(j, i), x = 1, dims = c(length(sg), length(sg))); m@x[] <- 1; m }
A <- list(independent = adj(pmax(lk$experimental, lk$database) >= THR), coexpression = adj(lk$coexpression >= THR))
idx <- function(g) match(g, sg)

## ---- edges ----
rp <- data.table(i = sample(length(sg), N_RANDOM_PAIRS, TRUE), j = sample(length(sg), N_RANDOM_PAIRS, TRUE))[i != j]
base <- vapply(A, function(M) mean(M[cbind(rp$i, rp$j)]), numeric(1))
edges <- rbindlist(lapply(names(E), function(l) {
  e <- E[[l]][, .(i = idx(source), j = idx(target), r = jackknife_r)][!is.na(i) & !is.na(j)]
  parts <- list(all = e)
  if (any(e$r < 0) && any(e$r > 0)) parts <- c(parts, list(positive_edges = e[r > 0], negative_edges = e[r < 0]))
  rbindlist(lapply(names(parts), function(p) { x <- parts[[p]]
    as.data.table(c(list(network = l, edges = p, n = nrow(x)),
      unlist(lapply(names(A), function(ch) { s <- mean(A[[ch]][cbind(x$i, x$j)])
        setNames(c(s, s / base[[ch]]), paste0(ch, c("_supported", "_enrichment"))) }), recursive = FALSE))) }))
}))

## ---- neighbourhoods ----
nb <- lapply(E, function(e) { x <- rbind(e[, .(g = source, p = target)], e[, .(g = target, p = source)])
  split(x$p, x$g) })
common <- Reduce(intersect, lapply(nb, function(n) names(n)[lengths(n) >= MIN_DEG]))
common <- intersect(common, sg)
qs <- sample(common, min(N_GENES, length(common)))
dens <- function(set, M) { o <- unique(na.omit(idx(set))); k <- length(o); if (k < 2) NA_real_ else sum(M[o, o]) / 2 / choose(k, 2) }
nbr <- rbindlist(lapply(names(E), function(l) {
  rbindlist(lapply(qs, function(g) { s <- nb[[l]][[g]]
    as.data.table(c(list(network = l, gene = g, partners = length(s)),
      unlist(lapply(names(A), function(ch) { d <- dens(s, A[[ch]])
        nd <- mean(vapply(seq_len(N_NULL), function(b) dens(sample(sg, length(s)), A[[ch]]), numeric(1)), na.rm = TRUE)
        setNames(c(d, nd), paste0(ch, c("_density", "_null"))) }), recursive = FALSE))) }))
}))
fwrite(edges, paste0(OUT, "_edges.tsv"), sep = "\t"); fwrite(nbr, paste0(OUT, "_neighbourhoods.tsv"), sep = "\t")
cat(sprintf("%d network genes in STRING; random-pair background: independent %.4f, coexpression %.4f\n",
            length(sg), base[["independent"]], base[["coexpression"]]))
print(edges[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 4) else v)])
cat(sprintf("\nneighbourhoods: %d genes with >= %d partners in every network\n", length(qs), MIN_DEG))
# enrichment over all sampled genes: mean observed density over mean density of the size-matched random sets
print(nbr[, .(partners = median(partners),
              independent = round(mean(independent_density, na.rm = TRUE) / mean(independent_null, na.rm = TRUE), 2),
              independent_any_link = round(mean(independent_density > 0, na.rm = TRUE), 2),
              coexpression = round(mean(coexpression_density, na.rm = TRUE) / mean(coexpression_null, na.rm = TRUE), 2),
              coexpression_any_link = round(mean(coexpression_density > 0, na.rm = TRUE), 2)), by = network])
