#!/usr/bin/env Rscript
# conserved_neighbourhood_string.R — functional coherence of the conserved neighbourhoods behind one species pair's
# co-expressologs, against STRING v12 (A. thaliana), for one correlation-sign variant of ComPlEx (abs, pos or neg).
#
# For a random sample of the variant's co-expressologs (Max.p.val < 0.1), the conserved neighbourhood of the
# species-1 member is rebuilt exactly as ComPlEx.R tests it: its neighbours (mutual rank at the 3% density threshold,
# the gene itself excluded) that are orthologues of the species-2 member's neighbours. Its genes are mapped to their
# orthogroups, and each orthogroup to its A. thaliana members (data/orthologs/orthologs.parquet). Two orthogroups are
# linked when any pair of their A. thaliana members is linked in STRING at score >= 400 in a channel:
#   independent:  max(experimental, database), evidence that does not come from expression
#   coexpression: STRING's co-expression channel (an external A. thaliana compendium, as in Fig. S4b)
# Coherence = the fraction of orthogroup pairs in the neighbourhood that are linked; enrichment = that fraction over
# its mean in random sets of the same number of genes drawn from the species-1 genes ComPlEx tested (N_NULL sets).
#
# Usage (from a directory with data/ as for ComPlEx.R):
#   Rscript conserved_neighbourhood_string.R <sp1> <sp2> <abs|pos|neg> <comparison table> <STRING dir> <out tsv> [N=500]
suppressMessages({ library(data.table); library(matrixStats); library(Matrix); library(arrow) })
a <- commandArgs(trailingOnly = TRUE)
s1 <- a[1]; s2 <- a[2]; SIGN <- a[3]; TAB <- a[4]; SD <- a[5]; OUT <- a[6]; N <- if (length(a) >= 7) as.integer(a[7]) else 500L
N_NULL <- 20L; THR <- 400; P <- 0.1; density_thr <- 0.03
set.seed(1)
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

## ---- ComPlEx's gene universe and networks for this pair ----
ortho <- as.data.table(ld(sprintf("data/ComPlEx/orthologs-%s-%s.RData", s1, s2)))
expressed <- function(sp) {
  x <- read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", sp))
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
  m <- m[rowSds(m) > 0 & rowCounts(m != 0) >= 3, , drop = FALSE]
  m[rowMeans(m) > 1.0, , drop = FALSE]
}
e1 <- expressed(s1); e2 <- expressed(s2)
ortho <- ortho[Species1 %in% rownames(e1) & Species2 %in% rownames(e2)]
e1 <- e1[rownames(e1) %in% ortho$Species1, , drop = FALSE]; e2 <- e2[rownames(e2) %in% ortho$Species2, , drop = FALSE]
neighbours <- function(e, genes) {
  net <- cor(t(e)); if (SIGN == "neg") net <- -net; if (SIGN == "abs") net <- abs(net)
  R <- t(apply(net, 1, rank)); rm(net); gc()
  net <- sqrt(R * t(R)); rm(R); gc(); diag(net) <- 0
  v <- sort(net[upper.tri(net)], decreasing = TRUE); thr <- v[round(density_thr * length(v))]; rm(v); gc()
  out <- lapply(genes, function(g) { r <- net[g, ]; names(r)[r >= thr] }); names(out) <- genes; rm(net); gc(); out
}
tab <- as.data.table(ld(TAB)); tab <- tab[as.numeric(Max.p.val) < P]
smp <- tab[sample(.N, min(N, .N))]
nb1 <- neighbours(e1, unique(smp$Species1)); nb2 <- neighbours(e2, unique(smp$Species2))

## ---- orthogroups -> A. thaliana -> STRING links between orthogroups ----
og1 <- unique(ortho[, .(gene = Species1, og = OrthoGroup)])
pq <- as.data.table(read_parquet("data/orthologs/orthologs.parquet"))
ath <- unique(pq[species == "Arabidopsis_thaliana", .(og = Ortholog_Group, locus = toupper(trimws(as.character(gene))))])
ath <- ath[grepl("^AT[1-5MC]G[0-9]{5}$", locus) & og %in% og1$og]
al <- fread(cmd = paste("gzip -dc", shQuote(file.path(SD, "3702.protein.aliases.v12.0.txt.gz"))), sep = "\t", header = FALSE,
            skip = 1, select = 1:2, col.names = c("prot", "alias"), quote = "")
al[, alias := toupper(alias)]; al <- al[grepl("^AT[1-5MC]G[0-9]{5}$", alias)][!duplicated(alias)]
ath[, prot := al$prot[match(locus, al$alias)]]; ath <- ath[!is.na(prot)]
lk <- fread(cmd = paste("gzip -dc", shQuote(file.path(SD, "3702.protein.links.detailed.v12.0.txt.gz"))), sep = " ",
            select = c("protein1", "protein2", "coexpression", "experimental", "database"))
lk <- lk[protein1 %in% ath$prot & protein2 %in% ath$prot]
ogs <- sort(unique(ath$og))
# a protein can sit in only one orthogroup here (first match); orthogroup adjacency from protein links
p2og <- ath[!duplicated(prot), setNames(match(og, ogs), prot)]
adj <- function(keep) { l <- lk[keep]; i <- p2og[l$protein1]; j <- p2og[l$protein2]; ok <- i != j
  m <- sparseMatrix(i = c(i[ok], j[ok]), j = c(j[ok], i[ok]), x = 1, dims = c(length(ogs), length(ogs))); m@x[] <- 1; m }
A <- list(independent = adj(pmax(lk$experimental, lk$database) >= THR), coexpression = adj(lk$coexpression >= THR))
gene_og <- setNames(match(og1$og, ogs), og1$gene)                        # NA: orthogroup without A. thaliana in STRING
dens <- function(genes, M) { o <- unique(na.omit(gene_og[genes])); k <- length(o)
  if (k < 2) return(c(NA_real_, k)); c(sum(M[o, o]) / 2 / choose(k, 2), k) }

## ---- per sampled co-expressolog ----
univ <- rownames(e1)
res <- rbindlist(lapply(seq_len(nrow(smp)), function(r) {
  g1 <- smp$Species1[r]; g2 <- smp$Species2[r]
  ov <- unique(intersect(nb1[[g1]], ortho$Species1[ortho$Species2 %in% nb2[[g2]]]))
  if (length(ov) != smp$Species1.neigh.overlap[r]) stop("overlap of ", g1, "-", g2, " differs from the comparison table")
  out <- list(Species1 = g1, Species2 = g2, overlap = length(ov))
  for (ch in names(A)) {
    d <- dens(ov, A[[ch]])
    nd <- vapply(seq_len(N_NULL), function(b) dens(sample(univ, length(ov)), A[[ch]])[1], numeric(1))
    out[[paste0(ch, "_ogs")]] <- d[2]; out[[paste0(ch, "_density")]] <- d[1]
    out[[paste0(ch, "_null")]] <- mean(nd, na.rm = TRUE); out[[paste0(ch, "_enrichment")]] <- d[1] / mean(nd, na.rm = TRUE)
  }
  as.data.table(out)
}))
res[, `:=`(sign = SIGN, pair = paste(s1, s2, sep = "-"))]
fwrite(res, OUT, sep = "\t")
cat(sprintf("%s %s-%s: %d co-expressologs (of %d) sampled; median conserved neighbourhood %d genes\n", SIGN, s1, s2,
            nrow(res), nrow(tab), as.integer(median(res$overlap))))
for (ch in names(A)) cat(sprintf("  %-12s median enrichment %.2f; enriched (>1) in %.0f%%\n", ch,
  median(res[[paste0(ch, "_enrichment")]], na.rm = TRUE), 100 * mean(res[[paste0(ch, "_enrichment")]] > 1, na.rm = TRUE)))
