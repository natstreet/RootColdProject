#!/usr/bin/env Rscript
# pair_level_null.R — pair-level null for the conserved core (ledger R070; audit C001, C002). In each of the 15 pairwise ComPlEx
# comparisons, the co-expressologs (Max.p.val < 0.1) are replaced by a uniform random draw of the same number from that comparison's
# tested pairs (the rows of its comparison table), which is the null the Methods describe. Per permutation:
#   all15        orthogroups with at least one co-expressolog in all 15 comparisons
#   clique       orthogroups with a complete clique (one gene per genotype, all 15 pairs co-expressologs)
#   clique_rule  of those, the orthogroups the size rule counts (fitted probability < 0.5, results/orthogroup_size_rule.tsv; the
#                observed fit is kept fixed, not refitted per permutation)
#   and the same three counts per decile of expressed genes (D080 expressed_genes of the size rule; orthogroups outside the eligible
#   set fall in decile NA).
# Mode "observed" uses the real co-expressologs and must reproduce the core counts (264, 137, 130 on the current tree).
# Usage (working tree): Rscript scripts/validation/pair_level_null.R <out tsv> <permutations, or "observed"> [seed]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); OUT <- a[1]; MODE <- a[2]; SEED <- if (length(a) >= 3) as.integer(a[3]) else 1L
set.seed(SEED)
GENO <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
files <- list.files("data/ComPlEx/PairwiseComparisons", "^comparison-.*-abspearsonMR0\\.03-table\\.RData$", full.names = TRUE)
stopifnot(length(files) == 15)
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
cmp <- lapply(files, function(f) {
  g <- strsplit(sub("^comparison-(.*)-abspearsonMR0\\.03-table\\.RData$", "\\1", basename(f)), "-")[[1]]
  x <- as.data.table(ld(f))[, .(g1 = paste0(g[1], ":", Species1), g2 = paste0(g[2], ":", Species2), og = OrthoGroup,
                                co = as.numeric(Max.p.val) < 0.1)]
  x[, `:=`(p1 = g[1], p2 = g[2])]; x })
stopifnot(setequal(unique(unlist(lapply(cmp, function(x) c(x$p1[1], x$p2[1])))), GENO))
sr <- fread("results/orthogroup_size_rule.tsv")
dec <- setNames(cut(sr$expressed_genes, unique(quantile(sr$expressed_genes, 0:10 / 10)), include.lowest = TRUE, labels = FALSE), sr$OrthoGroup)
rule_ok <- setNames(sr$p_core_given_size < 0.5, sr$OrthoGroup)

# complete clique: one gene per genotype, every pair an edge. Genes without a neighbour in some other genotype are pruned
# repeatedly; then backtracking with forward checking (choosing a gene restricts every other genotype to its neighbours).
has_clique <- function(edges) {                      # edges: data.table g1, g2 (genotype-tagged genes)
  nb <- split(c(edges$g2, edges$g1), c(edges$g1, edges$g2)); nb <- lapply(nb, unique)
  genes <- names(nb); cand <- split(genes, sub(":.*$", "", genes))
  if (!all(GENO %in% names(cand))) return(FALSE)
  cand <- cand[GENO]
  repeat {                                           # prune genes lacking a neighbour in any other genotype's candidates
    n0 <- sum(lengths(cand))
    for (h in GENO) cand[[h]] <- cand[[h]][vapply(cand[[h]], function(g) all(vapply(setdiff(GENO, h), function(k)
      any(nb[[g]] %in% cand[[k]]), TRUE)), TRUE)]
    if (any(lengths(cand) == 0)) return(FALSE)
    if (sum(lengths(cand)) == n0) break }
  rec <- function(cand) {
    if (!length(cand)) return(TRUE)
    h <- names(cand)[which.min(lengths(cand))]
    for (g in cand[[h]]) {
      nc <- lapply(cand[names(cand) != h], function(v) v[v %in% nb[[g]]])
      if (all(lengths(nc) > 0) && rec(nc)) return(TRUE) }
    FALSE }
  rec(cand)
}
# known-answer test: a 6-clique across the genotypes is found; removing one of its 15 edges leaves none
toy <- CJ(i = 1:6, j = 1:6)[i < j][, .(g1 = paste0(GENO[i], ":x"), g2 = paste0(GENO[j], ":x"))]
toy <- rbind(toy, data.table(g1 = "col0:y", g2 = "ost0:x"))
stopifnot(has_clique(toy), !has_clique(toy[-3]))

count_once <- function(sel) {                        # sel: list of 15 data.tables of co-expressolog pairs (g1, g2, og)
  per <- rbindlist(lapply(sel, function(x) data.table(og = unique(x$og))))[, .N, by = og]
  a15 <- per[N == 15, og]
  pairs <- rbindlist(lapply(sel, function(x) x[og %in% a15, .(g1, g2, og)]))
  cl <- if (length(a15)) a15[vapply(split(pairs, by = "og")[a15], has_clique, TRUE)] else character(0)
  cr <- cl[!is.na(rule_ok[cl]) & rule_ok[cl]]
  d <- function(v) table(factor(dec[v], 1:max(dec, na.rm = TRUE)), useNA = "always")
  list(tot = c(all15 = length(a15), clique = length(cl), clique_rule = length(cr)),
       dec = data.table(decile = names(d(a15)), all15 = as.integer(d(a15)), clique = as.integer(d(cl)), clique_rule = as.integer(d(cr))))
}
if (MODE == "observed") {
  r <- count_once(lapply(cmp, function(x) x[co == TRUE]))
  print(r$tot); fwrite(cbind(perm = 0L, r$dec), OUT, sep = "\t"); quit(save = "no")
}
NP <- as.integer(MODE); res <- vector("list", NP)
for (p in seq_len(NP)) {
  sel <- lapply(cmp, function(x) x[sample.int(nrow(x), sum(x$co))])
  r <- count_once(sel); res[[p]] <- cbind(perm = p, seed = SEED, r$dec)
  if (p %% 25 == 0) message(p, " permutations: last ", paste(names(r$tot), r$tot, collapse = ", "))
}
fwrite(rbindlist(res), OUT, sep = "\t")
