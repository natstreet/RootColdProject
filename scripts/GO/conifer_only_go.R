#!/usr/bin/env Rscript
# conifer_only_go.R — GO of the conifer-only co-expressologs WITHOUT the cold-DEG restriction of Table S2 (ledger R064, D089).
#
# Sets (as make_TableS2.R, co-expressolog threshold Max.p.val < THR, default 0.1): A = orthogroups with a Col-0/Ost-0
# co-expressolog, B = aspen/birch, C = spruce/pine; the conifer-only set is C minus (A or B). Study genes per conifer: that
# conifer's genes in the spruce-pine co-expressolog pairs of conifer-only orthogroups (Table S2 then keeps only cold DEGs;
# this script does not).
# 1. Gene level, per conifer: topGO classic Fisher, BP, background = that conifer's expressed genes (data/annotation/bg_*),
#    BH over all tested terms. For each term with padj < 0.05: annotated study genes, the distinct orthogroups they come
#    from, and the share from the largest orthogroup.
# 2. Family level, per conifer: unit = orthogroup; background = the orthogroups tested in the spruce-pine comparison
#    (data/ComPlEx/PairwiseComparisons/comparison-spruce-pine-...-table.RData), so every background family could have become
#    a conifer co-expressolog. A family carries a term when any of its genes carries it after topGO propagation: its
#    co-expressolog genes for the conifer-only families, its tested genes (genes in the comparison table) otherwise.
#    One-sided Fisher, BH over every term of the universe's GO graph (as topGO); and a size-matched null: the conifer-only families
#    redrawn within deciles of the number of tested genes, P computed exactly (size_matched_p(), go_table_utils.R), BH.
# Outputs (out dir): gene_level_<conifer>.tsv (every tested term), gene_level_family_count_<conifer>.tsv (significant terms),
#   gene_level_study_genes_<conifer>.tsv (the annotated study genes of each significant term, with their orthogroup),
#   family_level_<conifer>.tsv (every tested term, Fisher and size-matched), summary.txt
# Usage (working tree): Rscript scripts/GO/conifer_only_go.R <out dir> [threshold, default 0.1] [null draws, default 5000]
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
source("scripts/GO/go_table_utils.R")
a <- commandArgs(trailingOnly = TRUE); OUT <- a[1]; THR <- if (length(a) >= 2) as.numeric(a[2]) else 0.1
NULL_DRAWS <- if (length(a) >= 3) as.integer(a[3]) else 5000L
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); set.seed(1)
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
load("data/ComPlEx/co_expressologs.RData")
co <- as.data.table(co_expressologs)[MaxpVal < THR]
co[, pairkey := paste(pmin(Species1, Species2), pmax(Species1, Species2), sep = "|")]
A <- co[pairkey == "col0|ost0", unique(OrthologGroup)]; B <- co[pairkey == "aspen|birch", unique(OrthologGroup)]
C <- co[pairkey == "pine|spruce", unique(OrthologGroup)]; Conly <- setdiff(C, union(A, B))
cp <- co[pairkey == "pine|spruce" & OrthologGroup %in% Conly]
study <- list(spruce = unique(data.table(gene = c(cp[Species1 == "spruce", GeneSpecies1], cp[Species2 == "spruce", GeneSpecies2]),
                                         OrthoGroup = c(cp[Species1 == "spruce", OrthologGroup], cp[Species2 == "spruce", OrthologGroup]))),
              pine   = unique(data.table(gene = c(cp[Species1 == "pine", GeneSpecies1], cp[Species2 == "pine", GeneSpecies2]),
                                         OrthoGroup = c(cp[Species1 == "pine", OrthologGroup], cp[Species2 == "pine", OrthologGroup]))))
ANNf <- list(spruce = unstack(read.csv("data/annotation/Picab02_go_ids.csv")[, c("go_id", "external_gene_name")]),
             pine   = unstack(read.csv("data/annotation/Pinsy01_go_ids.csv")[, c("go_id", "external_gene_name")]))
BG <- list(spruce = ld("data/annotation/bg_spruce_PA.rda"), pine = ld("data/annotation/bg_pine_PS.rda"))
tab <- as.data.table(ld("data/ComPlEx/PairwiseComparisons/comparison-spruce-pine-abspearsonMR0.03-table.RData"))
tested <- list(spruce = unique(tab[, .(OrthoGroup, gene = Species1)]), pine = unique(tab[, .(OrthoGroup, gene = Species2)]))
term_members <- function(unit2go) {
  gl <- factor(rep(c(0L, 1L), length.out = length(unit2go))); names(gl) <- names(unit2go)
  g <- suppressMessages(new("topGOdata", ontology = "BP", allGenes = gl, annot = annFUN.gene2GO, gene2GO = unit2go, nodeSize = 1))
  list(members = genesInTerm(g), units = genes(g))
}
fisher_terms <- function(mem, units, st) {
  N <- length(units); n <- length(st)
  r <- rbindlist(lapply(names(mem), function(t) { m <- mem[[t]]; k <- length(m); x <- sum(m %in% st)
    data.table(GO.ID = t, Annotated = k, Significant = x, Expected = round(n * k / N, 2), p = phyper(x - 1, k, N - k, n, lower.tail = FALSE)) }))
  r[, padj := p.adjust(p, "BH")][Significant > 0][order(padj, p, GO.ID)]   # BH over every term of the graph, as topGO / Tables S1-S3
}
lines <- c(sprintf("conifer-only orthogroups (C minus A or B, Max.p.val < %g): %d", THR, length(Conly)))
for (cf in c("spruce", "pine")) {
  ## 1. gene level
  ann <- ANNf[[cf]][names(ANNf[[cf]]) %in% BG[[cf]]]; tm <- term_members(ann)
  st <- intersect(study[[cf]]$gene, tm$units)
  g1 <- fisher_terms(tm$members, tm$units, st); g1[, Term := unname(AnnotationDbi::Term(GO.db::GOTERM[GO.ID]))]
  fwrite(g1, file.path(OUT, sprintf("gene_level_%s.tsv", cf)), sep = "\t")
  fc <- rbindlist(lapply(g1[padj < 0.05, GO.ID], function(t) { h <- study[[cf]][gene %in% intersect(tm$members[[t]], st)]
    data.table(GO.ID = t, study_genes_annotated = uniqueN(h$gene), orthogroups = uniqueN(h$OrthoGroup),
               share_largest_orthogroup = round(max(table(unique(h)$OrthoGroup)) / uniqueN(h$gene), 3)) }))
  if (nrow(fc)) fc <- merge(g1[padj < 0.05, .(GO.ID, Term, padj)], fc, by = "GO.ID")[order(padj)]
  fwrite(fc, file.path(OUT, sprintf("gene_level_family_count_%s.tsv", cf)), sep = "\t")
  sg <- rbindlist(lapply(g1[padj < 0.05, GO.ID], function(t) study[[cf]][gene %in% intersect(tm$members[[t]], st)][, GO.ID := t]))
  fwrite(sg[order(GO.ID, OrthoGroup, gene), .(GO.ID, gene, OrthoGroup)], file.path(OUT, sprintf("gene_level_study_genes_%s.tsv", cf)), sep = "\t")
  ## 2. family level
  fg <- rbind(tested[[cf]][!OrthoGroup %in% Conly], study[[cf]][, .(OrthoGroup, gene)])
  f2go <- lapply(split(fg$gene, fg$OrthoGroup), function(gs) unique(unlist(ANNf[[cf]][intersect(gs, names(ANNf[[cf]]))], use.names = FALSE)))
  f2go <- f2go[lengths(f2go) > 0]; fm <- term_members(f2go); U <- fm$units; sf <- intersect(Conly, U)
  f1 <- fisher_terms(fm$members, U, sf)
  terms <- f1$GO.ID; M <- sparseMatrix(i = match(unlist(fm$members[terms]), U), j = rep(seq_along(terms), lengths(fm$members[terms])),
                                       x = 1, dims = c(length(U), length(terms)), dimnames = list(U, terms))
  sz <- fg[OrthoGroup %in% U, .N, by = OrthoGroup]; sz <- setNames(sz$N, sz$OrthoGroup)[U]
  bins <- setNames(cut(sz, unique(quantile(sz, 0:10 / 10)), include.lowest = TRUE, labels = FALSE), U)
  lv <- sort(unique(bins)); obs <- colSums(M[sf, , drop = FALSE])
  Ns <- as.numeric(table(factor(bins, lv))); ns <- as.numeric(table(factor(bins[sf], lv)))
  Ks <- matrix(sapply(lv, function(z) colSums(M[U[bins == z], , drop = FALSE])), nrow = length(terms))
  sm <- size_matched_p(obs, Ks, Ns, ns)                                       # exact (go_table_utils.R)
  f1[, `:=`(p_size_matched = sm$p, null_mean = round(sm$expected, 2))]
  f1[, padj_size_matched := p.adjust(p_size_matched, "BH", n = length(fm$members))]          # same m as Fisher
  f1[, Term := unname(AnnotationDbi::Term(GO.db::GOTERM[GO.ID]))]
  f1[, family_ids := vapply(GO.ID, function(t) paste(sort(intersect(fm$members[[t]], sf)), collapse = ","), "")]
  fwrite(f1, file.path(OUT, sprintf("family_level_%s.tsv", cf)), sep = "\t")
  lines <- c(lines, sprintf("%s: study genes %d (annotated %d) in %d conifer-only orthogroups; gene level: %d terms padj<0.05 (min padj %.2g); family level: %d of the %d annotated tested families are conifer-only (of %d conifer-only); Fisher padj<0.05 %d, size-matched %d, both %d",
    cf, nrow(unique(study[[cf]][, .(gene)])), length(st), uniqueN(study[[cf]]$OrthoGroup), g1[padj < 0.05, .N], min(g1$padj),
    length(sf), length(U), length(Conly), f1[padj < 0.05, .N], f1[padj_size_matched < 0.05, .N], f1[padj < 0.05 & padj_size_matched < 0.05, .N]))
}
writeLines(lines, file.path(OUT, "summary.txt")); cat(lines, sep = "\n")
