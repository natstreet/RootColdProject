#!/usr/bin/env Rscript
# tableS3_arabidopsis.R — Table S3 (D090): GO enrichment of each core cluster, one test per cluster, on the Arabidopsis (Col-0)
# annotation. The clusters are cross-species cliques, so the Col-0 clique members stand for each family; Arabidopsis has the
# most complete annotation. Two levels, each answering its own question:
#   gene level    "which functions are over-represented among the cluster's genes": unit = gene; study = the Col-0 genes in any
#                 complete clique of the cluster's orthogroups (D082; paralogues count separately); background = the Col-0 genes
#                 of the orthogroups eligible for the core that pass the ComPlEx expression filter in Col-0 (core orthogroups:
#                 their clique genes).
#   family level  "which kinds of gene family are over-represented": unit = orthogroup; a family carries a term when any of its
#                 Col-0 genes above carries it (topGO propagation); study = the cluster's families; background = the eligible
#                 families with an annotated Col-0 gene. Also: size-matched test (the cluster's families redrawn within deciles of
#                 the size rule's expressed_genes; exact, size_matched_p()) and the majority rule (a family carries a term only
#                 when at least half of its Col-0 genes do; D087).
# Eligible orthogroups (size_rule.R, D080): at least one gene passing the ComPlEx filter in each of the five species, an
# Arabidopsis gene counting when it passes in Col-0 or Ost-0 (6,889 on version 5).
# Tests: topGO classic Fisher (BP); BH over every term of the universe's GO graph (as Tables S1 and S2), per cluster; padj < 0.05.
# Also writes the annotation-coverage table (supplementary; D090): per core family and genotype, clique genes, genes with any
# GO term, and the Jaccard index of the family's direct GO terms with those of its Col-0 clique genes.
# Inputs: data/ComPlEx/cliques/{conserved_genes_COMPLETE_filtered,conserved_genes_COMPLETE_wider}.RDS, clique_HMcluster_list.rda,
#   manuscript_cluster_map.tsv; data/ComPlEx/Orthogroups_20240613.tsv; data/expression/col0_root_cold_stress_expression.txt;
#   data/annotation/{Athal,Potra02,Betpen}_go_ids.tsv, {Picab02,Pinsy01}_go_ids.csv; results/orthogroup_size_rule.tsv
# Outputs (out dir): TableS3_gene_level.tsv, TableS3_family_level.tsv (padj < 0.05 rows; with size-matched and majority),
#   all_terms_gene_level.tsv, all_terms_family_level.tsv, Table_S3.xlsx (Legend, gene level, family level),
#   annotation_coverage.tsv, annotation_coverage_summary.tsv
# Usage (working tree): Rscript scripts/GO/tableS3_arabidopsis.R <out dir>
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
source("scripts/GO/go_table_utils.R"); source("scripts/ComPlEx/clique_matrix.R")
OUT <- commandArgs(trailingOnly = TRUE)[1]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

core <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
hm <- as.data.table(ld("data/ComPlEx/cliques/clique_HMcluster_list.rda")); mm <- fread("data/ComPlEx/cliques/manuscript_cluster_map.tsv")
core[, cluster := mm$manuscript_cluster[match(hm$Cluster[match(cliqueID, hm$ID)], mm$analysis_cluster)]]
stopifnot(!anyNA(core$cluster))
CG <- as.data.table(clique_genes_by_genotype())
sr <- fread("results/orthogroup_size_rule.tsv"); elig <- sr$OrthoGroup; stopifnot(all(core$OrthoGroup %in% elig))
og <- read.delim("data/ComPlEx/Orthogroups_20240613.tsv", check.names = FALSE, colClasses = "character"); names(og)[1] <- "OrthoGroup"
og <- og[og$OrthoGroup %in% elig, ]
x <- read.delim("data/expression/col0_root_cold_stress_expression.txt"); m <- as.matrix(x[, -1])
pass <- x[[1]][matrixStats::rowSds(m) > 0 & rowSums(m != 0) >= 3 & rowMeans(m) > 1]
s <- strsplit(og$Arabidopsis_thaliana, ", ", fixed = TRUE)
fg <- data.table(OrthoGroup = rep(og$OrthoGroup, lengths(s)), gene = unlist(s))[gene != "" & gene %in% pass]
fg <- unique(rbind(fg[!OrthoGroup %in% core$OrthoGroup], CG[genotype == "col0", .(OrthoGroup, gene)]))
stopifnot(all(CG[genotype == "col0", gene] %in% pass))                         # clique genes pass the filter they were built on
athal <- read_go_mapping("data/annotation/Athal_go_ids.tsv")

term_members <- function(unit2go) {
  gl <- factor(rep(c(0L, 1L), length.out = length(unit2go))); names(gl) <- names(unit2go)
  g <- suppressMessages(new("topGOdata", ontology = "BP", allGenes = gl, annot = annFUN.gene2GO, gene2GO = unit2go, nodeSize = 1))
  list(members = genesInTerm(g), units = genes(g))
}
fisher_terms <- function(mem, units, st) {                                        # BH over every term of the graph
  N <- length(units); n <- length(st)
  r <- rbindlist(lapply(names(mem), function(t) { mm_ <- mem[[t]]; k <- length(mm_); xx <- sum(mm_ %in% st)
    data.table(GO.ID = t, Annotated = k, Significant = xx, Expected = round(n * k / N, 2), p = phyper(xx - 1, k, N - k, n, lower.tail = FALSE)) }))
  r[, padj := p.adjust(p, "BH")][Significant > 0][order(padj, p, GO.ID)]
}
TERM <- function(ids) unname(AnnotationDbi::Term(GO.db::GOTERM[ids]))

## gene level
ga <- athal[intersect(unique(fg$gene), names(athal))]; gm <- term_members(ga)
gene2og <- setNames(fg$OrthoGroup, fg$gene)
G1 <- rbindlist(lapply(sort(unique(core$cluster)), function(k) {
  st <- intersect(CG[genotype == "col0" & OrthoGroup %in% core[cluster == k, OrthoGroup], gene], gm$units)
  r <- fisher_terms(gm$members, gm$units, st)
  r[, `:=`(cluster = k, study_genes = length(st), background_genes = length(gm$units),
           orthogroups = vapply(GO.ID, function(t) uniqueN(gene2og[intersect(gm$members[[t]], st)]), 0L))]
}))
G1[, Term := TERM(GO.ID)]

## family level
f2go <- lapply(split(fg$gene, fg$OrthoGroup), function(gs) unique(unlist(athal[intersect(gs, names(athal))], use.names = FALSE)))
f2go <- f2go[lengths(f2go) > 0]; fm <- term_members(f2go); U <- fm$units; terms <- names(fm$members)
M <- sparseMatrix(i = match(unlist(fm$members), U), j = rep(seq_along(terms), lengths(fm$members)), x = 1,
                  dims = c(length(U), length(terms)), dimnames = list(U, terms))
bins <- setNames(cut(sr$expressed_genes, unique(quantile(sr$expressed_genes, 0:10 / 10)), include.lowest = TRUE, labels = FALSE), sr$OrthoGroup)[U]
fgu <- fg[OrthoGroup %in% U]; ngen <- fgu[, .N, by = OrthoGroup]; ngen <- setNames(ngen$N, ngen$OrthoGroup)[U]
Fg <- sparseMatrix(i = match(fgu$OrthoGroup, U), j = match(fgu$gene, unique(fgu$gene)), x = 1, dims = c(length(U), uniqueN(fgu$gene)),
                   dimnames = list(U, unique(fgu$gene)))
gt <- intersect(names(gm$members), terms); gl <- lapply(gm$members[gt], intersect, colnames(Fg))
Gt <- sparseMatrix(i = match(unlist(gl), colnames(Fg)), j = rep(seq_along(gt), lengths(gl)), x = 1, dims = c(ncol(Fg), length(gt)),
                   dimnames = list(colnames(Fg), gt))
share <- as.matrix(Fg %*% Gt) / ngen
maj <- lapply(setNames(terms, terms), function(t) if (t %in% gt) U[share[, t] >= 0.5] else character(0))   # same terms (and BH m) as Fisher
F1 <- rbindlist(lapply(sort(unique(core$cluster)), function(k) {
  sf <- intersect(core[cluster == k, OrthoGroup], U)
  r <- fisher_terms(fm$members, U, sf)
  lv <- sort(unique(bins)); Ns <- as.numeric(table(factor(bins, lv))); ns <- as.numeric(table(factor(bins[sf], lv)))
  Ks <- matrix(sapply(lv, function(z) colSums(M[U[bins == z], r$GO.ID, drop = FALSE])), nrow = nrow(r))
  sm <- size_matched_p(colSums(M[sf, r$GO.ID, drop = FALSE]), Ks, Ns, ns)
  r[, `:=`(padj_size_matched = p.adjust(sm$p, "BH", n = length(terms)), null_mean = round(sm$expected, 2))]
  rmj <- fisher_terms(maj, U, sf)
  r[, padj_majority := rmj$padj[match(GO.ID, rmj$GO.ID)]][is.na(padj_majority), padj_majority := 1]
  r[, `:=`(cluster = k, study_families = length(sf), background_families = length(U),
           family_ids = vapply(GO.ID, function(t) paste(sort(intersect(fm$members[[t]], sf)), collapse = ","), ""))]
}))
F1[, Term := TERM(GO.ID)]
F1[, `:=`(size_matched = padj_size_matched < 0.05, majority = padj_majority < 0.05)]
fwrite(G1, file.path(OUT, "all_terms_gene_level.tsv"), sep = "\t"); fwrite(F1, file.path(OUT, "all_terms_family_level.tsv"), sep = "\t")
GS <- G1[padj < 0.05, .(`Cluster (manuscript)` = cluster, GO.ID, Term, Annotated, Significant, Expected, classicFisher = p, padj,
                         orthogroups, study_genes, background_genes)]
FS <- F1[padj < 0.05, .(`Cluster (manuscript)` = cluster, GO.ID, Term, Annotated, Significant, Expected, classicFisher = p, padj,
                         padj_size_matched, size_matched, padj_majority, majority, study_families, background_families, family_ids)]
GS[, Species := "Col0"]; FS[, Species := "Col0"]                                 # readers filter on cluster and GO ID
fwrite(GS, file.path(OUT, "TableS3_gene_level.tsv"), sep = "\t"); fwrite(FS, file.path(OUT, "TableS3_family_level.tsv"), sep = "\t")
write_supp_xlsx(list(`Table S3 gene level` = as.data.frame(GS), `Table S3 family level` = as.data.frame(FS)), file.path(OUT, "Table_S3.xlsx"),
  sheet = NULL, sig_cols = c("classicFisher", "padj", "padj_size_matched", "padj_majority"),
  legend = paste("Table S3 GO enrichment (biological process) of the conserved-core clusters, one test per cluster on the Arabidopsis",
    "(Col-0) annotation of the cliques. Gene level: the Col-0 genes in the cluster's complete cliques against the Col-0 genes of the",
    "orthogroups eligible for the core; asks which functions are over-represented among the cluster's genes (paralogues count",
    "separately). Family level: orthogroups, an orthogroup carrying a term when any of its Col-0 genes does, against the eligible",
    "orthogroups; asks which kinds of gene family are over-represented. topGO classic Fisher's exact test, Benjamini-Hochberg within each",
    "cluster over every term of the GO graph, padj < 0.05. Family level also: size_matched (the same test against random sets of",
    "eligible orthogroups matched on size, exact) and majority (an orthogroup carries a term only when at least half of its Col-0 genes do)."))

## annotation coverage (supplementary)
ANN <- list(col0 = athal, ost0 = athal, aspen = read_go_mapping("data/annotation/Potra02_go_ids.tsv"),
            birch = read_go_mapping("data/annotation/Betpen_go_ids.tsv"),
            spruce = unstack(read.csv("data/annotation/Picab02_go_ids.csv")[, c("go_id", "external_gene_name")]),
            pine = unstack(read.csv("data/annotation/Pinsy01_go_ids.csv")[, c("go_id", "external_gene_name")]))
cov <- CG[OrthoGroup %in% core$OrthoGroup, {
  a <- ANN[[genotype[1]]]; gs <- unique(gene); terms_ <- unique(unlist(a[intersect(gs, names(a))]))
  ref <- unique(unlist(athal[intersect(CG[genotype == "col0" & OrthoGroup == .BY$OrthoGroup, gene], names(athal))]))
  .(clique_genes = length(gs), annotated_genes = sum(gs %in% names(a)), go_terms = length(terms_),
    jaccard_with_col0 = if (length(union(terms_, ref))) round(length(intersect(terms_, ref)) / length(union(terms_, ref)), 3) else NA_real_)
}, by = .(OrthoGroup, genotype)]
cov <- merge(cov, core[, .(OrthoGroup, cluster)], by = "OrthoGroup")[order(cluster, OrthoGroup, genotype)]
fwrite(cov, file.path(OUT, "annotation_coverage.tsv"), sep = "\t")
cs <- cov[, .(families = .N, clique_genes = sum(clique_genes), annotated_share = round(sum(annotated_genes) / sum(clique_genes), 3),
              families_without_annotation = sum(annotated_genes == 0), median_jaccard_with_col0 = median(jaccard_with_col0, na.rm = TRUE)), by = genotype]
fwrite(cs, file.path(OUT, "annotation_coverage_summary.tsv"), sep = "\t"); print(cs)
cat(sprintf("gene level: %d terms padj<0.05 over %d clusters; family level: %d (size-matched %d, majority %d, both %d)\n",
    nrow(GS), uniqueN(GS$`Cluster (manuscript)`), nrow(FS), sum(FS$size_matched), sum(FS$majority), sum(FS$size_matched & FS$majority)))
