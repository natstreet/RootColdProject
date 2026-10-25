#!/usr/bin/env Rscript
# tableS3_family_level.R — GO enrichment of the core clusters with the orthogroup (family) as the unit (D086, ledger R062).
#
# Step 1, family count of the gene-level Table S3 (every clique gene, make_TableS3.R): for each significant term of each
#   cluster x genotype, the study genes annotated to it (topGO propagation, as in that table), the number of distinct
#   orthogroups they come from and the share from the largest one; terms resting on one or two orthogroups are flagged.
# Step 2, family-level test, per genotype. Unit = orthogroup. Background = the orthogroups eligible for the core (the size
#   rule's eligible orthogroups, D080; results/orthogroup_size_rule.tsv). A family carries a GO term in a genotype when any
#   of its genes there carries it after topGO's propagation (biological process, GO.db): for a core family its clique
#   genes in that genotype (D082), for any other family its genes in that genotype that pass the ComPlEx expression
#   filter (size_rule.R). Gene annotation: the files of make_TableS3.R. Families with no annotated gene are not tested
#   (as topGO leaves out unannotated genes). Study set = the families of one cluster. One-sided Fisher's exact test, BH
#   over every term of the universe's GO graph (as topGO and Tables S1-S3), per cluster x genotype; padj < 0.05.
# Majority rule (D087): the same Fisher test with a stricter family annotation: a family carries a term in a genotype only when
#   at least half of its genes there (the same genes as above, annotated or not) carry it after propagation. Same universe.
# Step 3, size-matched null: per cluster x genotype, the cluster's families redrawn at random from the annotated background
#   within size strata (deciles of the size rule's expressed_genes, without replacement); P of at least the observed number
#   of term-carrying families, computed exactly as a sum of hypergeometric counts (size_matched_p() in go_table_utils.R;
#   no Monte Carlo floor), BH with the same number of terms as the Fisher test (terms without a study family have P = 1). Reported for the terms tested in step 2, with whether the two tests agree at padj < 0.05.
#
# Inputs (working tree): data/ComPlEx/cliques/{conserved_genes_COMPLETE_filtered,conserved_genes_COMPLETE_wider}.RDS,
#   clique_HMcluster_list.rda, manuscript_cluster_map.tsv; data/ComPlEx/Orthogroups_20240613.tsv;
#   data/expression/<genotype>_root_cold_stress_expression.txt; data/annotation/*go_ids*; results/orthogroup_size_rule.tsv;
#   <gene-level Table S3 tsv> (first argument).
# Outputs (second argument, a directory): step1_family_count.tsv, step1_cell_summary.tsv, TableS3_family_level.tsv,
#   Table_S3_family_level.xlsx (legend sheet), step3_size_matched.tsv, majority_rule.tsv (every tested term, both rules)
# Usage: Rscript scripts/GO/tableS3_family_level.R <gene-level TableS3 tsv> <out dir> [ignored; kept for the driver's call]
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
source("scripts/GO/go_table_utils.R"); source("scripts/ComPlEx/clique_matrix.R")
a <- commandArgs(trailingOnly = TRUE); S3 <- a[1]; OUT <- a[2]; NULL_DRAWS <- if (length(a) >= 3) as.integer(a[3]) else 5000L
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); set.seed(1)
GT <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
LAB <- c(col0 = "Col0", ost0 = "Ost0", aspen = "Aspen", birch = "Birch", spruce = "Spruce", pine = "Pine")
OGCOL <- c(col0 = "Arabidopsis_thaliana", ost0 = "Arabidopsis_thaliana", aspen = "Populus_tremula", birch = "Betula_pendula",
           spruce = "Picea_abies", pine = "Pinus_sylvestris")
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

## gene annotation (make_TableS3.R files; not restricted to a background here: families decide membership)
A <- "data/annotation"
ann_raw <- list(arab = read_go_mapping(file.path(A, "Athal_go_ids.tsv")), aspen = read_go_mapping(file.path(A, "Potra02_go_ids.tsv")),
  birch = read_go_mapping(file.path(A, "Betpen_go_ids.tsv")),
  spruce = unstack(read.csv(file.path(A, "Picab02_go_ids.csv"))[, c("go_id", "external_gene_name")]),
  pine = unstack(read.csv(file.path(A, "Pinsy01_go_ids.csv"))[, c("go_id", "external_gene_name")]))
ANN <- list(col0 = ann_raw$arab, ost0 = ann_raw$arab, aspen = ann_raw$aspen, birch = ann_raw$birch, spruce = ann_raw$spruce, pine = ann_raw$pine)

## clusters and clique genes (D082)
core <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
hm <- as.data.table(ld("data/ComPlEx/cliques/clique_HMcluster_list.rda"))
mm <- fread("data/ComPlEx/cliques/manuscript_cluster_map.tsv")
core[, cluster := mm$manuscript_cluster[match(hm$Cluster[match(cliqueID, hm$ID)], mm$analysis_cluster)]]
stopifnot(!anyNA(core$cluster))
CG <- as.data.table(clique_genes_by_genotype())                                    # OrthoGroup, genotype, gene

## eligible families and their genes per genotype
sr <- fread("results/orthogroup_size_rule.tsv"); elig <- sr$OrthoGroup
stopifnot(all(core$OrthoGroup %in% elig))
og <- read.delim("data/ComPlEx/Orthogroups_20240613.tsv", check.names = FALSE, colClasses = "character"); names(og)[1] <- "OrthoGroup"
og <- og[og$OrthoGroup %in% elig, ]
passes <- function(g) { x <- read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", g)); m <- as.matrix(x[, -1])
  x[[1]][matrixStats::rowSds(m) > 0 & rowSums(m != 0) >= 3 & rowMeans(m) > 1] }
fam_genes <- function(g) {
  p <- passes(g); s <- strsplit(og[[OGCOL[[g]]]], ", ", fixed = TRUE)
  d <- data.table(OrthoGroup = rep(og$OrthoGroup, lengths(s)), gene = unlist(s))[gene != "" & gene %in% p]
  rbind(d[!OrthoGroup %in% core$OrthoGroup], CG[genotype == g, .(OrthoGroup, gene)])  # core families: clique genes
}

## term membership with topGO's propagation: one topGOdata per genotype over families (or genes)
term_members <- function(unit2go) {
  gl <- factor(rep(c(0L, 1L), length.out = length(unit2go))); names(gl) <- names(unit2go)
  g <- suppressMessages(new("topGOdata", ontology = "BP", allGenes = gl, annot = annFUN.gene2GO, gene2GO = unit2go, nodeSize = 1))
  list(members = genesInTerm(g), units = genes(g))
}
fisher_terms <- function(mem, units, study) {
  st <- intersect(study, units); N <- length(units); n <- length(st)
  r <- rbindlist(lapply(names(mem), function(t) { m <- mem[[t]]; k <- length(m); x <- sum(m %in% st)
    data.table(GO.ID = t, Annotated = k, Significant = x, Expected = round(n * k / N, 2),
               p = phyper(x - 1, k, N - k, n, lower.tail = FALSE)) }))
  r[, padj := p.adjust(p, "BH")][Significant > 0][order(padj, p, GO.ID)]   # BH over every term of the graph, as topGO / Tables S1-S3
}

## step 1
s3 <- fread(S3); s1 <- list()
for (g in GT) {
  bgg <- ld(file.path(A, sprintf("bg_%s.rda", c(col0 = "col0_ATC", ost0 = "ost0_ATO", aspen = "aspen_PT", birch = "birch_BP",
                                                 spruce = "spruce_PA", pine = "pine_PS")[[g]])))
  tm <- term_members(ANN[[g]][names(ANN[[g]]) %in% bgg])
  for (k in sort(unique(core$cluster))) {
    st <- CG[genotype == g & OrthoGroup %in% core[cluster == k, OrthoGroup]]
    terms <- s3[`Cluster (manuscript)` == k & Species == LAB[[g]], GO.ID]
    for (t in terms) { hit <- st[gene %in% tm$members[[t]]]
      n_og <- uniqueN(hit$OrthoGroup); top <- if (nrow(hit)) max(table(hit$OrthoGroup)) / nrow(hit) else NA_real_
      s1[[length(s1) + 1]] <- data.table(cluster = k, genotype = g, GO.ID = t, study_genes_annotated = nrow(hit),
        orthogroups = n_og, share_largest_orthogroup = round(top, 3), rests_on_1_or_2 = n_og <= 2) }
  }
}
s1 <- rbindlist(s1); s1[, Term := s3$Term[match(GO.ID, s3$GO.ID)]]
fwrite(s1, file.path(OUT, "step1_family_count.tsv"), sep = "\t")
cell <- s1[, .(significant_terms = .N, terms_on_3plus_orthogroups = sum(orthogroups >= 3)), by = .(cluster, genotype)][order(cluster, genotype)]
fwrite(cell, file.path(OUT, "step1_cell_summary.tsv"), sep = "\t")
cat("step 1:", nrow(s1), "significant gene-level terms;", sum(s1$rests_on_1_or_2), "rest on 1 or 2 orthogroups\n")

## steps 2 and 3
bins <- setNames(cut(sr$expressed_genes, unique(quantile(sr$expressed_genes, 0:10 / 10)), include.lowest = TRUE, labels = FALSE), sr$OrthoGroup)
res2 <- list(); res3 <- list(); resm <- list()
for (g in GT) {
  fg <- fam_genes(g); a2 <- ANN[[g]]
  f2go <- lapply(split(fg$gene, fg$OrthoGroup), function(gs) unique(unlist(a2[intersect(gs, names(a2))], use.names = FALSE)))
  f2go <- f2go[lengths(f2go) > 0]
  tm <- term_members(f2go); U <- tm$units
  terms <- names(tm$members); M <- sparseMatrix(i = match(unlist(tm$members), U), j = rep(seq_along(terms), lengths(tm$members)),
                                                x = 1, dims = c(length(U), length(terms)), dimnames = list(U, terms))
  ## majority rule: gene-level propagation over the families' annotated genes, then share of each family's genes per term
  ga <- a2[intersect(unique(fg$gene), names(a2))]; ga <- ga[lengths(ga) > 0]
  gm <- term_members(ga); fgu <- unique(fg[OrthoGroup %in% U])
  ngen <- fgu[, .N, by = OrthoGroup]; ngen <- setNames(ngen$N, ngen$OrthoGroup)[U]
  Fg <- sparseMatrix(i = match(fgu$OrthoGroup, U), j = match(fgu$gene, unique(fgu$gene)), x = 1,
                     dims = c(length(U), uniqueN(fgu$gene)), dimnames = list(U, unique(fgu$gene)))
  gt <- intersect(names(gm$members), terms); gl <- lapply(gm$members[gt], intersect, colnames(Fg))
  Gt <- sparseMatrix(i = match(unlist(gl), colnames(Fg)), j = rep(seq_along(gt), lengths(gl)), x = 1,
                     dims = c(ncol(Fg), length(gt)), dimnames = list(colnames(Fg), gt))
  share <- as.matrix(Fg %*% Gt) / ngen
  maj <- lapply(setNames(gt, gt), function(t) U[share[, t] >= 0.5])
  for (k in sort(unique(core$cluster))) {
    study <- core[cluster == k, OrthoGroup]; st <- intersect(study, U)
    rm_ <- fisher_terms(maj, U, st)
    resm[[paste(g, k)]] <- rm_[, .(cluster = k, genotype = g, GO.ID, Significant_majority = Significant, Annotated_majority = Annotated,
                                   p_majority = p, padj_majority = padj)]
    r <- fisher_terms(tm$members, U, st)
    r[, `:=`(cluster = k, genotype = g, families_study = length(st), families_background = length(U),
             family_ids = vapply(GO.ID, function(t) paste(sort(intersect(tm$members[[t]], st)), collapse = ","), ""))]
    res2[[paste(g, k)]] <- r
    ## size-matched null
    obs <- colSums(M[st, r$GO.ID, drop = FALSE]); b <- bins[U]; lv <- sort(unique(b))
    Ns <- as.numeric(table(factor(b, lv))); ns <- as.numeric(table(factor(bins[st], lv)))
    Ks <- sapply(lv, function(z) colSums(M[U[b == z], r$GO.ID, drop = FALSE]))
    Ks <- matrix(Ks, nrow = nrow(r))                                         # terms x size strata
    sm <- size_matched_p(obs, Ks, Ns, ns)                                    # exact (go_table_utils.R)
    res3[[paste(g, k)]] <- data.table(cluster = k, genotype = g, GO.ID = r$GO.ID, observed = obs, null_mean = round(sm$expected, 2),
                                      p_size_matched = sm$p, padj_size_matched = p.adjust(sm$p, "BH", n = length(tm$members)), padj_fisher = r$padj)
    cat(sprintf("  %s cluster %d: %d study families (of %d annotated); terms tested %d; Fisher padj<0.05 %d; size-matched padj<0.05 %d\n",
                g, k, length(st), length(U), nrow(r), sum(r$padj < 0.05), sum(res3[[paste(g, k)]]$padj_size_matched < 0.05)))
  }
}
R2 <- rbindlist(res2); R2[, Term := unname(AnnotationDbi::Term(GO.db::GOTERM[GO.ID]))]
R3 <- rbindlist(res3); R3[, Term := R2$Term[match(GO.ID, R2$GO.ID)]]
R3[, agree := (padj_fisher < 0.05) == (padj_size_matched < 0.05)]
RM <- rbindlist(resm)
MJ <- merge(R2[, .(cluster, genotype, GO.ID, Term, padj_fisher = padj)], RM, by = c("cluster", "genotype", "GO.ID"), all.x = TRUE)
MJ[is.na(padj_majority), `:=`(Significant_majority = 0L, padj_majority = 1)]
fwrite(MJ[order(cluster, genotype, padj_fisher)], file.path(OUT, "majority_rule.tsv"), sep = "\t")
sig <- R2[padj < 0.05, .(`Cluster (manuscript)` = cluster, Species = LAB[genotype], genotype, GO.ID, Term, Annotated, Significant, Expected,
                         classicFisher = p, padj, families_study, families_background, family_ids)]
sig <- merge(sig, R3[, .(`Cluster (manuscript)` = cluster, genotype, GO.ID, padj_size_matched)],
             by = c("Cluster (manuscript)", "genotype", "GO.ID"), all.x = TRUE)
sig <- merge(sig, MJ[, .(`Cluster (manuscript)` = cluster, genotype, GO.ID, padj_majority)], by = c("Cluster (manuscript)", "genotype", "GO.ID"), all.x = TRUE)
sig[, call := fifelse(padj_size_matched < 0.05, "both", "Fisher only")]                       # D087
sig[, majority := padj_majority < 0.05]
sig <- sig[order(`Cluster (manuscript)`, Species, -(call == "both"), padj)][, genotype := NULL]
fwrite(sig, file.path(OUT, "TableS3_family_level.tsv"), sep = "\t")
fwrite(R3, file.path(OUT, "step3_size_matched.tsv"), sep = "\t")
write_supp_xlsx(as.data.frame(sig), file.path(OUT, "Table_S3_family_level.xlsx"), "Table S3", sig_cols = c("classicFisher", "padj", "padj_size_matched"),
  legend = paste("Table S3 GO enrichment of the conserved-core clusters with the orthogroup as the unit, per species (provisional",
    "cluster numbering). An orthogroup carries a term in a species when any of its genes there carries it (topGO propagation,",
    "biological process): its clique genes for core orthogroups, its genes passing the co-expression expression filter otherwise.",
    "Study set: the cluster's orthogroups; background: the orthogroups eligible for the core. One-sided Fisher's exact test,",
    "BH within each cluster and species; terms with padj < 0.05 are listed. padj_size_matched: the same terms against",
    "random sets of eligible orthogroups matched on size (exact; BH); call: 'both' when also significant there, otherwise 'Fisher only'.",
    "majority: significant (padj < 0.05) when an orthogroup carries a term only if at least half of its genes in that species carry it.",
    "Annotated and Significant count orthogroups."))
cat("family-level Table S3:", nrow(sig), "rows; size-matched disagreements (padj 0.05):", R3[agree == FALSE, .N], "of", nrow(R3), "tested terms\n")
