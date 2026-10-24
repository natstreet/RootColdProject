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
#   over the terms with at least one study family, per cluster x genotype; padj < 0.05.
# Step 3, size-matched null: per cluster x genotype, NULL_DRAWS draws of as many annotated background families as the
#   cluster has, matched on size (deciles of the size rule's expressed_genes, sampled within decile without replacement);
#   empirical P = (1 + draws with at least the observed number of study families carrying the term) / (1 + draws), BH over
#   the same terms. Reported for the terms tested in step 2, with whether the two tests agree at padj < 0.05.
#
# Inputs (working tree): data/ComPlEx/cliques/{conserved_genes_COMPLETE_filtered,conserved_genes_COMPLETE_wider}.RDS,
#   clique_HMcluster_list.rda, manuscript_cluster_map.tsv; data/ComPlEx/Orthogroups_20240613.tsv;
#   data/expression/<genotype>_root_cold_stress_expression.txt; data/annotation/*go_ids*; results/orthogroup_size_rule.tsv;
#   <gene-level Table S3 tsv> (first argument).
# Outputs (second argument, a directory): step1_family_count.tsv, step1_cell_summary.tsv, TableS3_family_level.tsv,
#   Table_S3_family_level.xlsx (legend sheet), step3_size_matched.tsv
# Usage: Rscript scripts/GO/tableS3_family_level.R <gene-level TableS3 tsv> <out dir> [null draws, default 5000]
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
    if (x == 0) return(NULL)
    data.table(GO.ID = t, Annotated = k, Significant = x, Expected = round(n * k / N, 2),
               p = phyper(x - 1, k, N - k, n, lower.tail = FALSE)) }))
  r[, padj := p.adjust(p, "BH")][order(padj, p, GO.ID)]
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
res2 <- list(); res3 <- list()
for (g in GT) {
  fg <- fam_genes(g); a2 <- ANN[[g]]
  f2go <- lapply(split(fg$gene, fg$OrthoGroup), function(gs) unique(unlist(a2[intersect(gs, names(a2))], use.names = FALSE)))
  f2go <- f2go[lengths(f2go) > 0]
  tm <- term_members(f2go); U <- tm$units
  terms <- names(tm$members); M <- sparseMatrix(i = match(unlist(tm$members), U), j = rep(seq_along(terms), lengths(tm$members)),
                                                x = 1, dims = c(length(U), length(terms)), dimnames = list(U, terms))
  for (k in sort(unique(core$cluster))) {
    study <- core[cluster == k, OrthoGroup]; st <- intersect(study, U)
    r <- fisher_terms(tm$members, U, st)
    r[, `:=`(cluster = k, genotype = g, families_study = length(st), families_background = length(U),
             family_ids = vapply(GO.ID, function(t) paste(sort(intersect(tm$members[[t]], st)), collapse = ","), ""))]
    res2[[paste(g, k)]] <- r
    ## size-matched null
    obs <- colSums(M[st, r$GO.ID, drop = FALSE]); b <- bins[U]; need <- table(bins[st])
    pool <- split(U, b); cnt <- matrix(0, NULL_DRAWS, nrow(r)); Mr <- M[, r$GO.ID, drop = FALSE]
    for (d in seq_len(NULL_DRAWS)) {
      pick <- unlist(lapply(names(need), function(z) sample(pool[[z]], need[[z]])))
      cnt[d, ] <- colSums(Mr[pick, , drop = FALSE])
    }
    emp <- (1 + colSums(sweep(cnt, 2, obs, ">="))) / (1 + NULL_DRAWS)
    res3[[paste(g, k)]] <- data.table(cluster = k, genotype = g, GO.ID = r$GO.ID, observed = obs, null_mean = round(colMeans(cnt), 2),
                                      empirical_p = emp, padj_size_matched = p.adjust(emp, "BH"), padj_fisher = r$padj)
    cat(sprintf("  %s cluster %d: %d study families (of %d annotated); terms tested %d; Fisher padj<0.05 %d; size-matched padj<0.05 %d\n",
                g, k, length(st), length(U), nrow(r), sum(r$padj < 0.05), sum(res3[[paste(g, k)]]$padj_size_matched < 0.05)))
  }
}
R2 <- rbindlist(res2); R2[, Term := unname(AnnotationDbi::Term(GO.db::GOTERM[GO.ID]))]
R3 <- rbindlist(res3); R3[, Term := R2$Term[match(GO.ID, R2$GO.ID)]]
R3[, agree := (padj_fisher < 0.05) == (padj_size_matched < 0.05)]
sig <- R2[padj < 0.05, .(`Cluster (manuscript)` = cluster, Species = LAB[genotype], GO.ID, Term, Annotated, Significant, Expected,
                         classicFisher = p, padj, families_study, families_background, family_ids)]
sig <- merge(sig, R3[, .(`Cluster (manuscript)` = cluster, Species = LAB[genotype], GO.ID, padj_size_matched)],
             by = c("Cluster (manuscript)", "Species", "GO.ID"), all.x = TRUE)[order(`Cluster (manuscript)`, Species, padj)]
fwrite(sig, file.path(OUT, "TableS3_family_level.tsv"), sep = "\t")
fwrite(R3, file.path(OUT, "step3_size_matched.tsv"), sep = "\t")
write_supp_xlsx(as.data.frame(sig), file.path(OUT, "Table_S3_family_level.xlsx"), "Table S3", sig_cols = c("classicFisher", "padj", "padj_size_matched"),
  legend = paste("Table S3 GO enrichment of the conserved-core clusters with the orthogroup as the unit, per species (provisional",
    "cluster numbering). An orthogroup carries a term in a species when any of its genes there carries it (topGO propagation,",
    "biological process): its clique genes for core orthogroups, its genes passing the co-expression expression filter otherwise.",
    "Study set: the cluster's orthogroups; background: the orthogroups eligible for the core. One-sided Fisher's exact test,",
    "BH within each cluster and species; terms with padj < 0.05 are listed. padj_size_matched: the same terms against",
    NULL_DRAWS, "size-matched random sets of eligible orthogroups (BH). Annotated and Significant count orthogroups."))
cat("family-level Table S3:", nrow(sig), "rows; size-matched disagreements (padj 0.05):", R3[agree == FALSE, .N], "of", nrow(R3), "tested terms\n")
