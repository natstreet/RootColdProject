#!/usr/bin/env Rscript
# paralogue_agreement.R — per core orthogroup and genotype (D083): do that genotype's clique genes (every gene that
# occurs in any of the orthogroup's complete cliques, as for the family-median profile, D082) agree in their cold response?
#   n_genes       clique genes of the genotype in the orthogroup
#   n_de_up       genes that are a DEG (deposit is_deg) with log2FoldChange > 0 at one or more time points
#   n_de_down     genes that are a DEG with log2FoldChange < 0 at one or more time points (a gene can count in both)
#   agreement     "none DE" (no DEG among them); "agree" (one gene, or every DEG in one direction only);
#                 "differ" (DEGs in both directions, counting a gene that is up at one time point and down at another)
#   median_r      median pairwise Pearson r of the genes' profiles (z-scored median expression over the five conditions,
#                 clique_matrix.R); NA with fewer than two profiled genes
#   cluster       the orthogroup's cluster (clique_clusters.R, analysis numbering)
# Inputs: data/ComPlEx/cliques/{conserved_genes_COMPLETE_filtered,conserved_genes_COMPLETE_wider}.RDS,
#         data/ComPlEx/cliques/clique_HMcluster_list.rda, data/expression/med_expression/, data/open/differential_expression.tsv
# Outputs: <out dir>/paralogue_agreement.tsv (one row per orthogroup x genotype), <out dir>/paralogue_agreement_summary.tsv
#          (counts of each agreement value per genotype and per cluster)
# Usage (working tree): Rscript scripts/analysis/paralogue_agreement.R [out dir, default results]
suppressPackageStartupMessages(library(data.table))
out <- if (length(a <- commandArgs(trailingOnly = TRUE))) a[1] else "results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("scripts/ComPlEx/clique_matrix.R")
core <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
G <- as.data.table(clique_genes_by_genotype())
e <- new.env(); load("data/ComPlEx/cliques/clique_HMcluster_list.rda", envir = e)
cl <- as.data.table(get(ls(e)[1], e))[, .(cliqueID = ID, cluster = Cluster)]
cl <- merge(cl, core[, .(cliqueID, OrthoGroup)], by = "cliqueID")
de <- fread("data/open/differential_expression.tsv")[is_deg == TRUE]
if (!all(clique_sp %in% unique(fread("data/open/differential_expression.tsv", select = "species")$species)))
  stop("differential_expression.tsv species codes do not cover ", paste(clique_sp, collapse = " "))
up <- unique(de[log2FoldChange > 0, .(genotype = species, gene)])[, up := TRUE]
dn <- unique(de[log2FoldChange < 0, .(genotype = species, gene)])[, down := TRUE]
G <- merge(merge(G, up, by = c("genotype", "gene"), all.x = TRUE), dn, by = c("genotype", "gene"), all.x = TRUE)
G[is.na(up), up := FALSE][is.na(down), down := FALSE]
prof <- setNames(lapply(clique_sp, function(s) gene_profiles(s, G[genotype == s, gene])), clique_sp)
mr <- function(s, genes) { g <- intersect(genes, rownames(prof[[s]])); if (length(g) < 2) return(NA_real_)
  r <- suppressWarnings(cor(t(prof[[s]][g, , drop = FALSE]))); median(r[upper.tri(r)], na.rm = TRUE) }
res <- G[, .(n_genes = .N, n_de_up = sum(up), n_de_down = sum(down),
             median_r = mr(genotype[1], gene)), by = .(OrthoGroup, genotype)]
res[, agreement := fifelse(n_de_up + n_de_down == 0, "none DE",
                    fifelse(n_genes == 1 | n_de_up == 0 | n_de_down == 0, "agree", "differ"))]
res <- merge(res, cl[, .(OrthoGroup, cluster)], by = "OrthoGroup", all.x = TRUE)
if (res[, uniqueN(OrthoGroup)] != nrow(core) || anyNA(res$cluster)) stop("orthogroups missing from the agreement table or clusters")
setcolorder(res, c("OrthoGroup", "cluster", "genotype")); setorder(res, cluster, OrthoGroup, genotype)
fwrite(res, file.path(out, "paralogue_agreement.tsv"), sep = "\t")
s1 <- res[, .N, by = .(level = genotype, agreement)][, by := "genotype"]
s2 <- res[, .N, by = .(level = as.character(cluster), agreement)][, by := "cluster"]
fwrite(rbind(s1, s2)[order(by, level, agreement), .(by, level, agreement, N)], file.path(out, "paralogue_agreement_summary.tsv"), sep = "\t")
cat(sprintf("%d orthogroups x %d genotypes = %d cells: %s\n", uniqueN(res$OrthoGroup), length(clique_sp), nrow(res),
            paste(names(table(res$agreement)), table(res$agreement), sep = " ", collapse = "; ")))
