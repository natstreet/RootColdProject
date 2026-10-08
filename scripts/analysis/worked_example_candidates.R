#!/usr/bin/env Rscript
# worked_example_candidates.R — candidates for the Fig. S1 worked example under D066 (ledger R042 item 3): sample-stable,
# orthologue-carried with an orthology-consistent representative, at most 20 expressed genes in the orthogroup, and the
# representative clique's members cold DEGs in at least 3 genotypes including spruce or pine; ranked by DE genotypes, then
# expressed genes (fewer first). Representative = the 133 core's (lowest MeanCliqueSum).
# Inputs (working tree): results/orthology_class.tsv, <interpretation dir>/kept_lost_orthogroups.tsv, data/DEGs/DEGs_<g>.RData,
#   data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS
# Output: <interpretation dir>/worked_example_candidates.tsv (all qualifying, ranked) and _relaxed.tsv (information only). Usage: Rscript worked_example_candidates.R <interpretation dir>
suppressPackageStartupMessages(library(data.table))
D <- commandArgs(trailingOnly = TRUE)[1]; G <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
oc <- fread("results/orthology_class.tsv"); kl <- fread(file.path(D, "kept_lost_orthogroups.tsv"))
rep <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
for (g in G) { d <- ld(sprintf("data/DEGs/DEGs_%s.RData", g)); rep[, (paste0("de_", g)) := get(g) %in% d] }
rep[, de_genotypes := rowSums(.SD), .SDcols = paste0("de_", G)][, de_conifer := de_spruce | de_pine]
rep[, de_in := apply(.SD, 1, function(r) paste(G[r], collapse = ";")), .SDcols = paste0("de_", G)]
x <- merge(rep[, .(OrthoGroup, representative = cliqueID, col0, de_genotypes, de_conifer, de_in)],
           oc[, .(OrthoGroup, expressed_genes, orthology_class, representative_orthology_consistent)], by = "OrthoGroup")
x <- merge(x, kl[, .(OrthoGroup, core_status, cluster, symbol)], by = "OrthoGroup")
cand <- x[core_status == "core" & orthology_class == "orthologue-carried" & representative_orthology_consistent == TRUE &
          expressed_genes <= 20 & de_genotypes >= 3 & de_conifer == TRUE][order(-de_genotypes, expressed_genes)]
fwrite(cand, file.path(D, "worked_example_candidates.tsv"), sep = "\t")
# relaxed, for information only (not the R042 criteria): DE in at least 2 genotypes (one a conifer) or up to 30 expressed genes
rel <- x[core_status == "core" & orthology_class == "orthologue-carried" & representative_orthology_consistent == TRUE &
         expressed_genes <= 30 & de_genotypes >= 2 & de_conifer == TRUE][order(-de_genotypes, expressed_genes)]
rel[, meets_r042 := OrthoGroup %in% cand$OrthoGroup]
fwrite(rel, file.path(D, "worked_example_candidates_relaxed.tsv"), sep = "\t")
cat("relaxed (DE >= 2 incl. a conifer, <= 30 expressed genes):", nrow(rel), "\n"); print(rel[, .(OrthoGroup, col0, symbol, cluster, expressed_genes, de_genotypes, de_in, meets_r042)])
cat("qualifying:", nrow(cand), "\n"); print(cand[, .(OrthoGroup, col0, symbol, cluster, expressed_genes, de_genotypes, de_in)])
