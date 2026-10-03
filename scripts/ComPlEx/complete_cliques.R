#!/usr/bin/env Rscript
# complete_cliques.R — the complete co-expressolog cliques: for each orthogroup, every maximal clique of at least six
# genes in the network of its co-expressologs (Max.p.val < 0.1), and the complete cliques (one gene in each of the six
# genotypes, all 15 pairs co-expressologs) in wide form.
#
# This is the method of completeCliques.R in the analysis repository the published objects came from
# (ColdStressRoots, R/ComPlEx_and_Processing/completeCliques.R, 2025-11-17), with paths made relative and its results
# collected by index rather than appended (same order). Kept as it was, including two details that shape the
# published objects: orthogroups need at least three co-expressolog pairs, the eight largest orthogroups are left out
# (discard list), and MeanCliqueSum is the mean Max.p.val over the clique's six genes, each gene taking the Max.p.val
# of its first listed pair in the clique (dplyr::distinct keeps the first row per gene). Cliques with exactly equal
# MeanCliqueSum are ordered by clique number (the original left their order to row order, so the conserved-core
# step, which takes the first, could pick either; this reproduces the published pick in 115 of 116 orthogroups).
#
# Input:  data/ComPlEx/all_expressed_genes.RData (scripts/ComPlEx/coexpressolog_tables.R)
# Output: <out>/clique_genes_filterable_COMPLETE.RDS  every pair of every clique (cliqueID, species_pair, OriginalCliqueSize)
#         <out>/conserved_genes_COMPLETE_wider.RDS    one row per complete clique (gene per genotype, MeanCliqueSum)
#         <out>/conserved_complete_groups.RDS         the orthogroups with a complete clique
# Usage (repository root): Rscript scripts/ComPlEx/complete_cliques.R [out dir, default data/ComPlEx/cliques]
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(stringr); library(igraph) })
a <- commandArgs(trailingOnly = TRUE)
out <- if (length(a)) a[1] else "data/ComPlEx/cliques"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
load("data/ComPlEx/all_expressed_genes.RData")                                     # expr_genes

discard_list <- c("OG0000000", "OG0000001", "OG0000003", "OG0000004", "OG0000006", "OG0000008", "OG0000013", "OG0000014")
pval <- 0.1
co_expressologs <- expr_genes %>%
  dplyr::rename(OrthoGroup = OrthologGroup) %>%
  dplyr::filter(!(OrthoGroup %in% discard_list)) %>%
  dplyr::filter(MaxpVal < pval) %>%
  dplyr::group_by(OrthoGroup) %>%
  dplyr::mutate(OG_size = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(OG_size)) %>%
  dplyr::filter(OG_size >= 3)
OGs_to_use <- unique(co_expressologs$OrthoGroup)

pieces <- list(); k <- 0
largest_clique_ok <- NULL
for (i in seq_along(OGs_to_use)) {
  g <- OGs_to_use[i]
  co_g <- co_expressologs %>% dplyr::filter(OrthoGroup == g)
  co_g <- co_g %>% dplyr::mutate(UGeneSpecies1 = paste0(co_g$Species1, "-", co_g$GeneSpecies1),
                                 UGeneSpecies2 = paste0(co_g$Species2, "-", co_g$GeneSpecies2))
  nodes <- data.frame(name = unique(c(co_g$UGeneSpecies1, co_g$UGeneSpecies2)))
  edges <- data.frame(from = co_g$UGeneSpecies1, to = co_g$UGeneSpecies2)
  net <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  largest_clique <- max_cliques(net, min = 6)
  if (length(largest_clique) > 0) largest_clique_ok <- largest_clique   # as the original: an orthogroup without
  if (is.null(largest_clique_ok)) next                                   # cliques reuses the previous list, whose
  for (c in seq_along(largest_clique_ok)) {                              # genes then match no pair here
    clique_name <- attr(largest_clique_ok[[c]], "names")
    k <- k + 1
    pieces[[k]] <- co_g %>%
      dplyr::filter(UGeneSpecies1 %in% clique_name & UGeneSpecies2 %in% clique_name) %>%
      dplyr::select(OrthoGroup, Species1, GeneSpecies1, Species2, GeneSpecies2, MaxpVal) %>%
      dplyr::mutate(cliqueID = paste0(g, "_", c))
  }
  if (i %% 500 == 0) cat(i, "of", length(OGs_to_use), "orthogroups\n")
}
clique_genes_unlisted <- as.data.frame(dplyr::bind_rows(pieces))
complete_cliques_filterable <- clique_genes_unlisted %>%
  dplyr::mutate(species_pair = paste0(Species1, Species2)) %>%
  dplyr::group_by(cliqueID) %>%
  dplyr::mutate(OriginalCliqueSize = n()) %>%
  dplyr::ungroup()
saveRDS(complete_cliques_filterable, file.path(out, "clique_genes_filterable_COMPLETE.RDS"))

caas <- complete_cliques_filterable %>%
  dplyr::filter(OriginalCliqueSize == 15) %>%
  dplyr::select(Species1, GeneSpecies1, Species2, GeneSpecies2, OrthoGroup, cliqueID, MaxpVal)
spc1 <- caas %>% dplyr::select(Species1, GeneSpecies1, OrthoGroup, cliqueID, MaxpVal) %>%
  dplyr::mutate(UGeneSpecies1 = paste0(caas$Species1, "-", caas$GeneSpecies1)) %>%
  dplyr::rename(Genes = GeneSpecies1, UGenes = UGeneSpecies1, Species = Species1)
spc2 <- caas %>% dplyr::select(Species2, GeneSpecies2, OrthoGroup, cliqueID, MaxpVal) %>%
  dplyr::mutate(UGeneSpecies2 = paste0(caas$Species2, "-", caas$GeneSpecies2)) %>%
  dplyr::rename(Genes = GeneSpecies2, UGenes = UGeneSpecies2, Species = Species2)
caas_wider <- rbind(spc1, spc2) %>%
  dplyr::group_by(cliqueID) %>%
  dplyr::distinct(UGenes, .keep_all = TRUE) %>%
  dplyr::mutate(MeanCliqueSum = mean(MaxpVal)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(UGenes, MaxpVal)) %>%
  tidyr::pivot_wider(names_from = Species, values_from = Genes) %>%
  dplyr::select(OrthoGroup, cliqueID, MeanCliqueSum, col0, ost0, aspen, birch, spruce, pine) %>%
  dplyr::group_by(OrthoGroup) %>%
  dplyr::arrange(OrthoGroup, MeanCliqueSum, as.integer(sub("^.*_", "", cliqueID))) %>%   # exact ties: lower clique number
  dplyr::ungroup() %>%
  dplyr::mutate(conservedType = "complete")
saveRDS(caas_wider, file.path(out, "conserved_genes_COMPLETE_wider.RDS"))
saveRDS(caas_wider %>% dplyr::distinct(OrthoGroup) %>% dplyr::pull(OrthoGroup), file.path(out, "conserved_complete_groups.RDS"))
cat(sprintf("%d orthogroups searched; %d cliques (%d pairs); %d complete cliques in %d orthogroups -> %s\n",
            length(OGs_to_use), length(unique(complete_cliques_filterable$cliqueID)), nrow(complete_cliques_filterable),
            nrow(caas_wider), length(unique(caas_wider$OrthoGroup)), out))
