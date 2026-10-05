#!/usr/bin/env Rscript
# lineage_gene_sets.R — the conserved and lineage-specific co-expressolog gene sets (deposit:
# data/ComPlEx/Common_OG_genes_<sp>_005.RData and Unique_OG_genes_<sp>_005.RData).
#
# Ported from Co-expressologAnalysis.Rmd in the analysis repository the published objects came from (ColdStressRoots,
# R/ComPlEx_and_Processing, E. van Zalen, 2025-11-17), chunks "sets of co-expressologs", "enrich clade genes" and
# "extract genes". The co-expressologs within each lineage (Arabidopsis: Col-0 and Ost-0; broadleaf trees: aspen and
# birch; conifers: spruce and pine) at Max.p.val < 0.05 give three sets of orthogroups.
#   Common (conserved): the genes of the co-expressolog pairs, Max.p.val < 0.05, of orthogroups found in all three
#                       lineages.
#   Unique (lineage-specific): the genes of a lineage's co-expressolog pairs in orthogroups found in that lineage only.
# Each set is a vector of gene identifiers per genotype, as deposited.
#
# Input:  data/ComPlEx/co_expressologs.RData (scripts/ComPlEx/coexpressolog_tables.R)
# Output: <out>/Common_OG_genes_<sp>_005.RData, <out>/Unique_OG_genes_<sp>_005.RData (object names as deposited)
# Usage (repository root): Rscript scripts/ComPlEx/lineage_gene_sets.R [out dir, default data/ComPlEx]
suppressMessages(library(dplyr))
a <- commandArgs(trailingOnly = TRUE); out <- if (length(a)) a[1] else "data/ComPlEx"
load("data/ComPlEx/co_expressologs.RData")                                   # co_expressologs
lineage <- function(s1, s2) co_expressologs %>%
  dplyr::filter(Species1 == s1 & Species2 == s2) %>%
  dplyr::filter(MaxpVal < 0.05) %>%
  dplyr::rename(OrthoGroup = OrthologGroup) %>%
  dplyr::mutate(UGeneSpecies1 = paste0(Species1, "-", GeneSpecies1), UGeneSpecies2 = paste0(Species2, "-", GeneSpecies2))
ara <- lineage("col0", "ost0"); broad <- lineage("aspen", "birch"); gymno <- lineage("spruce", "pine")
set1 <- unique(ara$OrthoGroup); set2 <- unique(broad$OrthoGroup); set3 <- unique(gymno$OrthoGroup)
n_abc <- unique(intersect(set1, unique(intersect(set2, set3))))
n_ab <- setdiff(unique(intersect(set1, set2)), n_abc)
n_bc <- setdiff(unique(intersect(set2, set3)), n_abc)
n_ac <- setdiff(unique(intersect(set1, set3)), n_abc)
n_a <- setdiff(set1, c(n_ab, n_bc, n_ac, n_abc)); n_b <- setdiff(set2, c(n_ab, n_bc, n_ac, n_abc))
n_c <- setdiff(set3, c(n_ab, n_bc, n_ac, n_abc))
genes_of <- function(tab, ogs) tab %>% dplyr::filter(OrthoGroup %in% ogs) %>%
  dplyr::select(UGeneSpecies1, UGeneSpecies2) %>% unlist() %>% unique()
uniq <- list(ara = genes_of(ara, n_a), broad = genes_of(broad, n_b), conif = genes_of(gymno, n_c))
common <- rbind(ara, broad, gymno) %>% dplyr::filter(MaxpVal < 0.05) %>% genes_of(n_abc)
pick <- function(v, sp) gsub(".*-", "", v[grepl(paste0(sp, "-"), v)])
sets <- list(Unique = list(col0 = pick(uniq$ara, "col0"), ost0 = pick(uniq$ara, "ost0"), aspen = pick(uniq$broad, "aspen"),
                           birch = pick(uniq$broad, "birch"), spruce = pick(uniq$conif, "spruce"), pine = pick(uniq$conif, "pine")),
             Common = setNames(lapply(c("col0", "ost0", "aspen", "birch", "spruce", "pine"), function(s) pick(common, s)),
                               c("col0", "ost0", "aspen", "birch", "spruce", "pine")))
objname <- list(Unique = c(col0 = "A_col0", ost0 = "A_ost0", aspen = "B_Asp", birch = "B_Bir", spruce = "C_Spr", pine = "C_Pin"),
                Common = c(col0 = "ABC_col0", ost0 = "ABC_ost0", aspen = "ABC_Asp", birch = "ABC_Bir", spruce = "ABC_Spr", pine = "ABC_Pin"))
for (k in names(sets)) for (s in names(sets[[k]])) {
  e <- new.env(); assign(objname[[k]][[s]], sets[[k]][[s]], envir = e)
  save(list = objname[[k]][[s]], envir = e, file = file.path(out, sprintf("%s_OG_genes_%s_005.RData", k, s)))
}
cat(sprintf("orthogroups: Arabidopsis %d, broadleaf %d, conifer %d; in all three %d\n", length(set1), length(set2), length(set3), length(n_abc)))
cat("genes, conserved (Common):", paste(names(sets$Common), lengths(sets$Common), sep = "=", collapse = " "), "\n")
cat("genes, lineage-specific (Unique):", paste(names(sets$Unique), lengths(sets$Unique), sep = "=", collapse = " "), "\n")
