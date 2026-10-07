#!/usr/bin/env Rscript
# fig4_captions_from_terms.R — the Fig. 4 caption table (cluster, GO.ID) from the approved caption terms
# (scripts/figures/Fig4_terms_approved.tsv: NRS's choice, DECISIONS.md D008; arabidopsis_only marks the terms the legend
# names as Arabidopsis-only). Each term is matched to Table S3 by exact name within its manuscript cluster; the script
# stops if a term is not in Table S3 for that cluster, and reports the genotypes in which each term is listed, so a
# change from the approved state can be seen before the figure is drawn.
# Inputs:  scripts/figures/Fig4_terms_approved.tsv, tables/TableS3_cliqueClusterGO_classic.tsv (scripts/GO/make_TableS3.R)
# Outputs: Fig4_captions_approved.tsv (cluster, GO.ID, Term, arabidopsis_only, genotypes)
# Usage (working directory): Rscript scripts/figures/fig4_captions_from_terms.R
suppressPackageStartupMessages(library(dplyr))
terms <- read.delim("scripts/figures/Fig4_terms_approved.tsv")
s4 <- read.delim("tables/TableS3_cliqueClusterGO_classic.tsv", check.names = FALSE)
names(s4)[names(s4) == "Cluster (manuscript)"] <- "cluster"
hit <- s4 %>% semi_join(terms, by = c("cluster", "Term")) %>% group_by(cluster, Term) %>%
  summarise(GO.ID = paste(unique(GO.ID), collapse = ","), genotypes = paste(Species, collapse = ", "), .groups = "drop")
out <- terms %>% left_join(hit, by = c("cluster", "Term"))
if (anyNA(out$GO.ID)) stop("approved term not in Table S3 for its cluster: ", paste(out$Term[is.na(out$GO.ID)], collapse = "; "))
if (any(grepl(",", out$GO.ID))) stop("a term name maps to more than one GO ID")
out <- out %>% select(cluster, GO.ID, Term, arabidopsis_only, genotypes)
write.table(out, "Fig4_captions_approved.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
