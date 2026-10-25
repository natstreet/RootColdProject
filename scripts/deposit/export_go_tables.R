#!/usr/bin/env Rscript
# export_go_tables.R — the GO enrichment tables of the data deposit (data/open/), written from the supplementary
# tables so that the two agree row for row and value for value.
#
# Run make_TableS1.R, make_TableS2.R and tableS3_arabidopsis.R (via the driver) first; this script reads their workbooks.
#
# Inputs:
#   tables/Table_S1.xlsx   GO enrichment of the commonly differentially expressed orthogroups
#   tables/Table_S2.xlsx   GO enrichment of the conserved co-expressolog clade sets
#   tables/Table_S3.xlsx   GO enrichment of the conserved-core clique clusters on the Col-0 annotation (manuscript
#                          numbering; tableS3_arabidopsis.R, D090): sheets "Table S3 family level" and "Table S3 gene level"
# Outputs (tab-separated, one header row):
#   <outdir>/go_enrichment_common_degs.tsv
#   <outdir>/go_enrichment_conserved_sets.tsv
#   <outdir>/go_enrichment_clusters.tsv              (family level)
#   <outdir>/go_enrichment_clusters_gene_level.tsv
#
# Usage: Rscript scripts/deposit/export_go_tables.R [outdir]      # default outdir: tables/open

suppressMessages({ library(readxl); library(dplyr) })

out <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(out)) out <- "tables/open"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

GENOTYPE <- c(Col0 = "col0", Ost0 = "ost0", Aspen = "aspen", Birch = "birch", Spruce = "spruce", Pine = "pine")
read_table <- function(n, sheet = sprintf("Table S%d", n)) read_excel(sprintf("tables/Table_S%d.xlsx", n), sheet = sheet)
genotype <- function(x) { g <- unname(GENOTYPE[x]); stopifnot(!anyNA(g)); g }
wr <- function(d, f) {
  write.table(d, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat(sprintf("%s: %d rows\n", f, nrow(d)))
}

s1 <- read_table(1)
wr(transmute(s1, species = "col0", go_id = GO.ID, term = Term, annotated = Annotated, significant = Significant,
             expected = Expected, p_classic_fisher = `classicFisher P`, padj),
   "go_enrichment_common_degs.tsv")

s2 <- read_table(2)
wr(transmute(s2, species = genotype(Species), conserved_set = `Co-expressolog set`, go_id = GO.ID, term = Term,
             annotated = Annotated, significant = Significant, expected = Expected,
             p_classic_fisher = `classicFisher P`, padj),
   "go_enrichment_conserved_sets.tsv")

s3 <- read_table(3, "Table S3 family level")
wr(transmute(s3, species = genotype(Species), cluster = `Cluster (manuscript)`, go_id = GO.ID, term = Term,
             annotated_families = Annotated, significant_families = Significant, expected = Expected,
             p_classic_fisher = classicFisher, padj, padj_size_matched, size_matched, padj_majority, majority,
             study_families, background_families, family_ids),
   "go_enrichment_clusters.tsv")

s3g <- read_table(3, "Table S3 gene level")
wr(transmute(s3g, species = genotype(Species), cluster = `Cluster (manuscript)`, go_id = GO.ID, term = Term,
             annotated_genes = Annotated, significant_genes = Significant, expected = Expected,
             p_classic_fisher = classicFisher, padj, orthogroups, study_genes, background_genes),
   "go_enrichment_clusters_gene_level.tsv")
