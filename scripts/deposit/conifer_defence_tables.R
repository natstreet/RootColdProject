#!/usr/bin/env Rscript
# conifer_defence_tables.R — deposit tables for the conifer-only defence result (ledger R064, R066, R068; D089, D091, D094).
#   summary   per conifer: the gene-level GO:0009626 (plant-type hypersensitive response) result among the conifer-only co-expressologs
#             (annotated study genes, orthogroups, padj), the family-level result (study / background families, padj, size-matched padj),
#             and the study genes in NLR-class families with the R066 domain sets and with the extended Pfam NLR models.
#   families  per family: spruce and pine genes, domain counts, class (R066 sets and extended), and the study genes of each conifer.
# Inputs (working tree): results/conifer_only/ (conifer_only_go.R), results/hr_family_domains/ (hr_family_domains.R: family_domains.tsv,
#   hr_genes_in_nlr_families.tsv). The study genes of the two must agree; the script stops otherwise.
# Outputs: tables/conifer_only_defence_summary.tsv, tables/conifer_only_defence_families.tsv
# Usage (working tree): Rscript scripts/deposit/conifer_defence_tables.R
suppressPackageStartupMessages(library(data.table))
TERM <- "GO:0009626"; C1 <- "results/conifer_only"; H <- "results/hr_family_domains"
hr <- fread(file.path(H, "hr_genes_in_nlr_families.tsv")); fd <- fread(file.path(H, "family_domains.tsv"))
S <- rbindlist(lapply(c("spruce", "pine"), function(cf) {
  sg <- fread(file.path(C1, sprintf("gene_level_study_genes_%s.tsv", cf)))[GO.ID == TERM]
  h <- hr[conifer == cf]
  if (!setequal(sg$gene, h$gene)) stop(cf, ": hr_family_domains study genes differ from conifer_only_go.R; rerun hr_family_domains.R")
  g <- fread(file.path(C1, sprintf("gene_level_family_count_%s.tsv", cf)))[GO.ID == TERM]
  f <- fread(file.path(C1, sprintf("family_level_%s.tsv", cf)))[GO.ID == TERM]
  stopifnot(nrow(g) == 1, nrow(f) == 1)
  data.table(conifer = cf, term = TERM, gene_level_study_genes = g$study_genes_annotated, gene_level_orthogroups = g$orthogroups,
             gene_level_padj = g$padj, family_level_study_families = f$Significant, family_level_background_families = f$Annotated,
             family_level_padj = f$padj, family_level_padj_size_matched = f$padj_size_matched,
             study_genes_in_nlr_families = sum(h$nlr_family), study_genes_in_nlr_families_ext = sum(h$nlr_family_ext), study_genes = nrow(h)) }))
fwrite(S, "tables/conifer_only_defence_summary.tsv", sep = "\t")
F <- merge(fd[, !"description"], dcast(hr, OrthoGroup ~ conifer, value.var = "gene", fun.aggregate = length), by = "OrthoGroup", all.x = TRUE)
setnames(F, c("spruce", "pine"), c("spruce_study_genes", "pine_study_genes"), skip_absent = TRUE)
fwrite(F, "tables/conifer_only_defence_families.tsv", sep = "\t"); print(S)
