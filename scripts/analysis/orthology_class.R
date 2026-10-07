#!/usr/bin/env Rscript
# orthology_class.R — the orthology class of each orthogroup with a complete clique (DECISIONS.md D043): the N1
# (seed-plant) test classifies, it does not filter. An orthogroup is orthologue-carried when at least one of its
# complete cliques has all six members in one N1 hierarchical orthogroup (HOG), i.e. descended from one gene of the
# seed-plant ancestor; otherwise paralogue-carried. The class is a call from gene trees, with error in both directions
# (a wrong tree can split true orthologues, an unresolved old duplication can join paralogues); it probably overstates
# the paralogue-carried cases. orthology_basis records which trees the N1 HOGs of each orthogroup were inferred from.
#
# Inputs:  <hog prefix>_per_orthogroup.tsv and <hog prefix>_representatives.tsv (scripts/analysis/hog_clique_consistency.R
#          on the N1 table), results/orthogroup_size_rule.tsv (scripts/ComPlEx/size_rule.R),
#          data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS and conserved_genes_COMPLETE_filtered_flagged.RDS (core_cliques.R),
#          <msa tree dir>: the multiple-alignment gene trees used for the N1 HOGs (<orthogroup>_tree.txt)
# Output:  results/orthology_class.tsv, one row per orthogroup with a complete clique: OrthoGroup, counted_as_conserved,
#          expressed_genes, p_core_given_size, n_complete_cliques, n_orthology_consistent_cliques, orthology_class,
#          representative_clique, representative_orthology_consistent, best_consistent_clique and its six members
#          (best_col0 ... best_pine), orthology_basis (msa_tree / dendroblast_tree)
# Usage (repository root): Rscript scripts/analysis/orthology_class.R <hog prefix> <msa tree dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 2) stop("usage: orthology_class.R <hog prefix> <msa tree dir>")
G <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
per  <- fread(paste0(a[1], "_per_orthogroup.tsv"))
reps <- fread(paste0(a[1], "_representatives.tsv"))
rule <- fread("results/orthogroup_size_rule.tsv")[has_complete_clique == TRUE]
core <- as.data.table(rbind(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"),
                            readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered_flagged.RDS")))
if (!setequal(rule$OrthoGroup, per$orthogroup) || !setequal(rule$OrthoGroup, core$OrthoGroup))
  stop("the HOG table, the size rule and the core do not cover the same orthogroups")
chk <- merge(core[, .(orthogroup = OrthoGroup, cliqueID)], reps[, .(orthogroup, representative_clique)], by = "orthogroup")
if (any(chk$cliqueID != chk$representative_clique)) stop("representative cliques differ from core_cliques.R's pick")

d <- merge(rule[, .(OrthoGroup, counted_as_conserved, expressed_genes, p_core_given_size)],
           per[, .(OrthoGroup = orthogroup, n_complete_cliques = complete_cliques,
                   n_orthology_consistent_cliques = hog_consistent_cliques)], by = "OrthoGroup")
d[, orthology_class := ifelse(n_orthology_consistent_cliques > 0, "orthologue-carried", "paralogue-carried")]
r <- reps[, c(list(OrthoGroup = orthogroup, representative_clique = representative_clique,
                   representative_orthology_consistent = representative_consistent,
                   best_consistent_clique = best_consistent_clique), setNames(.SD, paste0("best_", G))), .SDcols = G]
d <- merge(d, r, by = "OrthoGroup")
d[, orthology_basis := ifelse(file.exists(file.path(a[2], paste0(OrthoGroup, "_tree.txt"))), "msa_tree", "dendroblast_tree")]
setorder(d, -counted_as_conserved, -expressed_genes, OrthoGroup)
fwrite(d, "results/orthology_class.tsv", sep = "\t", na = "")

cat("orthogroups with a complete clique:", nrow(d), "| counted as conserved:", sum(d$counted_as_conserved), "\n")
print(d[, .N, by = .(counted_as_conserved, orthology_class, orthology_basis)], row.names = FALSE)
cat("counted, orthologue-carried, representative mixes paralogues:",
    d[counted_as_conserved == TRUE & orthology_class == "orthologue-carried" & representative_orthology_consistent == FALSE, .N], "\n")
