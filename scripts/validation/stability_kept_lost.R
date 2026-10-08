#!/usr/bin/env Rscript
# stability_kept_lost.R — what the sample-stability rule (D052, D059) keeps and loses, per orthogroup, with evidence that
# does not come from the co-expression test itself (ledger R041; NRS: a decision grounded in biology).
#
# Per core orthogroup (133, conserved_genes_COMPLETE_filtered.RDS), using the clique that represents it (the new
# representative for orthogroups that stay core, the old one for sample-sensitive orthogroups):
#   status, manuscript cluster, complete cliques, passing cliques, orthology class (orthology_class.tsv)
#   de_members        members that are cold DEGs in their own genotype (DEGs_<genotype>, union over time points; 0-6)
#   peak_sign_agree   members whose largest-|log2FC| time point has the majority sign (direction of the cold response; 0-6)
#   mean_expr         mean over members of the member's mean expression (the ComPlEx input values)
#   median_base_p     median unmodified Max.p.val of the clique's 15 pairs
#   col0 gene, Arabidopsis symbol and name (from the comparison table annotation)
# Summary: kept vs lost for each quantity (medians, Wilcoxon or Fisher P; descriptive, no correction).
# Inputs (working tree): results/loo_stability/{loo_orthogroups.tsv, loo_pairs_worst.tsv.gz}, data/ComPlEx/cliques/
#   {conserved_genes_COMPLETE_filtered.RDS, conserved_genes_COMPLETE_wider.RDS, clique_HMcluster_list.rda,
#   manuscript_cluster_map.tsv}, results/orthology_class.tsv, data/DEGs/{DEGs_<g>, log2FoldChange_<code>}.RData,
#   data/expression/<g>_root_cold_stress_expression.txt, data/ComPlEx/PairwiseComparisons/comparison-col0-aspen-*.RData
# Outputs (<out dir>): kept_lost_orthogroups.tsv, kept_lost_summary.tsv
# Usage (working tree): Rscript stability_kept_lost.R <out dir>
suppressPackageStartupMessages(library(data.table))
OUT <- commandArgs(trailingOnly = TRUE)[1]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
G <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
CODE <- c(col0 = "col0", ost0 = "ost0", aspen = "pt", birch = "bp", spruce = "pa", pine = "ps")
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
og <- fread("results/loo_stability/loo_orthogroups.tsv")
old <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
w <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS"))
og[, clique := ifelse(core_status == "core", representative, old$cliqueID[match(OrthoGroup, old$OrthoGroup)])]
cl <- as.data.table(ld("data/ComPlEx/cliques/clique_HMcluster_list.rda")); cl[, OrthoGroup := sub("_.*", "", ID)]
mp <- fread("data/ComPlEx/cliques/manuscript_cluster_map.tsv")
og[, cluster := mp$manuscript_cluster[match(cl$Cluster[match(OrthoGroup, cl$OrthoGroup)], mp$analysis_cluster)]]
oc <- fread("results/orthology_class.tsv"); og[, orthology_class := oc$orthology_class[match(OrthoGroup, oc$OrthoGroup)]]
mem <- melt(w[cliqueID %in% og$clique, c("cliqueID", G), with = FALSE], id.vars = "cliqueID", variable.name = "genotype", value.name = "gene")
mem[, genotype := as.character(genotype)]
for (g in G) {
  degs <- ld(sprintf("data/DEGs/DEGs_%s.RData", g)); lf <- ld(sprintf("data/DEGs/log2FoldChange_%s.RData", CODE[[g]]))
  top <- vapply(lf$TopL2FC, function(v) if (length(v) && !is.na(v[1])) unname(v[1]) else NA_real_, numeric(1)); names(top) <- rownames(lf)
  ex <- fread(sprintf("data/expression/%s_root_cold_stress_expression.txt", g)); mexp <- setNames(rowMeans(as.matrix(ex[, -1])), ex[[1]])
  mem[genotype == g, `:=`(de = gene %in% degs, top_lfc = top[gene], mean_expr = mexp[gene])]
}
per <- mem[, .(de_members = sum(de), peak_sign_agree = { s <- sign(top_lfc[!is.na(top_lfc) & top_lfc != 0]); if (length(s)) max(sum(s > 0), sum(s < 0)) else 0L },
               mean_expr = mean(mean_expr, na.rm = TRUE)), by = cliqueID]
og <- merge(og, per, by.x = "clique", by.y = "cliqueID", all.x = TRUE)
pw <- fread(cmd = "gzip -dc results/loo_stability/loo_pairs_worst.tsv.gz")
pairs <- rbindlist(lapply(seq_len(length(G) - 1), function(i) rbindlist(lapply((i + 1):length(G), function(j)
  w[cliqueID %in% og$clique, .(cliqueID, s1 = G[i], s2 = G[j], Species1 = get(G[i]), Species2 = get(G[j]))]))))
pairs <- merge(pairs, pw[, .(s1, s2, Species1, Species2, baseline)], by = c("s1", "s2", "Species1", "Species2"), all.x = TRUE)
og <- merge(og, pairs[, .(median_base_p = median(baseline, na.rm = TRUE)), by = cliqueID], by.x = "clique", by.y = "cliqueID", all.x = TRUE)
ct <- as.data.table(ld("data/ComPlEx/PairwiseComparisons/comparison-col0-aspen-abspearsonMR0.03-table.RData"))
ann <- unique(ct[, .(gene = Species1, symbol = Symbol, name = Name)])[, .SD[1], by = gene]
og[, col0_gene := w$col0[match(clique, w$cliqueID)]]
og <- merge(og, ann, by.x = "col0_gene", by.y = "gene", all.x = TRUE)
setcolorder(og, c("OrthoGroup", "core_status", "cluster", "complete_cliques", "passing_cliques", "orthology_class", "de_members",
                  "peak_sign_agree", "mean_expr", "median_base_p", "col0_gene", "symbol", "name"))
fwrite(og[order(core_status, cluster, OrthoGroup)], file.path(OUT, "kept_lost_orthogroups.tsv"), sep = "\t")
k <- og[core_status == "core"]; l <- og[core_status == "sample-sensitive"]
num <- function(v) data.table(quantity = v, kept_median = median(k[[v]], na.rm = TRUE), lost_median = median(l[[v]], na.rm = TRUE),
                              p = wilcox.test(k[[v]], l[[v]], exact = FALSE)$p.value)
sm <- rbindlist(lapply(c("complete_cliques", "de_members", "peak_sign_agree", "mean_expr", "median_base_p"), num))
ft <- function(lab, x) { t <- table(og$core_status, x); data.table(quantity = lab, kept_median = NA, lost_median = NA, p = if (min(dim(t)) > 1) fisher.test(t)$p.value else NA_real_,
                         detail = paste(apply(t, 1, function(r) paste(names(r), r, sep = "=", collapse = ",")), rownames(t), sep = " @ ", collapse = " | ")) }
sm <- rbind(sm[, detail := NA_character_], ft("orthology_class", og$orthology_class), ft("cluster", og$cluster), ft("DE_in_at_least_3_genotypes", og$de_members >= 3))
fwrite(sm, file.path(OUT, "kept_lost_summary.tsv"), sep = "\t"); print(sm, row.names = FALSE)
