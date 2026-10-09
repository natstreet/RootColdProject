#!/usr/bin/env Rscript
# manuscript_values.R — the value list for the manuscript's sample-stability and validation statements (ledger R042 item 2),
# each read at run time from a producer output; nothing is typed in. One row per value: name, value, source (file and filter).
# Inputs (working tree, written by the driver's producers): results/core_stability_label.tsv, results/loo_stability/loo_gate.tsv,
#   results/stability/{kept_lost_orthogroups,kept_lost_summary,kept_lost_logistic,why_spruce_signal,timepoint_attribution,
#   matched_n_summary}.tsv, results/stability/design/design_model.tsv, results/stability/{stable,sensitive}_string_{obs.txt,null.csv},
#   results/external_validation/{external_validation,external_validation_control}.tsv, results/Figure5_values.tsv and
#   results/Figure5c_values.tsv (make_Fig5.R), results/orthogroup_sizes_all.tsv and results/orthogroup_size_rule.tsv (near-core),
#   data/TranscriptionFactors/TF_col0.RData, the col0 comparison tables, results/cbf_check/ (cold_motif_check.R)
# Output: results/manuscript_values.tsv. Usage (working tree): Rscript scripts/analysis/manuscript_values.R
suppressPackageStartupMessages(library(data.table))
V <- list(); add <- function(name, value, source) V[[length(V) + 1]] <<- data.table(name = name, value = as.character(value), source = source)
S <- "results/stability"
lab <- fread("results/core_stability_label.tsv")
add("core_orthogroups", nrow(lab), "results/core_stability_label.tsv: rows")
for (st in c("stable", "sample-sensitive")) add(paste0("core_", st), lab[core_status == st, .N], paste0("core_stability_label.tsv: core_status == ", st))
gate <- fread("results/loo_stability/loo_gate.tsv")
add("gate_comparisons_reproduced", gate[pairs_unmatched == 0 & max_rel_diff < 1e-9, .N], "loo_gate.tsv: pairs_unmatched == 0 & max_rel_diff < 1e-9")
kl <- fread(file.path(S, "kept_lost_orthogroups.tsv"))
for (cl in sort(unique(kl$cluster))) { add(sprintf("cluster%d_orthogroups", cl), kl[cluster == cl, .N], "kept_lost_orthogroups.tsv: by cluster")
  add(sprintf("cluster%d_stable", cl), kl[cluster == cl & core_status == "core", .N], "kept_lost_orthogroups.tsv: by cluster, core_status == core") }
add("de3_stable", kl[core_status == "core" & de_members >= 3, .N], "kept_lost_orthogroups.tsv: de_members >= 3, stable")
add("de3_sensitive", kl[core_status == "sample-sensitive" & de_members >= 3, .N], "kept_lost_orthogroups.tsv: de_members >= 3, sample-sensitive")
sm <- fread(file.path(S, "kept_lost_summary.tsv"))
add("de3_fisher_p", signif(sm[quantity == "DE_in_at_least_3_genotypes", p], 3), "kept_lost_summary.tsv: DE_in_at_least_3_genotypes")
for (q in c("complete_cliques", "mean_expr", "median_base_p", "de_members")) {
  add(paste0(q, "_median_stable"), signif(sm[quantity == q, kept_median], 3), paste0("kept_lost_summary.tsv: ", q))
  add(paste0(q, "_median_sensitive"), signif(sm[quantity == q, lost_median], 3), paste0("kept_lost_summary.tsv: ", q))
  add(paste0(q, "_wilcoxon_p"), signif(sm[quantity == q, p], 3), paste0("kept_lost_summary.tsv: ", q)) }
add("orthologue_carried_detail", sm[quantity == "orthology_class", detail], "kept_lost_summary.tsv: orthology_class")
lg <- fread(file.path(S, "kept_lost_logistic.tsv"))
for (i in 2:nrow(lg)) { add(paste0("logistic_", lg$term[i], "_estimate"), signif(lg$estimate[i], 3), "kept_lost_logistic.tsv")
  add(paste0("logistic_", lg$term[i], "_p"), signif(lg$p[i], 3), "kept_lost_logistic.tsv") }
sg <- fread(file.path(S, "why_spruce_signal.tsv"))[gene_set == "all_filtered"]
for (g in sg$genotype) add(paste0("time_adjR2_median_", g), round(sg[genotype == g, median_time_adjR2], 3), "why_spruce_signal.tsv: all_filtered")
ta <- fread(file.path(S, "timepoint_attribution.tsv"))[core_status == "sample-sensitive", .(rate = sum(losses) / sum(pair_tests)), by = left_out_genotype]
for (g in ta$left_out_genotype) add(paste0("loss_rate_sensitive_leftout_", g), round(ta[left_out_genotype == g, rate], 3), "timepoint_attribution.tsv: sample-sensitive, by left-out genotype")
mn <- fread(file.path(S, "matched_n_summary.tsv"))[set == "random"]
for (i in seq_len(nrow(mn))) add(sprintf("neighbourhood_worst_jaccard_%s_%s", mn$genotype[i], mn$design[i]), mn$median_worst_jaccard[i], "matched_n_summary.tsv: random genes")
dm <- fread(file.path(S, "design/design_model.tsv"))
for (i in seq_len(nrow(dm))) add(sprintf("design_%s_%s_%s", dm$measure[i], dm$type[i], gsub("[()]", "", dm$term[i])), round(dm$estimate[i], 3), "design_model.tsv")
for (st in c("stable", "sensitive")) { o <- as.integer(readLines(file.path(S, paste0(st, "_string_obs.txt")))); nl <- fread(file.path(S, paste0(st, "_string_null.csv")))$null
  add(paste0("string_", st, "_within_pairs"), o, paste0(st, "_string_obs.txt")); add(paste0("string_", st, "_null_mean"), round(mean(nl), 1), paste0(st, "_string_null.csv"))
  add(paste0("string_", st, "_fold"), round(o / mean(nl), 2), "obs / null mean"); add(paste0("string_", st, "_p"), signif((sum(nl >= o) + 1) / (length(nl) + 1), 3), "empirical, (k + 1) / (N + 1)") }
ctl <- fread("results/external_validation/external_validation_control.tsv")
add("external_control_fold", round(ctl$fold, 2), "external_validation_control.tsv"); add("external_control_p", signif(ctl$p, 3), "external_validation_control.tsv")
ev <- fread("results/external_validation/external_validation.tsv")[members_definition == "representative"]
for (st in c("stable", "sensitive", "all")) for (cl in c("cluster_2", "all_clusters")) for (se in c("calixto_up", "calixto_all_de", "zhao_cbf", "jia_cbf")) {
  r <- ev[stratum == st & cluster == cl & set == se]
  add(sprintf("ext_%s_%s_%s", st, cl, se), sprintf("%d/%d; expected %.2f; fold %.2f; q %.2g; %s", r$overlap, r$members, r$expected, r$fold, r$q, r$call),
      "external_validation.tsv: representative members") }
add("external_universe", unique(ev$universe), "external_validation.tsv: universe")
# Fig. 5 (R043): panel a/b values per cluster, genotype and time point; panel c means per genotype
f5 <- fread("results/Figure5_values.tsv")
for (i in seq_len(nrow(f5))) add(sprintf("fig5_cluster%d_%s_%s", f5$cluster[i], f5$genotype[i], f5$timepoint[i]), f5$median_log2fc[i], "Figure5_values.tsv")
f5c <- fread("results/Figure5c_values.tsv")
for (i in seq_len(nrow(f5c))) { add(sprintf("fig5c_%s_cluster_a_mean_24h_10d", f5c$genotype[i]), round(f5c$mean_24h_to_10d_cluster_a[i], 3), "Figure5c_values.tsv")
  add(sprintf("fig5c_%s_cluster_b_mean_24h_10d", f5c$genotype[i]), round(f5c$mean_24h_to_10d_cluster_b[i], 3), "Figure5c_values.tsv") }
# near-core (D050) and its orthogroups with an Arabidopsis TF among the Col-0 genes ComPlEx tested (R043)
sz <- fread("results/orthogroup_sizes_all.tsv"); sr <- fread("results/orthogroup_size_rule.tsv")
nc <- setdiff(sz[eligible == TRUE & all15 == TRUE, OrthoGroup], sr[has_complete_clique == TRUE, OrthoGroup])
add("near_core_orthogroups", length(nc), "orthogroup_sizes_all.tsv all15, without a complete clique (orthogroup_size_rule.tsv)")
ldd <- function(f) { e <- new.env(); load(f, envir = e); as.data.table(get(ls(e)[1], e)) }
tfc <- unique(ldd("data/TranscriptionFactors/TF_col0.RData")$Gene_ID)
tst <- unique(rbindlist(lapply(list.files("data/ComPlEx/PairwiseComparisons", "^comparison-col0-.*-table[.]RData$", full.names = TRUE),
  function(f) ldd(f)[OrthoGroup %in% nc, .(OrthoGroup, gene = Species1)])))
add("near_core_with_col0_tf", tst[gene %in% tfc, uniqueN(OrthoGroup)], "near-core orthogroups with a tested Col-0 gene in TF_col0.RData")
# CBF / cold-panel promoter check (D068)
cc <- fread("results/cbf_check/cold_motif_control.tsv")
add("cbf_control_odds_ratio", sprintf("%.2f (95%% CI %.2f-%.2f); %d/%d vs %d/%d", cc$odds_ratio, cc$ci_low, cc$ci_high, cc$fg_with, cc$fg_genes, cc$bg_with, cc$bg_genes), "cold_motif_control.tsv")
add("cbf_control_zhao_presence_agreement", cc$presence_agreement, "cold_motif_control.tsv")
ce <- fread("results/cbf_check/cold_motif_enrichment.tsv")
for (g in unique(ce$genotype)) { r <- ce[genotype == g & motif == "DRE_CRT"]
  add(paste0("cbf_core_dre_", g), sprintf("%d/%d; OR %.2f (95%% CI %.2f-%.2f); q %.2g; %s", r$fg_with, r$fg_genes, r$odds_ratio, r$ci_low, r$ci_high, r$q, r$call), "cold_motif_enrichment.tsv") }
add("cbf_panel_min_q", signif(min(ce$q), 3), "cold_motif_enrichment.tsv: all 30 tests")
add("cbf_panel_enriched_tests", ce[call == "enriched", .N], "cold_motif_enrichment.tsv")
cb <- fread("results/cbf_check/cbf_genes_status.tsv")
add("cbf1_3_orthogroup_status", paste(unique(paste(cb$OrthoGroup, cb$core_status)), collapse = "; "), "cbf_genes_status.tsv")
# A1 (N159), A2 (N159, R045) and the aspen 3 d DEG count (R045)
a1 <- fread("results/a1_rewiring_duplication.tsv")
add("a1_paralogue_opposed", sprintf("%d/%d", a1$paralogue_opposed, a1$paralogue_total), "a1_rewiring_duplication.tsv")
add("a1_orthologue_opposed", sprintf("%d/%d", a1$orthologue_opposed, a1$orthologue_total), "a1_rewiring_duplication.tsv")
add("a1_fisher", sprintf("OR %.2f (95%% CI %.2f-%.2f), P %.2g", a1$fisher_or, a1$fisher_ci_low, a1$fisher_ci_high, a1$fisher_p), "a1_rewiring_duplication.tsv")
add("a1_adjusted", sprintf("coefficient %.2f, P %.2g", a1$adjusted_coef, a1$adjusted_p), "a1_rewiring_duplication.tsv")
add("a1_cluster1_opposed", sprintf("%d/%d", a1$cluster1_opposed, a1$cluster1_with_direction), "a1_rewiring_duplication.tsv"); add("a1_verdict", a1$verdict, "a1_rewiring_duplication.tsv")
a2 <- fread("results/a2_series_tests.tsv")
for (i in seq_len(nrow(a2))) add(sprintf("a2_%s_%s_%s", a2$sp[i], a2$tissue[i], a2$trt[i]),
  sprintf("cluster 1 median %.2f, cluster 2 median %.2f; T1 q %.2g, T2 q %.2g; pattern %s", a2$c1_median[i], a2$c2_median[i], a2$T1_q[i], a2$T2_q[i], a2$pattern[i]), "a2_series_tests.tsv")
add("a2_series_with_pattern", sprintf("%d of %d", sum(a2$pattern), nrow(a2)), "a2_series_tests.tsv")
ad <- fread("results/aspen_timepoint_degs.tsv")
for (i in seq_len(nrow(ad))) add(paste0("aspen_degs_", ad$timepoint[i]), ad$study_degs[i], "aspen_timepoint_degs.tsv: study DEG call (equal to padj < 0.05 alone at every time point)")
stopifnot(all(ad$study_degs == ad$padj_below_0.05))
out <- rbindlist(V); fwrite(out, "results/manuscript_values.tsv", sep = "\t"); cat(nrow(out), "values written\n")
