#!/usr/bin/env Rscript
# manuscript_values.R — the value list for the manuscript's sample-stability and validation statements (ledger R042 item 2),
# each read at run time from a producer output; nothing is typed in. One row per value: name, value, source (file and filter).
# Inputs (working tree, written by the driver's producers): results/core_stability_label.tsv, results/loo_stability/loo_gate.tsv,
#   results/stability/{kept_lost_orthogroups,kept_lost_summary,kept_lost_logistic,why_spruce_signal,timepoint_attribution,
#   matched_n_summary}.tsv, results/stability/design/design_model.tsv, results/stability/{stable,sensitive}_string_{obs.txt,null.csv},
#   results/external_validation/{external_validation,external_validation_control}.tsv
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
out <- rbindlist(V); fwrite(out, "results/manuscript_values.tsv", sep = "\t"); cat(nrow(out), "values written\n")
