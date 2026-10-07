#!/usr/bin/env Rscript
# loo_stability_score.R — score the sample-stability test of the conserved core (DECISIONS.md D048, D052, D059; ledger
# R039) from the leave-one-sample-out ComPlEx runs.
#
# Input runs (scripts/validation/loo_task.R, one per task of tasks.tsv): for each of the 15 comparisons, one unmodified
# run (sample "none") and one run per sample of either genotype left out; each keeps Species1, Species2, OrthoGroup and
# Max.p.val (as ComPlEx.R writes it, 3 significant figures) for the orthogroups with a complete clique.
# Rules, fixed before any run (N115, R039, D059):
#   gate     every unmodified run reproduces the deposited comparison table's Max.p.val for every kept pair (relative
#            difference < 1e-9) and lists the same pairs; otherwise the script stops
#   worst    for a pair in comparison (A, B): the largest Max.p.val over the runs leaving out any one sample of A or of
#            B; a pair absent from a run (not tested there) counts as 1
#   member   a clique member passes if its pair with the clique partner is below 0.1 in the worst run in at least 4 of
#            its 5 comparisons
#   clique   passes if all six members pass
#   core     a core orthogroup stays if at least one complete clique passes; its representative is then the passing
#            clique with the lowest MeanCliqueSum (ties by clique number); otherwise it is sample-sensitive
# Inputs (working tree): data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS, _filtered.RDS,
#   data/ComPlEx/PairwiseComparisons/comparison-<s1>-<s2>-abspearsonMR0.03-table.RData; <run dir> with tasks.tsv and
#   out/task_NNN.tsv.gz
# Outputs (<out dir>): loo_gate.tsv, loo_pairs_worst.tsv.gz, loo_members.tsv.gz, loo_cliques.tsv.gz,
#   loo_orthogroups.tsv (core_status, representative, what failed), loo_summary.tsv
# Usage (working tree): Rscript scripts/validation/loo_stability_score.R <run dir> <out dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); RUN <- a[1]; OUT <- a[2]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P_THR <- 0.1; MIN_ROBUST <- 4; REL_TOL <- 1e-9
G <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
tk <- fread(file.path(RUN, "tasks.tsv"))
rd <- function(id) { f <- file.path(RUN, "out", sprintf("task_%03d.tsv.gz", id)); if (!file.exists(f)) stop("missing output of task ", id)
  x <- fread(cmd = paste("gzip -dc", shQuote(f)), colClasses = c(Species1 = "character", Species2 = "character", Max.p.val = "character")); x[, p := as.numeric(Max.p.val)]; x }
ld <- function(f) { e <- new.env(); load(f, envir = e); as.data.table(get(ls(e)[1], e)) }

## gate: unmodified runs against the deposited tables
gate <- rbindlist(lapply(tk[sample == "none", task], function(id) {
  x <- tk[task == id]; b <- rd(id)
  d <- ld(sprintf("data/ComPlEx/PairwiseComparisons/comparison-%s-%s-abspearsonMR0.03-table.RData", x$s1, x$s2))[OrthoGroup %in% unique(b$OrthoGroup)]
  d[, p := as.numeric(Max.p.val)]; m <- merge(b[, .(Species1, Species2, p_run = p)], d[, .(Species1, Species2, p_dep = p)], by = c("Species1", "Species2"), all = TRUE)
  data.table(comparison = paste(x$s1, x$s2, sep = "-"), pairs_run = nrow(b), pairs_deposited = nrow(d), pairs_unmatched = m[is.na(p_run) | is.na(p_dep), .N],
             max_rel_diff = m[!is.na(p_run) & !is.na(p_dep), max(abs(p_run - p_dep) / pmax(p_dep, 1e-300))]) }))
fwrite(gate, file.path(OUT, "loo_gate.tsv"), sep = "\t"); print(gate, row.names = FALSE)
if (nrow(gate) != 15 || any(gate$pairs_unmatched > 0) || any(gate$max_rel_diff >= REL_TOL)) stop("the gate fails: unmodified runs do not reproduce the deposited tables")

## worst leave-one-out Max.p.val per pair and comparison, over both genotypes' samples (D059)
worst <- rbindlist(lapply(split(tk[sample != "none"], by = c("s1", "s2")), function(t) {
  r <- rbindlist(lapply(t$task, function(id) rd(id)[, .(Species1, Species2, p)]))
  w <- r[, .(worst = max(p), runs = .N), by = .(Species1, Species2)]
  w[runs < nrow(t), worst := 1]                                    # absent from a run: not tested there
  w[, `:=`(s1 = t$s1[1], s2 = t$s2[1], loo_runs = nrow(t))][] }))
base <- rbindlist(lapply(tk[sample == "none", task], function(id) rd(id)[, .(Species1, Species2, baseline = p, s1 = tk[task == id, s1], s2 = tk[task == id, s2])]))
worst <- merge(worst, base, by = c("s1", "s2", "Species1", "Species2"), all.x = TRUE)
fwrite(worst, file.path(OUT, "loo_pairs_worst.tsv.gz"), sep = "\t")

## members, cliques, orthogroups
core <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
w <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS"))[OrthoGroup %in% core$OrthoGroup]
w[, ord := .I]                                                     # input order: MeanCliqueSum, then clique number
pairs <- rbindlist(lapply(seq_len(length(G) - 1), function(i) rbindlist(lapply((i + 1):length(G), function(j)
  w[, .(cliqueID, OrthoGroup, s1 = G[i], s2 = G[j], Species1 = get(G[i]), Species2 = get(G[j]))]))))
pairs <- merge(pairs, worst[, .(s1, s2, Species1, Species2, worst, baseline)], by = c("s1", "s2", "Species1", "Species2"), all.x = TRUE)
if (pairs[is.na(baseline), .N] > 0) stop(pairs[is.na(baseline), .N], " clique pairs are missing from the unmodified runs")
pairs[is.na(worst), worst := 1][, robust := worst < P_THR]
mem <- rbindlist(lapply(G, function(g) pairs[s1 == g | s2 == g, .(genotype = g, robust_comparisons = sum(robust),
  failed_with = paste(ifelse(s1 == g, s2, s1)[!robust], collapse = ";")), by = .(cliqueID, OrthoGroup)]))
mem[, passes := robust_comparisons >= MIN_ROBUST]
fwrite(mem, file.path(OUT, "loo_members.tsv.gz"), sep = "\t")
cl <- mem[, .(passes = all(passes), failing_members = paste(genotype[!passes], collapse = ";")), by = .(cliqueID, OrthoGroup)]
cl <- merge(cl, w[, .(cliqueID, MeanCliqueSum, ord)], by = "cliqueID")[order(ord)]
fwrite(cl[, !"ord"], file.path(OUT, "loo_cliques.tsv.gz"), sep = "\t")
og <- cl[, .(complete_cliques = .N, passing_cliques = sum(passes),
             representative = if (any(passes)) cliqueID[passes][1] else "",
             old_representative_passes = passes[1], old_representative_failing = failing_members[1]), by = OrthoGroup]
og[, core_status := ifelse(passing_cliques > 0, "core", "sample-sensitive")]
og[, representative_changed := core_status == "core" & representative != core$cliqueID[match(OrthoGroup, core$OrthoGroup)]]
fwrite(og, file.path(OUT, "loo_orthogroups.tsv"), sep = "\t")
sm <- data.table(quantity = c("core_orthogroups_tested", "stay_core", "sample_sensitive", "representative_changed", "cliques_tested", "cliques_passing"),
                 value = c(nrow(og), og[core_status == "core", .N], og[core_status == "sample-sensitive", .N], og[representative_changed == TRUE, .N], nrow(cl), cl[passes == TRUE, .N]))
fwrite(sm, file.path(OUT, "loo_summary.tsv"), sep = "\t"); print(sm, row.names = FALSE)
print(og[core_status == "sample-sensitive", .(OrthoGroup, complete_cliques, old_representative_failing)], row.names = FALSE)
