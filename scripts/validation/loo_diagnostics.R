#!/usr/bin/env Rscript
# loo_diagnostics.R — what drives failures in the sample-stability test (D052, D059; scorer loo_stability_score.R).
# For every clique pair that is below 0.1 in the unmodified run but not in the worst leave-one-out run (a "lost" pair),
# each leave-one-out run that loses it is attributed to its left-out sample, and the way it is lost is recorded:
#   absent  the pair is not in the run (ComPlEx did not test it, e.g. a gene fell below the expression filters)
#   p_rise  the pair was tested with Max.p.val >= 0.1
# Outputs (<out dir>): loo_lost_by_sample.tsv (comparison, left-out genotype and sample, lost pairs, of which absent),
#   loo_lost_summary.tsv (per comparison: clique pairs, lost pairs, lost only through one sample, lost through absence only),
#   loo_baseline_p.tsv (per comparison: unmodified-run Max.p.val quantiles of the pairs below 0.1, lost vs kept)
# Usage (working tree): Rscript loo_diagnostics.R <run dir> <scorer out dir> <out dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); RUN <- a[1]; SC <- a[2]; OUT <- a[3]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
tk <- fread(file.path(RUN, "tasks.tsv"))
rd <- function(id) fread(cmd = paste("gzip -dc", shQuote(file.path(RUN, "out", sprintf("task_%03d.tsv.gz", id)))),
                         colClasses = c(Species1 = "character", Species2 = "character", Max.p.val = "character"))[, .(Species1, Species2, p = as.numeric(Max.p.val))]
w <- fread(cmd = paste("gzip -dc", shQuote(file.path(SC, "loo_pairs_worst.tsv.gz"))))
mem <- fread(cmd = paste("gzip -dc", shQuote(file.path(SC, "loo_members.tsv.gz"))))
lost <- w[baseline < 0.1 & worst >= 0.1]
res <- rbindlist(lapply(split(tk[sample != "none"], by = c("s1", "s2")), function(t) {
  L <- lost[s1 == t$s1[1] & s2 == t$s2[1]]
  rbindlist(lapply(seq_len(nrow(t)), function(i) {
    r <- rd(t$task[i]); m <- merge(L[, .(Species1, Species2)], r, by = c("Species1", "Species2"), all.x = TRUE)
    m[, .(comparison = paste(t$s1[1], t$s2[1], sep = "-"), left_out_genotype = t$left_out_genotype[i], sample = t$sample[i],
          Species1, Species2, how = fifelse(is.na(p), "absent", fifelse(p >= 0.1, "p_rise", "kept")))][how != "kept"] })) }))
by_s <- res[, .(lost_pairs = .N, absent = sum(how == "absent")), by = .(comparison, left_out_genotype, sample)][order(-lost_pairs)]
fwrite(by_s, file.path(OUT, "loo_lost_by_sample.tsv"), sep = "\t")
per_pair <- res[, .(n_samples = uniqueN(paste(left_out_genotype, sample)), all_absent = all(how == "absent")), by = .(comparison, Species1, Species2)]
sm <- merge(w[, .(clique_pairs = .N, baseline_below = sum(baseline < 0.1, na.rm = TRUE)), by = .(comparison = paste(s1, s2, sep = "-"))],
            per_pair[, .(lost_pairs = .N, lost_by_one_sample = sum(n_samples == 1), lost_only_by_absence = sum(all_absent)), by = comparison], by = "comparison", all.x = TRUE)
fwrite(sm, file.path(OUT, "loo_lost_summary.tsv"), sep = "\t")
q <- w[baseline < 0.1, .(pairs = .N, q25 = quantile(baseline, 0.25), median = median(baseline), q75 = quantile(baseline, 0.75),
                          share_above_0.01 = mean(baseline >= 0.01)), by = .(comparison = paste(s1, s2, sep = "-"), status = fifelse(worst >= 0.1, "lost", "kept"))][order(comparison, status)]
fwrite(q, file.path(OUT, "loo_baseline_p.tsv"), sep = "\t")
print(sm, row.names = FALSE); print(head(by_s, 15), row.names = FALSE); print(q, row.names = FALSE)
