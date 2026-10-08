#!/usr/bin/env Rscript
# stability_timepoint_attribution.R — for the 15 pairs of each representative clique (kept and sample-sensitive core
# orthogroups; stability_kept_lost.R), which left-out samples make a pair that is below 0.1 in the unmodified run reach
# 0.1 or drop out, grouped by the left-out sample's time point (01 control, 02 6 h, 03 24 h, 04 3 d, 05 10 d) and genotype
# (ledger R041: are transient, time-point-specific responses the fragile ones?).
# Output (<interpretation dir>): timepoint_attribution.tsv (status, cluster, left-out genotype, time point, pair losses,
#   leave-one-out runs of that time point, losses per run)
# Usage (working tree): Rscript stability_timepoint_attribution.R <run dir> <interpretation dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); RUN <- a[1]; D <- a[2]
G <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
kl <- fread(file.path(D, "kept_lost_orthogroups.tsv")); w <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS"))
pairs <- rbindlist(lapply(seq_len(length(G) - 1), function(i) rbindlist(lapply((i + 1):length(G), function(j)
  w[cliqueID %in% kl$clique, .(cliqueID, s1 = G[i], s2 = G[j], Species1 = get(G[i]), Species2 = get(G[j]))]))))
pairs <- merge(pairs, kl[, .(cliqueID = clique, core_status, cluster)], by = "cliqueID")
tk <- fread(file.path(RUN, "tasks.tsv"))
rd <- function(id) fread(cmd = paste("gzip -dc", shQuote(file.path(RUN, "out", sprintf("task_%03d.tsv.gz", id)))),
                         colClasses = c(Species1 = "character", Species2 = "character", Max.p.val = "character"))[, .(Species1, Species2, p = as.numeric(Max.p.val))]
base <- rbindlist(lapply(tk[sample == "none", task], function(id) rd(id)[, `:=`(s1 = tk[task == id, s1], s2 = tk[task == id, s2])]))
pairs <- merge(pairs, base[, .(s1, s2, Species1, Species2, base_p = p)], by = c("s1", "s2", "Species1", "Species2"), all.x = TRUE)[base_p < 0.1]
ev <- rbindlist(lapply(tk[sample != "none", task], function(id) {
  t <- tk[task == id]; P <- pairs[s1 == t$s1 & s2 == t$s2]; if (!nrow(P)) return(NULL)
  m <- merge(P, rd(id), by = c("Species1", "Species2"), all.x = TRUE)
  m[, .(core_status, cluster, left_out_genotype = t$left_out_genotype, time_point = sub(".*\\.", "", t$sample), lost = is.na(p) | p >= 0.1)] }))
runs <- unique(tk[sample != "none", .(left_out_genotype, time_point = sub(".*\\.", "", sample), sample)])[, .(loo_samples = .N), by = .(left_out_genotype, time_point)]
res <- ev[, .(pair_tests = .N, losses = sum(lost)), by = .(core_status, cluster, left_out_genotype, time_point)]
res <- merge(res, runs, by = c("left_out_genotype", "time_point"))[, loss_rate := losses / pair_tests]
setorder(res, core_status, cluster, left_out_genotype, time_point)
fwrite(res, file.path(D, "timepoint_attribution.tsv"), sep = "\t")
print(dcast(res[core_status == "sample-sensitive", .(losses = sum(losses), tests = sum(pair_tests)), by = .(cluster, time_point)][, rate := round(losses / tests, 3)],
            cluster ~ time_point, value.var = "rate"))
print(dcast(res[core_status == "core", .(losses = sum(losses), tests = sum(pair_tests)), by = .(cluster, time_point)][, rate := round(losses / tests, 3)],
            cluster ~ time_point, value.var = "rate"))
print(dcast(res[, .(losses = sum(losses), tests = sum(pair_tests)), by = .(core_status, left_out_genotype)][, rate := round(losses / tests, 3)],
            left_out_genotype ~ core_status, value.var = "rate"))
