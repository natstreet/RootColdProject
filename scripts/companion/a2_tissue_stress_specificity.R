#!/usr/bin/env Rscript
# a2_tissue_stress_specificity.R — A2 (ledger N159, N160, R045): is the cell-division-down / translation-maintained pattern of the
# conserved core root- and cold-specific? Written by SCI; committed by CODE with the paths as arguments and outputs in results/.
# Companion data: van Zalen et al. (AbioticStressConifers data/dds, deposited with that study). Criteria fixed in N159 before the run:
# per series, the core representative spruce or pine members of manuscript clusters 1 and 2; DESeq2 Wald log2FC (MLE, unshrunk) of each
# level vs its control; per gene the mean over the series' levels; T1 one-sided one-sample Wilcoxon cluster 1 < 0; T2 one-sided rank-sum
# cluster 2 > cluster 1; BH over the 12 series within each test; pattern present if both q < 0.05; positive control spruce root 5 C.
# Usage (working tree): Rscript scripts/companion/a2_tissue_stress_specificity.R <companion dds dir>
# Outputs: results/a2_series_tests.tsv, results/a2_per_level_medians.tsv
suppressPackageStartupMessages({ library(DESeq2); library(data.table); library(BiocParallel) })
register(SerialParam())
A  <- commandArgs(trailingOnly = TRUE)[1]; T <- "."
ld1 <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
rep <- as.data.table(readRDS(file.path(T, "data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS")))
hm  <- as.data.table(ld1(file.path(T, "data/ComPlEx/cliques/clique_HMcluster_list.rda"))); mp <- fread(file.path(T, "data/ComPlEx/cliques/manuscript_cluster_map.tsv"))
hm[, OrthoGroup := sub("_.*$", "", ID)]; hm[, ms := mp$manuscript_cluster[match(Cluster, mp$analysis_cluster)]]
rep <- merge(rep, hm[, .(OrthoGroup, ms)], by = "OrthoGroup"); stopifnot(nrow(rep) == 133)
series <- rbindlist(list(
  data.table(sp = "spruce", tissue = "root",   obj = "SCR", var = "Condition", ctrl = "20C_0h", trt = "5C",   lev = list(paste0("5C_", c("6h","24h","3d","10d")))),
  data.table(sp = "spruce", tissue = "root",   obj = "SCR", var = "Condition", ctrl = "20C_0h", trt = "-5C",  lev = list(paste0("neg5C_", c("6h","24h","3d","10d")))),
  data.table(sp = "spruce", tissue = "needle", obj = "SCN", var = "Condition", ctrl = "20C_0h", trt = "5C",   lev = list(paste0("5C_", c("6h","24h","3d","10d")))),
  data.table(sp = "spruce", tissue = "needle", obj = "SCN", var = "Condition", ctrl = "20C_0h", trt = "-5C",  lev = list(paste0("neg5C_", c("6h","24h","3d","10d")))),
  data.table(sp = "spruce", tissue = "root",   obj = "SDR", var = "Level", ctrl = "80%", trt = "drought", lev = list(c("60%","40%","30%","30%7d","Collapse","C2d"))),
  data.table(sp = "spruce", tissue = "needle", obj = "SDN", var = "Level", ctrl = "80%", trt = "drought", lev = list(c("60%","40%","30%","30%7d","Collapse","C2d"))),
  data.table(sp = "pine",   tissue = "root",   obj = "PCR", var = "Condition", ctrl = "20C_0h", trt = "5C",   lev = list(paste0("5C_", c("6h","24h","3d","10d")))),
  data.table(sp = "pine",   tissue = "root",   obj = "PCR", var = "Condition", ctrl = "20C_0h", trt = "-5C",  lev = list(paste0("neg5C_", c("6h","24h","3d","10d")))),
  data.table(sp = "pine",   tissue = "needle", obj = "PCN", var = "Condition", ctrl = "20C_0h", trt = "5C",   lev = list(paste0("5C_", c("6h","24h","3d","10d")))),
  data.table(sp = "pine",   tissue = "needle", obj = "PCN", var = "Condition", ctrl = "20C_0h", trt = "-5C",  lev = list(paste0("neg5C_", c("6h","24h","3d","10d")))),
  data.table(sp = "pine",   tissue = "root",   obj = "PDR", var = "Condition", ctrl = "FC80", trt = "drought", lev = list(c("FC60","FC40","FC30","FC30d7","Collapsed","Collapsed2d"))),
  data.table(sp = "pine",   tissue = "needle", obj = "PDN", var = "Condition", ctrl = "FC80", trt = "drought", lev = list(c("FC60","FC40","FC30","FC30d7","Collapsed","Collapsed2d")))))
fits <- list(); per_level <- list(); out <- list()
for (i in seq_len(nrow(series))) {
  s <- series[i]; key <- s$obj
  if (is.null(fits[[key]])) { d <- ld1(file.path(A, sprintf("dds_%s.rda", key))); fits[[key]] <- DESeq(d, quiet = TRUE) }
  d <- fits[[key]]; genes <- list(c1 = intersect(rep[ms == 1][[s$sp]], rownames(d)), c2 = intersect(rep[ms == 2][[s$sp]], rownames(d)))
  L <- sapply(s$lev[[1]], function(l) results(d, contrast = c(s$var, l, s$ctrl))$log2FoldChange); rownames(L) <- rownames(d)
  for (l in colnames(L)) per_level[[length(per_level) + 1]] <- data.table(sp = s$sp, tissue = s$tissue, trt = s$trt, level = l,
      c1_median = median(L[genes$c1, l], na.rm = TRUE), c2_median = median(L[genes$c2, l], na.rm = TRUE))
  m1 <- rowMeans(L[genes$c1, , drop = FALSE], na.rm = TRUE); m2 <- rowMeans(L[genes$c2, , drop = FALSE], na.rm = TRUE)
  out[[i]] <- data.table(sp = s$sp, tissue = s$tissue, trt = s$trt, n_c1 = sum(!is.na(m1)), n_c2 = sum(!is.na(m2)),
    c1_median = median(m1, na.rm = TRUE), c2_median = median(m2, na.rm = TRUE),
    T1_P = wilcox.test(m1, mu = 0, alternative = "less")$p.value, T2_P = wilcox.test(m2, m1, alternative = "greater")$p.value)
  message(s$sp, " ", s$tissue, " ", s$trt, " done")
}
res <- rbindlist(out); res[, `:=`(T1_q = p.adjust(T1_P, "BH"), T2_q = p.adjust(T2_P, "BH"))]; res[, pattern := T1_q < 0.05 & T2_q < 0.05]
dir.create("results", showWarnings = FALSE)
fwrite(res, "results/a2_series_tests.tsv", sep = "\t"); fwrite(rbindlist(per_level), "results/a2_per_level_medians.tsv", sep = "\t")
ctl <- res[sp == "spruce" & tissue == "root" & trt == "5C"]; cat("POSITIVE CONTROL (spruce root 5C) pattern:", ctl$pattern, "\n")
print(res[, .(sp, tissue, trt, n_c1, n_c2, c1 = round(c1_median, 2), c2 = round(c2_median, 2), T1_q = signif(T1_q, 2), T2_q = signif(T2_q, 2), pattern)])
