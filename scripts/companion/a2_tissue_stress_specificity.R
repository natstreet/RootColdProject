#!/usr/bin/env Rscript
# a2_tissue_stress_specificity.R — A2 (ledger N159, N160, R045): is the cell-division-down / translation-maintained pattern of the
# conserved core root- and cold-specific? Written by SCI; committed by CODE with the paths as arguments and outputs in results/.
# Companion data: van Zalen et al. (AbioticStressConifers data/dds, deposited with that study). Criteria fixed in N159 before the run:
# per series, the core representative spruce or pine members of the cell-division and translation clusters (selected by role, D095); DESeq2 Wald log2FC (MLE, unshrunk) of each
# level vs its control; per gene the mean over the series' levels; T1 one-sided one-sample Wilcoxon cell-division cluster < 0; T2 one-sided rank-sum
# translation cluster > cell-division cluster; BH over the 12 series within each test; pattern present if both q < 0.05; positive control spruce root 5 C.
# Usage (working tree): Rscript scripts/companion/a2_tissue_stress_specificity.R <companion dds dir> [--include-unexpressed|--default] [out dir]   (default excludes unexpressed genes, D096)
# Outputs: results/a2_series_tests.tsv, results/a2_per_level_medians.tsv
suppressPackageStartupMessages({ library(DESeq2); library(data.table); library(BiocParallel) })
register(SerialParam())
ARGS <- commandArgs(trailingOnly = TRUE); A <- ARGS[1]; T <- "."
EXCL <- !(length(ARGS) >= 2 && ARGS[2] == "--include-unexpressed")  # D096: a gene with zero counts in both groups of a contrast gets NA
OUTD <- if (length(ARGS) >= 3) ARGS[3] else "results"               #   instead of DESeq2's log2FC of exactly 0 (--include-unexpressed: as N159)
ld1 <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
rep <- as.data.table(readRDS(file.path(T, "data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS")))
hm  <- as.data.table(ld1(file.path(T, "data/ComPlEx/cliques/clique_HMcluster_list.rda"))); mp <- fread(file.path(T, "data/ComPlEx/cliques/manuscript_cluster_map.tsv"))
hm[, OrthoGroup := sub("_.*$", "", ID)]; hm[, ms := mp$manuscript_cluster[match(Cluster, mp$analysis_cluster)]]
n_core <- nrow(rep); rep <- merge(rep, hm[, .(OrthoGroup, ms)], by = "OrthoGroup"); stopifnot(nrow(rep) == n_core)   # every core orthogroup has a cluster
source("scripts/ComPlEx/cluster_roles.R")                                       # clusters by role (D095); output columns carry their numbers
K <- c(cd = cluster_of("cell_division"), tr = cluster_of("translation"), ep = cluster_of("epidermal_root_hair"))
cn <- function(fmt, r) sprintf(fmt, K[[r]])
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
  d <- fits[[key]]; genes <- lapply(K, function(k) intersect(rep[ms == k][[s$sp]], rownames(d)))      # cd, tr, ep (ep: D088, D093)
  L <- sapply(s$lev[[1]], function(l) results(d, contrast = c(s$var, l, s$ctrl))$log2FoldChange); rownames(L) <- rownames(d)
  if (EXCL) { cn0 <- counts(d); grp <- colData(d)[[s$var]]
    for (l in colnames(L)) L[rowSums(cn0[, grp %in% c(l, s$ctrl), drop = FALSE]) == 0, l] <- NA }
  for (l in colnames(L)) per_level[[length(per_level) + 1]] <- as.data.table(c(list(sp = s$sp, tissue = s$tissue, trt = s$trt, level = l),
      setNames(lapply(names(K), function(r) median(L[genes[[r]], l], na.rm = TRUE)), sapply(names(K), function(r) cn("c%d_median", r)))))
  m <- lapply(genes, function(g) rowMeans(L[g, , drop = FALSE], na.rm = TRUE))
  out[[i]] <- as.data.table(c(list(sp = s$sp, tissue = s$tissue, trt = s$trt),
    setNames(list(sum(!is.na(m$cd)), sum(!is.na(m$tr)), median(m$cd, na.rm = TRUE), median(m$tr, na.rm = TRUE)),
             c(cn("n_c%d", "cd"), cn("n_c%d", "tr"), cn("c%d_median", "cd"), cn("c%d_median", "tr"))),
    list(T1_P = wilcox.test(m$cd, mu = 0, alternative = "less")$p.value,                 # T1: cell-division cluster below 0
         T2_P = wilcox.test(m$tr, m$cd, alternative = "greater")$p.value),               # T2: translation cluster above cell division
    setNames(list(sum(!is.na(m$ep)), median(m$ep, na.rm = TRUE), wilcox.test(m$ep, mu = 0, alternative = "less")$p.value,
                  wilcox.test(m$tr, m$ep, alternative = "greater")$p.value),                       # D088: the same tests for the epidermal cluster
             c(cn("n_c%d", "ep"), cn("c%d_median", "ep"), cn("T1c%d_P", "ep"), cn("T2c%d_P", "ep")))))
  message(s$sp, " ", s$tissue, " ", s$trt, " done")
}
res <- rbindlist(out); res[, `:=`(T1_q = p.adjust(T1_P, "BH"), T2_q = p.adjust(T2_P, "BH"))]; res[, pattern := T1_q < 0.05 & T2_q < 0.05]
res[, (cn("T1c%d_q", "ep")) := p.adjust(get(cn("T1c%d_P", "ep")), "BH")][, (cn("T2c%d_q", "ep")) := p.adjust(get(cn("T2c%d_P", "ep")), "BH")]
res[, (cn("pattern_c%d", "ep")) := get(cn("T1c%d_q", "ep")) < 0.05 & get(cn("T2c%d_q", "ep")) < 0.05]
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
fwrite(res, file.path(OUTD, "a2_series_tests.tsv"), sep = "\t"); fwrite(rbindlist(per_level), file.path(OUTD, "a2_per_level_medians.tsv"), sep = "\t")
ctl <- res[sp == "spruce" & tissue == "root" & trt == "5C"]; cat("POSITIVE CONTROL (spruce root 5C) pattern:", ctl$pattern, "\n")
print(res[, c("sp", "tissue", "trt", cn("c%d_median", "cd"), cn("c%d_median", "tr"), "T1_q", "T2_q", "pattern"), with = FALSE])
